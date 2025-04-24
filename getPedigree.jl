mutable struct PedNode
    seqID::Int64
    sire::String
    dam::String
    f::Float64
end

"""
    getPedigree(pedfile::Union{AbstractString,DataFrames.DataFrame}; header=false, separator=',', missingstrings=["0"])

Reads and processes a pedigree file or DataFrame to construct a detailed pedigree for individuals, where each individual is represented by a `PedNode` structure containing information on the individual's unique sequence ID, sire ID, dam ID, and inbreeding coefficient. 

The function performs the following steps:
1. If a file path is provided (`pedfile` is a string), it reads the file as a CSV file into a DataFrame. If a DataFrame is already provided, it is used directly.
2. Strips any leading or trailing spaces from the columns of the DataFrame.
3. Initializes a mapping of individual IDs to `PedNode` objects, with missing data filled as "missing" for sire and dam, and the inbreeding coefficient set to -1.
4. Identifies and fills missing pedigree information where applicable, such as missing sire or dam IDs.
5. Constructs the ancestry information for each individual by tracing back through the sire and dam lines, ensuring no cycles or incorrect ancestry links.
6. If any inconsistencies are found in an individual's pedigree (e.g., an individual appears as its own ancestor), the function prints a warning and sets the problematic ancestor to "missing".
7. Orders the individuals by their ancestry, ensuring that individuals are listed only once all of their ancestors have been accounted for.
8. Creates a new DataFrame with columns for sequence IDs, sire IDs, dam IDs, and original individual IDs.

### Arguments:
- `pedfile::Union{AbstractString, DataFrames.DataFrame}`: The pedigree file (as a path to a CSV file) or DataFrame containing pedigree data.
- `header::Bool`: Flag indicating whether the file contains a header row (default is `false`).
- `separator::Char`: The delimiter used in the CSV file (default is `','`).
- `missingstrings::Vector{String}`: List of strings that represent missing values in the dataset (default is `["0"]`).

### Returns:
A DataFrame containing:
- `seq_IDs`: Sequence IDs assigned to individuals.
- `sire_IDs`: Sequence IDs of the sires.
- `dam_IDs`: Sequence IDs of the dams.
- `ori_IDs`: Original individual IDs from the pedigree input.
"""
function getPedigree(pedfile::Union{AbstractString,DataFrames.DataFrame};header=false,separator=',',missingstrings=["0"])
    if typeof(pedfile) <: AbstractString
        printstyled("The delimiter in ",split(pedfile,['/','\\'])[end]," is \'",separator,"\'.\n",bold=false,color=:green)
        df  = CSV.read(pedfile,DataFrame,types=[String,String,String,String,String],
                        delim=separator,header=header,missingstring=missingstrings)
    elseif typeof(pedfile) == DataFrames.DataFrame
        df  = pedfile
    else
        error("Please provide a file path or a dataframe.")
    end
    df[!,1]=strip.(string.(df[!,1]))
    df[!,2]=strip.(string.(df[!,2]))
    df[!,3]=strip.(string.(df[!,3]))
    df[!,4]=strip.(string.(df[!,4]))

    idMap = Dict{AbstractString,PedNode}()
    
    # 填充缺失的系谱信息
    n = size(df,1)
    for i in unique(vcat(df[:,1], df[:, 2], df[:, 3]))
        idMap[i]=PedNode(0,"missing","missing", -1.0)
    end
    for i in 1:size(df, 1)
        idMap[df[i, 1]]=PedNode(0,df[i,2],df[i,3], -1.0)
    end
    if "missing" in keys(idMap)
        delete!(idMap, "missing")
    end

    # 错误校对
    anc = Dict()
    for i in keys(idMap)
        anc[i] = Dict()
        if idMap[i].sire != "missing"
            anc[i][idMap[i].sire] = 1
        end
        if idMap[i].dam != "missing"
            anc[i][idMap[i].dam] = 1
        end
        d = 1000
        while d > 0
            id_key = keys(anc[i])
            id_key_len1 = length(id_key)
            for j in id_key
                if j in keys(idMap)
                    anc[i][idMap[j].sire] = 0
                    anc[i][idMap[j].dam] = 0
                end
            end
            if "missing" in keys(anc[i])
                delete!(anc[i], "missing")
            end
            id_key = keys(anc[i])
            id_key_len2 = length(id_key)
            d = id_key_len2 - id_key_len1
            if d == 0
                break
            end
        end
    end

    for i in keys(anc)
        if i in keys(anc[i])
            println("The pedigree of individual $(i) is incorrect.")
            if idMap[i].sire in keys(anc[i])
                idMap[i].sire = "missing"
                println("Change the sire ID of individual $(i) to missing.")
            end
            if idMap[i].dam in keys(anc[i])
                idMap[i].dam =  "missing"
                println("Change the dam ID of individual $(i) to missing.")
            end
        end
    end

    # 系谱重新排序
    sort_idMap = Dict{AbstractString,PedNode}()
    output = Dict()
    output["missing"] = 1
    seqID = 0
    ped_len = length(keys(idMap))
    while ped_len > 0
        for ikey in keys(idMap)
            if idMap[ikey].sire in keys(output) && idMap[ikey].dam in keys(output)
                seqID += 1
                sort_idMap[ikey] = PedNode(seqID, idMap[ikey].sire, idMap[ikey].dam, -1.0)
                output[ikey] = 1
                delete!(idMap, ikey)
            end
        end
        ped_len = length(keys(idMap))
    end

    n = length(sort_idMap)
    seq_IDs = 1:n
    sires = Array{Int}(undef, n)
    dams = Array{Int}(undef, n)
    ori_IDs = Array{String}(undef, n)
    for i in keys(sort_idMap)
        if sort_idMap[i].sire == "missing"
            sires[sort_idMap[i].seqID] = 0
        else
            sires[sort_idMap[i].seqID] = sort_idMap[sort_idMap[i].sire].seqID
        end
        if sort_idMap[i].dam == "missing"
            dams[sort_idMap[i].seqID] = 0
        else
            dams[sort_idMap[i].seqID] = sort_idMap[sort_idMap[i].dam].seqID
        end
        ori_IDs[sort_idMap[i].seqID] = i
    end
    
    ped = DataFrame(;seq_IDs = seq_IDs, sire_IDs = sires, dam_IDs = dams, ori_IDs = ori_IDs)

    return ped
end
