"""
Geno

A struct that represents genotype data from a genetic study. It contains:

- `geno::Matrix` - The genotype data
- `fam::DataFrame` - Family/phenotype data
- `map::DataFrame` - SNP annotation data
- `n` - Number of individuals
- `m` - Number of SNPs

The Geno struct is initialized with the `Geno` constructor.
This struct and associated functions provide a way to represent and process genotype data from genetic studies.
"""
mutable struct Geno
    geno::Matrix
    fam::DataFrame
    map::DataFrame
    n
    m
    function Geno(geno, fam, map)
        new(geno, fam, map, size(geno, 1), size(geno, 2))
    end
end


"""
    readPlink(file::String) -> Geno

Read PLINK bed/bim/fam files and return a Geno object.

# Arguments
- `file::String`: The base filename (without extensions) of the PLINK files.

# Returns
A `Geno` object containing the following fields:
- `geno`: The genotype data.
- `fam`: The individual information.
- `bim`: The SNP information.

# Example
```julia
geno = readPlink("data/my_plink_data")
```
"""
function readPlink(file)
    bed_file = string(file, ".bed")
    bim_file = string(file, ".bim")
    fam_file = string(file, ".fam")
    bim = CSV.read(bim_file, DataFrame; header = false)
    rename!(bim, [:CHROM, :SNP, :Genetic_Distance, :POS, :A1, :A2])
    bim = bim[:, [1,2,4,5,6]]
    fam = CSV.read(fam_file, DataFrame; header = false, missingstring="NA")
    rename!(fam[:, 1:6], [:Family_ID, :Individual_ID, :Paternal_ID, :Maternal_ID, :Sex, :Trait])
    # fam = DataFrame(;:ID => fam[:, 2])
    bed = read_bed(bed_file, size(fam, 1), size(bim, 1))
    return Geno(bed, fam, bim)
end

function read_bed(genfil, nind, nloci)

    code = Dict(0 => 2, 2 => 1, 1 => 3, 3 => 0) # an array for mapping genotype values

    # open a binary file, read plink magic numbers and mode, check if it is SNP-major mode
    f = open(genfil, "r")
    b1 = read(f, UInt8) # plink magic number 1
    if b1 != parse(UInt8,"01101100",base =2)
        println("Binary genotype file may not be a plink file")
        return
    end

    b1 = read(f, UInt8) # plink magic number 2
    if b1 != parse(UInt8,"00011011",base =2)
        println("Binary genotype file may not be a plink file")
        return
    end

    b1 = read(f, UInt8) # mode
    if b1 != 0x01
        println("SNP file should not be in individual mode")
        return
    end

    len = filesize(f) - 3
    buffer = Vector{Int8}(undef, len*4)

    for i in 1:len
        p = read(f, UInt8)
        # r = i ÷ nind
        for x in 0:3
           buffer[4 * (i - 1) + x + 1] = code[(p >> (2 * x))&0x03]
        end
    end

    close(f)
    file = "$(chop(genfil, tail = 4)).tmp.geno.bin"
    println("Mapping genotype data to $(file).")
    s = open(file, "w+")
    write(s, reshape(buffer, Int(len*4/nloci), nloci)[1:nind, :])
    close(s)
    s = open(file, "r+")
    geno = mmap(s, Matrix{Int8}, (nind, nloci))
    close(s)
    # geno = reshape(buffer, Int(len*4/nloci), nloci)[1:nind, :]
    # return the geno matrix
    return geno
end


"""
    solveMME(X, Z, Kinv, λ, y) -> (b, u)

Solves the mixed model equations (MME) for a linear mixed model of the form:
    y = Xb + Zu + e,
where `b` is the vector of fixed effects and `u` is the vector of random effects with 
    u ~ N(0, K * σ²_u), and e ~ N(0, I * σ²_e). 
The MME incorporates a penalty term λ = σ²_e / σ²_u.

# Arguments
- `X`: Design matrix for fixed effects.
- `Z`: Design matrix for random effects.
- `Kinv`: Inverse of the covariance matrix for random effects (typically the inverse of a genomic relationship matrix).
- `λ`: Regularization parameter (σ²_e / σ²_u).
- `y`: Vector of phenotypic values.

# Returns
- `(b, u)`: Tuple containing the estimated fixed effects (`b`) and random effects (`u`).

# Notes
- Constructs and solves the mixed model equations using conjugate gradient iteration.
- Sparse matrix techniques are used for efficiency, especially in large-scale genomic datasets.
"""
function solveMME(X, Z, Kinv, λ, y)
    W = sparse(hcat(X, Z))
    LHS = W'W
    LHS[size(X, 2)+1:end, size(X, 2)+1:end] = LHS[size(X, 2)+1:end, size(X, 2)+1:end] + λ .* Kinv
    LHS = Matrix(LHS)
    RHS = W'y
    x = ones(size(LHS, 1))
    IterativeSolvers.cg!(x, LHS, RHS)
    return x[1:size(X, 2)], x[size(X, 2)+1:end]
end



"""
    pQsel(y, Kinv, geno, meta, λ; n_pre=100, r2_threshold=0.99, n_fold=10) -> Vector{Int}

Selects pseudo quantitative trait nucleotides (pQTNs) based on GWAS p-values and cross-validated prediction accuracy.

# Arguments
- `y`: A vector of phenotypic values.
- `Kinv`: The inverse of the relationship matrix (e.g., inverse of the genomic relationship matrix).
- `geno`: A structure Geno.
- `meta`: A DataFrame-like structure containing GWAS results, with a column `p_values` holding the p-values for each SNP.
- `λ`: Regularization parameter used in solving the mixed model equations.
- `n_pre`: (Optional) Number of top candidate SNPs to consider for selection (default is 100).
- `r2_threshold`: (Optional) Maximum allowed squared correlation (LD) between selected SNPs to ensure low collinearity (default is 0.99).
- `n_fold`: (Optional) Number of folds used in cross-validation (default is 10).

# Returns
- A vector of indices corresponding to selected SNPs (pQTNs) that improve prediction accuracy when added to the model.

# Description
The function performs a stepwise selection of SNPs from the top `n_pre` most significant SNPs ranked by GWAS p-values, ensuring that added SNPs are not in high linkage disequilibrium (LD) with previously selected ones.

A nested cross-validation loop is used to assess the contribution of each candidate SNP to the prediction accuracy of a linear mixed model. A SNP is retained in the final pQTN set if it improves prediction accuracy in at least 90% of cross-validation folds by a margin of at least 0.001.

This method allows the construction of a set of informative SNPs (pQTNs) that can be used as fixed effects to enhance genomic prediction accuracy, particularly in models such as ssGWABLUP_pQTNs.
"""

function pQsel(y, Kinv, geno, meta, λ; n_pre=100, r2_threshold = 0.99, n_fold = 10)
    index = [sortperm(meta.p_values)[1]]
    i = 1
    while length(index) < n_pre
        i += 1
        if all(cor(geno.geno[:,sortperm(meta.p_values)[i]], geno.geno[:,index]) .^ 2 .< r2_threshold)
            push!(index, sortperm(meta.p_values)[i])
        end
    end

    function cv(y, Kinv, λ, X, Z)
        fold_size = round(Int, length(y)/n_fold)
        accs = []
        for fold in 1:n_fold
            vid = (fold-1)*fold_size+1 : min(fold*fold_size, length(y))
            id = setdiff(1:length(y), vid)
            b, u = solveMME(X[id, :], Z[id, :], Kinv, λ, y[id, :])
            yhat = X * b + u
            push!(accs, cor(yhat[vid], y[vid]))
        end
        return accs
    end
    X = ones(length(y), 1)
    Z = Matrix(I, length(y), length(y))
    pQs = []

    old_accs = cv(y, Kinv, λ, hcat(X, geno.geno[:, pQs]), Z)
    for pQ in index
        new_accs = cv(y, Kinv, λ, hcat(X, geno.geno[:, vcat(pQs, pQ)]), Z)
        if sum(new_accs .- old_accs .>= 0.001) >= 0.9*n_fold
            old_accs = copy(new_accs)
            push!(pQs, pQ)
        end
    end
    return pQs
end
