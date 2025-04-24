using Mmap, DataFrames, CSV, LinearAlgebra, Optim, Statistics, StatsBase, HypothesisTests, IterativeSolvers, SparseArrays
"""
    calA(ped::DataFrame)

Calculates a kinship matrix `A` for a given pedigree DataFrame, where each entry represents the relationship between pairs of individuals. The matrix is symmetric, with entries determined based on the sire and dam information provided in the pedigree.

The kinship matrix is calculated as follows:
1. The matrix `A` is initialized as a square matrix of zeros with dimensions equal to the number of individuals (`n`).
2. For each individual pair (i, j), the relationship is calculated:
   - If `i == j` (i.e., the individual is compared to itself), the diagonal entry is set to 1. If either the sire or dam of the individual is missing (`0`), the kinship value is set to 1, otherwise, it is calculated based on the kinship between the sire and dam of the individual.
   - For `i != j`, the relationship is calculated by averaging the kinship values between the two individuals' sires and dams, with special cases where one or both parents are missing.
3. The matrix is symmetric, meaning the value at `A[i, j]` is the same as `A[j, i]`.

### Arguments:
- `ped::DataFrame`: A DataFrame containing pedigree data with at least three columns:
  - Column 2: The sire ID for each individual (0 if missing).
  - Column 3: The dam ID for each individual (0 if missing).

### Returns:
A square matrix `A` of size `n x n`, where `n` is the number of individuals in the pedigree. Each element `A[i, j]` represents the kinship coefficient between individuals `i` and `j`.
"""
function calA(ped::DataFrame)
    n = size(ped, 1)
    A = zeros(n, n)
    for i in 1:n
        for j in i:n
            sj = ped[j, 2]
            dj = ped[j, 3]
            if i == j
                if sj == 0 || dj == 0
                    A[i, j] = 1
                else
                    A[i, j] = 1 + 0.5 * A[sj, dj]
                end
            else
                if sj != 0 && dj !=0
                    A[i,j] = 0.5*(A[i, sj]+A[i, dj])
                elseif sj == 0 && dj !=0
                    A[i,j]= 0.5*A[i, dj]
                elseif sj != 0 && dj ==0
                    A[i,j] = 0.5*A[i, sj]
                end
                A[j, i] = A[i,j]
            end
            # println("$(i) and $(j) : $(A[i,j])")
        end
    end
    return A
end


function cal_aij!(A, ped, i, j)
    # if i=="missing" || j=="missing"           # zero
    #     return 0.0
    # end
    if i==0 || j==0
        return 0.0
    end
    # old,yng = idMap[i].seqID < idMap[j].seqID ? (i,j) : (j,i)
    # oldID = idMap[old].seqID
    # yngID = idMap[yng].seqID
    oldID, yngID = i < j ? (i, j) : (j, i)
    # sireOfYng = idMap[yng].sire
    # damOfYng  = idMap[yng].dam
    sireOfYng = ped[yngID, "sire_IDs"]
    damOfYng = ped[yngID, "dam_IDs"]
    if A[yngID, oldID] > 0
        return A[yngID, oldID]
    end

    if oldID == yngID
        aij = 1 + 0.5 * cal_aij!(A, ped, sireOfYng, damOfYng)
        A[yngID, oldID] = aij
        # A[oldID, yngID] = aij
        return aij
    end
    aOldDamYoung  = (oldID == 0|| damOfYng == 0) ? 0.0 : cal_aij!(A, ped,oldID,damOfYng)
    aOldSireYoung = (oldID == 0 || sireOfYng== 0) ? 0.0 : cal_aij!(A, ped,oldID,sireOfYng)
    aij = 0.5*(aOldSireYoung + aOldDamYoung)
    A[yngID, oldID] = aij
    # A[oldID, yngID] = aij
    return aij
end

"""
    calInbreeding!(ped::DataFrame)

This function calculates the inbreeding coefficients for each individual in the pedigree dataset.
It computes the kinship matrix (A) for all individuals and updates the `f` column in the `ped` DataFrame with the inbreeding coefficients.

### Steps:
1. A sparse matrix `A` is initialized to store the pairwise kinship coefficients between all individuals.
2. The function iterates over each individual in the `ped` DataFrame and calculates the kinship coefficient between each individual and themselves using the helper function `cal_aij!`.
3. The diagonal of the matrix `A` is updated with the inbreeding coefficients (calculated as 1 - relatedness), and these values are stored in the `f` column of the `ped` DataFrame.
4. The updated `ped` DataFrame with inbreeding coefficients is returned.

### Arguments:
- `ped::DataFrame`: A DataFrame containing pedigree information. Each row represents an individual and includes columns for sire, dam, and other pedigree-related data.

### Returns:
- `ped::DataFrame`: The updated pedigree DataFrame with the `f` column containing inbreeding coefficients for each individual.

"""
function calInbreeding!(ped::DataFrame)
    A = spzeros(size(ped, 1), size(ped, 1))
    for i in 1:size(ped, 1)
        cal_aij!(A, ped, i, i)
    end
    ped.f = Vector(diag(A)) .- 1
    return ped
end

function HAi(ped::DataFrame)
    ii = Int64[]
    jj = Int64[]
    vv = Float64[]
    for ind in 1:size(ped, 1)
        sire = ped.sire_IDs[ind]
        dam  = ped.dam_IDs[ind]
        sirePos = ped.sire_IDs[ind]
        damPos  = ped.dam_IDs[ind]
        myPos   = ind
        if sirePos>0 && damPos>0
            # d = sqrt(4.0/(2 - idMap[sire].f - idMap[dam].f))
            d = sqrt(4.0/(2 - ped.f[sire] - ped.f[dam]))
            push!(ii,myPos)
            push!(jj,sirePos)
            push!(vv,-0.5*d)
            push!(ii,myPos)
            push!(jj,damPos)
            push!(vv,-0.5*d)
            push!(ii,myPos)
            push!(jj,myPos)
            push!(vv,d)
         elseif sirePos>0
            d = sqrt(4.0/(3 - ped.f[sire]))
            push!(ii,myPos)
            push!(jj,sirePos)
            push!(vv,-0.5*d)
            push!(ii,myPos)
            push!(jj,myPos)
            push!(vv,d)
          elseif damPos>0
            d = sqrt(4.0/(3 - ped.f[dam]))
            push!(ii,myPos)
            push!(jj,damPos)
            push!(vv,-0.5*d)
            push!(ii,myPos)
            push!(jj,myPos)
            push!(vv,d)
        else
            d = 1.0
            push!(ii,myPos)
            push!(jj,myPos)
            push!(vv,d)
        end
    end
    return (ii,jj,vv)
end

"""
    calAInverse(ped::DataFrame)

This function calculates the inverse of the kinship matrix for a given pedigree dataset.

### Arguments:
- `ped::DataFrame`: A DataFrame containing pedigree information for a set of individuals. This data is used to calculate the kinship matrix.

### Returns:
- `Ai`: The inverse of the kinship matrix, which is a sparse matrix containing the pairwise relatedness information of the individuals in the pedigree.

"""
function calAInverse(ped::DataFrame)
    ii,jj,vv = HAi(ped)
    hAi      = sparse(ii,jj,vv)
    Ai       = hAi'hAi
    return Ai
end

"""
    calG(geno::Geno; weights::Vector{Float64}=nothing)

This function calculates the genomic relationship matrix (GRM) from a genotype dataset. There are two versions of the function depending on whether a `weights` vector is provided.

### Steps:
1. The allele frequencies are computed from the genotypic data by taking the mean of the genotypes across individuals.
2. The genotypic data is centered and scaled by subtracting the allele frequencies and dividing by the standard deviation of the allele frequencies.
3. If weights are provided, the genotypic data is additionally scaled by the weights.
4. The relationship matrix `K` is computed by multiplying the standardized genotypic matrix `Z` with its transpose.
5. The resulting relationship matrix `K` is scaled by dividing by the number of markers `geno.m`.
6. The function returns the symmetric genomic relationship matrix `K`.

### Arguments:
- `geno::Geno`: An object containing the genotypic data, including the genotypes of individuals and the number of markers (`geno.m`).
- `weights::Vector{Float64}` (optional): A vector of weights applied to the genotypic data before calculating the relationship matrix. If not provided, the genotypic data is used without weights.

### Returns:
- `K`: The genomic relationship matrix (GRM), a symmetric matrix representing the pairwise genetic relatedness between individuals. If `weights` are provided, the matrix is weighted accordingly.

"""
function calG(geno::Geno)
    p = mean(geno.geno, dims=1)' ./ 2
    Z = Matrix((geno.geno' .- 2p)
    ZZ = Z'Z
    K = ZZ ./ sum(2 .* p .* (1 .- p))
    return Symmetric(K)
end

function calG(geno::Geno, weights::Vector{Float64})
    p = mean(geno.geno, dims=1)' ./ 2
    sw = sqrt.(weights)
    Z = Matrix((geno.geno' .- 2p) .* sw
    ZZ = Z'Z
    K = ZZ ./ sum(2 .* p .* (1 .- p))
    return Symmetric(K)
end

"""
    calH(A, A22inv, G, gid)

This function calculates a matrix \( H \) based on the provided genomic relationship matrix \( G \), the genomic matrix \( A \), the inverse of the submatrix \( A_{22} \), and the set of IDs (`gid`). It performs a transformation of the matrix structure and generates a modified matrix \( H \) which accounts for genetic information from the given individuals.

### Steps:
1. **Reordering the Individuals**: 
   The function begins by reordering the individuals. It combines the indices of those without genotype information (`setdiff(1:size(A, 1), gid)`) and those with genotype information (`gid`) into a new index array `ids`. This reordering is critical for partitioning the matrix `A` into relevant submatrices.

2. **Partitioning the Matrix \( A \)**:
   The matrix \( A \) is partitioned into several submatrices:
   - `A12`: A matrix representing the relationship between the individuals with and without genotype information.
   - `A22`: A matrix corresponding to the individuals with genotype information.

3. **Calculating the Matrix \( W \)**:
   The matrix \( W \) is constructed using the product of `A12` and `A22inv`, along with the identity matrix \( I \) of the same size as \( G \). This transformation accounts for the genetic relationships between the individuals in the study.

4. **Matrix Modification**:
   The final matrix \( H \) is calculated as:
   \[
   H = A_{sort} + W \cdot (G - A_{22}) \cdot W'
   \]
   where \( A_{sort} \) is the reordered matrix, and \( G - A_{22} \) represents the difference between the genomic relationship matrix \( G \) and the matrix \( A_{22} \).

5. **Returning the Result**:
   The final matrix \( H \) is reordered using `sortperm(ids)` to ensure the correct ordering of individuals. The result is a matrix reflecting the genetic relationships, adjusted for the individuals with known genotypes.

### Arguments:
- `A`: The matrix representing the genetic relationships between individuals.
- `A22inv`: The inverse of the submatrix \( A_{22} \) from the matrix \( A \), which represents the relationship among individuals with genotype data.
- `G`: The genomic relationship matrix representing genetic relatedness.
- `gid`: A vector containing the indices of individuals with genotype information.

### Returns:
- `H`: A matrix that reflects the adjusted genetic relationships among the individuals, considering both those with and without genotype information.

"""
function calH(A, A22inv, G, gid)
    ids = vcat(setdiff(1:size(A, 1), gid), gid)
    n_nogeno = size(A, 1) - length(gid)
    A_sort = A[ids, ids]
    A12 = A_sort[1:n_nogeno, n_nogeno+1:end]
    A22 = A_sort[n_nogeno+1:end, n_nogeno+1:end]
    GminusA22 = G - A22
    W = [A12 * A22inv
        Matrix(I, size(G))]

    H = A_sort + W * GminusA22 * W'
    return H[sortperm(ids), sortperm(ids)]
end

"""
    calHInverse(Ginv, Ainv, A22inv, gid)

This function computes the inverse of matrix \( H \), denoted as \( H^{-1} \), using the inverse of the genomic relationship matrix \( G \), the inverse of the matrix \( A \), and the inverse of the submatrix \( A_{22} \). The function updates the matrix \( H^{-1} \) specifically for the individuals with genotype data.

### Steps:
1. **Copying the Inverse Matrix**:
   The function begins by copying the inverse matrix \( A^{-1} \) to create a matrix `Hi`.

2. **Updating the Submatrix for Genotyped Individuals**:
   The function then modifies the submatrix corresponding to the individuals with genotype information. Specifically, for the individuals in the `gid` set, it updates their values using the formula:
   \[
   H^{-1}_{gid, gid} = G^{-1} - A_{22}^{-1}
   \]
   This adjustment reflects the inverse relationship between the genotyped individuals' genetic data and their phenotypic correlations.

3. **Returning the Sparse Matrix**:
   The final matrix \( H^{-1} \) is returned as a sparse matrix, ensuring memory efficiency when working with large datasets.

### Arguments:
- `Ginv`: The inverse of the genomic relationship matrix \( G \).
- `Ainv`: The inverse of the matrix \( A \).
- `A22inv`: The inverse of the submatrix \( A_{22} \) from matrix \( A \).
- `gid`: A vector containing the indices of individuals with genotype information.

### Returns:
- `Hi`: The inverse matrix of \( H \), which has been adjusted for the individuals with genotype data. The result is returned as a sparse matrix to improve memory usage.

"""
function calHInverse(Ginv, Ainv, A22inv, gid)
    Hi = copy(Matrix(Ainv))
    Hi[gid, gid] .+= Ginv - A22inv
    return sparse(Hi)
end
