"""
    calWeights(dat::DataFrame, PI::Float64, n::Int) -> DataFrame

Calculate SNP weights based on local averaging of likelihood ratio (LR) statistics derived from GWAS results.

# Arguments
- `dat::DataFrame`: A DataFrame containing GWAS results, where `dat.effects` holds the estimated SNP effects and `dat.se` contains their corresponding standard errors.
- `PI::Float64`: The prior probability that a SNP is associated with the trait (typically a small value such as 1e-4).
- `n::Int`: Window size for local averaging; determines the number of neighboring SNPs (on either side) included in the smoothing process.

# Returns
- A DataFrame with two columns:
  - `Weights`: The calculated SNP weights, representing the posterior probability that each SNP is associated with the trait.
  - `SmoothedLRvalues`: The locally averaged LR statistics for each SNP.

# Description
The function computes a local average of squared likelihood ratio (LR) statistics, defined as 0.5*(effect/SE)^2, across a specified window around each SNP.
This smoothed LR value is then converted to a posterior inclusion probability using a Bayesian formula that incorporates the prior `PI`.
These posterior probabilities are returned as SNP weights, which can be used to construct a weighted genomic relationship matrix.
"""
function calWeights(dat::DataFrame, PI::Float64, n::Int)
    lr=0.5*(dat.effects./dat.se).^2
    nleft=round(Int,(n-1)/2)
    if(nleft<0)
        nleft=0
    end
    avg=zeros(size(lr,1),2);
    for i=1:size(lr,1)
        lrs = lr[max(1,i-nleft):min(i+nleft,size(lr,1))]
        avg[i,2]=mean(lrs[.!ismissing.(lrs)])
        avg[i,1]=PI*exp(avg[i,2])/(PI*exp(avg[i,2])+1-PI)
    end
    return DataFrame(;Weights = avg[:, 1], SmoothedLRvalues = avg[:, 2])
end
