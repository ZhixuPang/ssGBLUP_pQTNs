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
