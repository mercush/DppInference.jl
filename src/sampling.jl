
# Functions for sampling from a DPP
function sample_exact(dpp::DeterminantalPointProcess)
    Λ,Q = dpp.Lfact
    mask = rand.(Bernoulli.(Λ./(Λ.+1)))
    return randprojDPP(Q[:,mask])
end

function sample_mcmc(dpp::DeterminantalPointProcess)
    # TODO: https://proceedings.neurips.cc/paper_files/paper/2013/file/20d135f0f28185b84a4cf7aa51f29500-Paper.pdf
end

function randprojDPP(Y::Matrix)
    N = size(Y,2)
    𝓘 = fill(0,N)
    for k=1:N
        p = mean(abs.(Y).^2, dims=2)
        𝓘[k] = rand(Categorical(p[:]))
        Y=(Y*qr(Y[𝓘[k],:]).Q )[:,2:end] 
    end
    return sort(𝓘)
end


# Functions for sampling from a kDPP
function sample_exact(kdpp::kDeterminantalPointProcess)
    Λ, Q = kdpp.Lfact
    e = elementary_symmetric_polynomial(kdpp.k, Λ)
    mask = fill(0,kdpp.N)
    l = kdpp.k
    for n = kdpp.N:-1:1
        if l==0
            break
        end
        if rand.(Bernoulli(Λ[n]*e[n,l]/e[n+1,l+1]))
            mask[n] = 1
            l = l-1
        end
    end
    return randprojDPP(Q[:, map(Bool, mask)])
end