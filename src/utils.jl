function restrictI(N::Int, A::Vector{Int})
    IA = zeros(N, N)
    for i in A
        IA[i,i] = 1
    end
    return IA
end

function elementary_symmetric_polynomial(k::Int, Λ::Vector)
    N = length(Λ)
    e = zeros(N+1, k+1) # e^n_l = e[n+1,l+1]
    for n in 0:N
        e[n+1,1] = 1
    end
    for l in 1:k
        for n in 1:N
            e[n+1,l+1] = e[n,l+1] + Λ[n]*e[n,l]
        end
    end
    return e
end

function _iter_subset(N::Int)
    if N == 0
        return [[]]
    end
    prev_iter = _iter_subset(N-1)
    return cat([[x; [N]] for x in prev_iter], [x for x in prev_iter], dims=1)
end

function _subset_to_int(A::Vector)
    res = 1
    for i in A
        res += 2^(i-1)
    end
    return res
end