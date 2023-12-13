# NLL loss function with respect to L and K
function mle_loss_K(Y::Vector{Vector{Int}}, K::Symmetric)
    N = size(K, 1)
    return sum(map(y->log(abs(det(K-restrictI(N, complement(N, y))))), Y))
end

function mle_loss_L(Y::Vector{Vector{Int}}, L::Symmetric)
    n = size(Y, 1)
    return sum(map(y->log(det(L[y,y])), Y))-n*log(det(I+L))
end

# Derivative of MLE with respect to K
function dLdK(Y::Vector{Vector{Int}}, K::Symmetric)
    N = size(K, 1)
    return Symmetric(sum(map(y->(K-restrictI(N, complement(N, y)))^-1, Y)))
end

# KL divergence of two DPPs
function kl_div(dppP::DeterminantalPointProcess, dppQ::DeterminantalPointProcess)
    # dppP.N == dppQ.N
    s = 0
    for set in _iter_subset(dppP.N)
        s += pmf(dppP, set)*log(pmf(dppP, set)/pmf(dppQ, set))
    end
    return s
end

# Helper functions for FP learning
function fixed_point_map(Y::Vector{Vector{Int}}, L::Symmetric, a::Float64)
    Z = zeros(size(L))
    n = size(Y, 1)
    for y in Y
        if y != []
            Z[y,y] = Z[y,y] + inv(Symmetric(L[y,y])) # TODO: condition number of matrix L is too high. implement inverse function that is more robust to this.
        end
    end
    res = Symmetric(L+a/n*L*Z*L-L/(I+L)*L)
    display(eigen(res))
    return res
end

# Helper functions for EM learning
function avgPk(Λ, V, j)
    R = diag(√(Λ./(1-Λ)))
    U = V'
    for y in Y
        U_y = U[:, y]
        Z_Y = U_y * U_y'
        Q_Y = R*Z_Y*R
        _, V_hat = eigen(Q_Y)
        p_K = j -> sum(V_hat[:,j].*V_hat[:,j])/length(Y)
    end
end

function deltaF(Λ, V)
    res = 0
    for y in Y
        U_y = U[:, y]
        B_Y = Matrix(I, N, N)[:, y]
        H_Y = U_y'*R*R*U_y
        res += 2 * B_Y * H_Y^-1 * U_y' * R * R
    end
end

# Learning algorithms
function learn_ka(Y::Vector{Vector{Int}}, N::Int, n_iter::Int)
    L = L_wishart_init(N)
    K = Symmetric(L/(I+L))
    for i in 1:n_iter
        G = dLdK(Y, K)
        η = 1
        Q = zeros(N, N)
        while true
            Q = K + η*G
            Λ, V = eigen(Q)
            Λ = min.(max.(Λ, 0), 1)
            Q = Symmetric(V*diagm(Λ)*V')
            η = η/2
            mle_loss_K(Y, Q) < mle_loss_K(Y, K) || break
        end
        K = Q
    end
    return K
end

function learn_fp(Y::Vector{Vector{Int}}, N::Int, a::Float64, n_iter::Int)
    L = L_wishart_init(N)
    for _ in 1:n_iter
        L = fixed_point_map(Y, L, a)
        println(mle_loss_L(Y, L))
    end
    return L
end

function learn_em(Y::Vector{Vector{Int}}, N::Int, n_iter::Int)
    # TODO: https://arxiv.org/pdf/1411.1088.pdf
    L = L_wishart_init(N)
    K = L/(I+L)
    Λ, V = eigen(K)
    for _ in 1:n_iter
        Λ_prime = zeros(N)
        for j in 1:N
            Λ_prime[j] = avePk(Λ, V, j)
        end
        G = deltaF(Λ, F)
        while true
            V_prime = V*exp(η*(V'*G-G'*V))
            η = η/2
            loss(V', Λ_prime) > loss(V, Λ_prime) || break
        end
        δ = 
        Λ = Λ_prime
        V = V_prime
        η = 2*η
    end
end
