
function loss(Λ, V, Y)

end

# With dpp as a starting point, learns the parameters for a distribution Y. Terminating condition given by c.
function learnEM(N, Y::Vector{Set{Int64}}, c)
    # TODO: https://arxiv.org/pdf/1411.1088.pdf
    K =1/2*Matrix(I,N,N)
    Λ, V = eigen(K)
    function avgPk(Λ, V, j)
        R = diag(√(Λ./(1-Λ)))
        U = V'
        for y in Y
            U_y = U[:, y]
            Z_Y = U_y * U_y'
            Q_Y = R*Z_Y*R
            _, V_hat = eigen(Q_Y)
            p_K(j) -> sum(V_hat[:,j].*V_hat[:,j])/length(Y)
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
    while true
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
        δ < c || break
    end
end

function learnFP(dpp::DeterminantalPointProcess, Y::Vector{Set{Int64}}, c)
    # TODO: https://arxiv.org/pdf/1508.00792.pdf
end