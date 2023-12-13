function build_training_set(L::DeterminantalPointProcess, n::Int)
    # Returns n samples generated from L
    res = Vector{Int}[]
    for _ in 1:n
        push!(res, sample_exact(L))
    end
    return res
end