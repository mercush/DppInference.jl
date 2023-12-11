# Functions for getting the pmf of a DPP
function _unnormalized_pmf(dpp::DeterminantalPointProcess, A::Vector)
    return det(dpp.L[A, A])
end

function pmf(dpp::DeterminantalPointProcess, A::Vector)
    return _unnormalized_pmf(dpp, A) / dpp.Z
end

# Functions for getting the pmf of a kDPP
function _unnormalized_pmf(kdpp::kDeterminantalPointProcess, A::Vector)
    return det(kdpp.L[A, A])
end

function pmf(kdpp::kDeterminantalPointProcess, A::Vector)
    if length(A) != kdpp.k
        return 0
    end
    return _unnormalized_pmf(kdpp, A) / kdpp.Z
end

# Functions for getting the MAP of a DPP
function enumerate_map(dpp::DeterminantalPointProcess)
    map = []
    ma = 1
    for set in _iter_subset(dpp.N)
        c = _unnormalized_pmf(dpp, set)
        if c > ma
            ma = c
            map = set
        end
    end
    return (map, ma/dpp.Z)
end

function rand_approx_map(dpp::DeterminantalPointProcess, iters=500)
    map = 0 
    ma = 1
    for _ in 1:iters
        r = randsubseq(1:dpp.N, 0.5)
        c = _unnormalized_pmf(dpp, r)
        if c > ma
            ma = c
            map = r
        end
    end
    return (map, ma/dpp.Z)
end

function mcmc_map(dpp::DeterminantalPointProcess, sample_type="exact", iters=500)
    seen = Dict()
    ct = 0
    map = []
    for _ in 1:500
        if sample_type == "exact"
            r = sample_exact(dpp)
        elseif sample_type == "mcmc"
            r = sample_mcmc(dpp)
        end
        r_int =_subset_to_int(r)
        seen[r_int] = get(seen, r_int, 0) + 1
        c = seen[r_int]
        if c > ct
            ct = c
            map = r
        end
    end
    return (map, ct/iters)
end

function greedy_map(dpp::DeterminantalPointProcess)
    # TODO: https://arxiv.org/pdf/1709.05135.pdf
end

# Functions for getting the MAP of a kDPP
function enumerate_map(kdpp::kDeterminantalPointProcess)
    # Search for subsets only among the ones that are length k
    map = []
    ma = 1
    for set in _iter_subset(kdpp.N)
        c = _unnormalized_pmf(kdpp, set)
        if c > ma
            ma = c
            map = set
        end
    end
    return (map, ma/kdpp.Z)
end

function rand_approx_map(kdpp::kDeterminantalPointProcess, iters=500)
    map = 0 
    ma = 1
    for _ in 1:iters
        r = randsubseq(1:kdpp.N, 0.5)
        c = _unnormalized_pmf(kdpp, r)
        if c > ma
            ma = c
            map = r
        end
    end
    return (map, ma/kdpp.Z)
end

function mcmc_map(kdpp::kDeterminantalPointProcess, sample_type="exact", iters=500)
    seen = Dict()
    ma = 0
    for _ in 1:500
        if sample_type == "exact"
            r = sample_exact(kdpp)
        elseif sample_type == "mcmc"
            r = sample_mcmc(kdpp)
        end
        seen[r] = seen.get(r, 0) + 1
        c = seen[r]
        if c > ma
            ma = c
            map = r
        end
    end
    return (map, ma/iters)
end