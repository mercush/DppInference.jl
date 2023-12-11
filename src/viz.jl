function viz_hist_exact(dpp::DeterminantalPointProcess)
    v = zeros(2^dpp.N)
    for set in _iter_subset(dpp.N)
        set_int = _subset_to_int(set)
        v[set_int] = pmf(dpp, set)
    end
    bar((0:2^dpp.N-1),v)
end

function viz_hist_mcmc(dpp::DeterminantalPointProcess, sample_type="exact", iters=5000)
    data = []
    for _ in 1:iters
        if sample_type == "exact"
            r = sample_exact(dpp)
        elseif sample_type == "approx"
            r = sample_mcmc(dpp)
        end
        r_int = _subset_to_int(r)
        append!(data, r_int)
    end
    histogram(data, label="DPP PMF")
end

