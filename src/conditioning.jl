
# Functions for doing conditional inference on an DPP
function _borodin_transform(dpp::DeterminantalPointProcess, A::Vector{Int})
    Abar = complement(dpp.N, A)
    return Symmetric((((dpp.L+restrictI(dpp.N, Abar))^-1)[Abar, Abar])^-1-I)
end

function condition_exclude(dpp::DeterminantalPointProcess, A::Vector{Int})
    Abar = complement(dpp.N, A)
    Lp = Symmetric(dpp.L[Abar, Abar])
    return DeterminantalPointProcess(Lp)
end

function condition_include(dpp::DeterminantalPointProcess, A::Vector{Int})
    LA = _borodin_transform(dpp, A)
    return DeterminantalPointProcess(LA)
end

# sample from kDPP
function _borodin_transform(kdpp::kDeterminantalPointProcess, A::Vector{Int})
    Abar = collect(setdiff(Set(1:kdpp.N), A))
    return Symmetric((((L+restrictI(dpp.N, Abar))^-1)[Abar, Abar])^-1-I)
end

function condition_include(kdpp::kDeterminantalPointProcess, A::Vector{Int})
    LA = _borodin_transform(kdpp, A)
    return kDeterminantalPointProcess(k, LA)
end