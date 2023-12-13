abstract type PointProcess end;

mutable struct DeterminantalPointProcess <: PointProcess
    L::Symmetric
    Lfact::Eigen
    K::Symmetric
    Z::Float64
    N::Int

    function DeterminantalPointProcess(L::Symmetric)
        Lfact = eigen(L)
        K = Symmetric(L/(L+I))
        Z = det(L+I)
        N = size(L, 1)
        new(L, Lfact, K, Z, N)
    end
end

mutable struct kDeterminantalPointProcess <: PointProcess
    k::Int
    L::Symmetric
    Lfact::Eigen
    Z::Float64
    N::Int

    function kDeterminantalPointProcess(k::Int, L::Symmetric)
        Lfact = eigen(L)
        N = size(L, 1)
        Z = elementary_symmetric_polynomial(k, Lfact.values)[N+1,k+1]
        new(k, L, Lfact, Z, N)
    end
end

const DPP = DeterminantalPointProcess
const kDPP = kDeterminantalPointProcess