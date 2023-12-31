module DppInference

using LinearAlgebra
using Distributions
using Random
using Plots

export
    DetrminantalPointProcesses, DPP,
    kDeterminantalPointProcess, kDPP,
    
    # prob
    pmf,
    enumerate_map,
    rand_approx_map,
    mcmc_map,
    
    # conditioning
    condition_exclude,
    condition_include,
    
    # sampling
    sample_exact,

    # visualizing
    viz_hist_exact,
    viz_hist_mcmc,

    # learning
    learn_ka,
    learn_fp,
    learn_em,
    kl_div,

    # synthetic data generation
    build_training_set

# types
include("types.jl")

# methods
include("conditioning.jl")
include("sampling.jl")
include("viz.jl")
include("prob.jl")
include("learn.jl")

# utils
include("utils.jl")

# synthetic data
include("syntheticdatagen/build_training_set.jl")
include("syntheticdatagen/borodin.jl")

end