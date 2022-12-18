module DyadicKDE

using LinearAlgebra
using Distributions
using Statistics
using Convex
using Mosek
using MosekTools
using COSMO

include("kernels.jl")
include("estimation.jl")
include("counterfactual.jl")
include("helpers.jl")

export DyadicKernelDensityEstimator
export estimate_ROT_bandwidth

export CounterfactualDyadicKernelDensityEstimator
export fit

export make_dyadic_data
export mean
export get_RIMSE
export get_coverage
export get_average_width

end
