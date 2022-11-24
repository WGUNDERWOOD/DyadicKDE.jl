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
export DyadicKernelDensityEstimator
export estimate_ROT_bandwidth

include("helpers.jl")
export make_data
export mean
export get_RIMSE
export get_ucb_coverage
export get_pci_coverage
export get_bci_coverage
export get_ucb_average_width
export get_pci_average_width
export get_bci_average_width

include("counterfactual.jl")
export CounterfactualDyadicKernelDensityEstimator
export fit

end
