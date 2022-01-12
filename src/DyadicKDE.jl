module DyadicKDE

using LinearAlgebra
using Distributions
using Statistics
using Convex
using Mosek
using MosekTools
using COSMO

export DyadicKernelDensityEstimator
export fit
#export estimate_fhat
#export estimate_conditional_expectation
#export estimate_Sigmahat
#export estimate_Sigmahatplus
#export estimate_uniform_confidence_band
#export estimate_pointwise_confidence_band
#export estimate_bonferroni_confidence_band
export estimate_ROT_bandwidth

include("estimation.jl")

end
