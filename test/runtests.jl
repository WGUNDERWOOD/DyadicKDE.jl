using DyadicKDE
using Test
using LinearAlgebra

@testset "DyadicKDE.jl" begin

    # make evals
    n_evals = 10
    evals = collect(range(-2.0, stop=2.0, length=n_evals))

    # make dyadic data
    n_data = 20
    p = [0.25, 0, 0.75]
    r = rand(n_data)
    c = cumsum(p)
    data_A = -1.0 * (r .<= c[1]) + 1.0 * (r .>= c[2])
    data_V = UpperTriangular(randn((n_data, n_data)))
    data = UpperTriangular(data_A .* data_A' .* data_V)

    # fit estimator
    bandwidth = estimate_ROT_bandwidth(data, "epanechnikov_order_2")
    kernel_name = "epanechnikov_order_2"
    significance_level = 0.05
    n_resample = 100
    sdp_solver = "cosmo"

    estimator = DyadicKernelDensityEstimator(
        kernel_name, bandwidth, significance_level,
        n_resample, sdp_solver, evals, data)

    fit(estimator)

end
