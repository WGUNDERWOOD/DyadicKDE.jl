using DyadicKDE
using Test
using Random

@testset "small_example" begin

    # specify parameters
    n_data = 50
    kernel_names = ["epanechnikov_order_2", "epanechnikov_order_4"]
    evals = collect(range(-2.0, stop=2.0, length=10))
    sdp_solver = "cosmo"
    n_resample = 1000
    significance_level = 0.05
    p = [0.25, 0.0, 0.75]

    # make data and get bandwidth
    data = make_data(n_data, p)
    h_ROT = estimate_ROT_bandwidth(data, "epanechnikov_order_2")

    # run experiments
    Random.seed!(314159)
    for kernel_name in kernel_names

        est = DyadicKernelDensityEstimator(
            kernel_name, h_ROT, significance_level,
            n_resample, sdp_solver, evals, data, Dict())

        fit(est)

    end

end
