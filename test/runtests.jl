using DyadicKDE
using Test
using Random

@testset "Confidence bands" begin

    Random.seed!(314159)

    # specify parameters
    n_data = 80
    kernel_names = ["epanechnikov_order_2", "epanechnikov_order_4"]
    evals = collect(range(-2.0, stop=2.0, length=10))
    sdp_solver = "cosmo"
    n_resample = 1000
    significance_level = 0.5
    p = [0.25, 0.0, 0.75]

    for rep in 1:10

        # make data and get bandwidth
        data = make_data(n_data, p)
        h_ROT = estimate_ROT_bandwidth(data, "epanechnikov_order_2")

        # run experiments
        for kernel_name in kernel_names

            est = DyadicKernelDensityEstimator(
                kernel_name, h_ROT, significance_level,
                n_resample, sdp_solver, evals, data, Dict())

            fit(est)

            @test all(
                0
                .<= est.bci[1,:]
                .<= est.ucb[1,:]
                .<= est.pci[1,:]
                .<= est.fhat
                .<= est.pci[2,:]
                .<= est.ucb[2,:]
                .<= est.bci[2,:]
            )

        end

    end

end
