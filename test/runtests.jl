using DyadicKDE
using Test
using Random



@testset "ROT bandwidth" begin

    Random.seed!(314159)
    n_data = 100
    p = [0.25, 0.0, 0.75]

    for rep in 1:5

        data = make_data(n_data, p)
        h_ROT = estimate_ROT_bandwidth(data, "epanechnikov_order_2")
        @test 0.5 <= h_ROT <= 0.7
    end
end



@testset "Estimator" begin

    Random.seed!(314159)
    n_data = 100
    kernel_names = ["epanechnikov_order_2", "epanechnikov_order_4"]
    evals = collect(range(-2.0, stop=2.0, length=10))
    sdp_solver = "cosmo"
    n_resample = 1000
    significance_level = 0.5
    bandwidth = 0.5
    ps = [[0.5, 0.0, 0.5],
          [0.25, 0.0, 0.75],
          [0.2, 0.2, 0.6]]

    for rep in 1:5

        for p in ps

            data = make_data(n_data, p)

            for kernel_name in kernel_names

                est = DyadicKernelDensityEstimator(
                    kernel_name, bandwidth, significance_level,
                    n_resample, sdp_solver, evals, data,
                    Dict("p" => p))

                fit(est)

                @test all(0 .<= est.bci[1,:])
                @test all(est.bci[1,:] .<= est.ucb[1,:])
                @test all(est.ucb[1,:] .<= est.pci[1,:])
                @test all(est.pci[1,:] .<= est.fhat)
                @test all(est.fhat .<= est.pci[2,:])
                @test all(est.pci[2,:] .<= est.ucb[2,:])
                @test all(est.ucb[2,:] .<= est.bci[2,:])

                ucb_average_width = get_ucb_average_width(est)
                pci_average_width = get_pci_average_width(est)
                bci_average_width = get_bci_average_width(est)

                @test pci_average_width <= ucb_average_width
                @test ucb_average_width <= bci_average_width

                ucb_coverage = get_ucb_coverage(est)
                pci_coverage = get_pci_coverage(est)
                bci_coverage = get_bci_coverage(est)

                @test pci_coverage <= ucb_coverage
                @test ucb_coverage <= bci_coverage
            end
        end
    end
end
