using DyadicKDE
using Test
using Random
using Suppressor
using Aqua



function trapezium_integrate(f::Vector{Float64}, x::Vector{Float64})

    @assert length(f) == length(x)
    n_areas = length(x) - 1
    areas = fill(NaN, n_areas)

    for i in 1:n_areas
        areas[i] = 0.5 * (f[i+1] + f[i]) * (x[i+1] - x[i])
    end

    area = sum(areas)
    return area
end



# Aqua tests
Aqua.test_ambiguities(DyadicKDE)
Aqua.test_unbound_args(DyadicKDE)
Aqua.test_undefined_exports(DyadicKDE)
Aqua.test_project_extras(DyadicKDE)
Aqua.test_stale_deps(DyadicKDE, ignore=[:Aqua, :Suppressor])
Aqua.test_deps_compat(DyadicKDE)
Aqua.test_project_toml_formatting(DyadicKDE)



@testset verbose = true "Kernels" begin

    h = 0.1
    w_min = 0.0
    w_max = 1.0
    tol = 1e-4
    len = 4000

    ws = collect(range(0, stop=1, length=len))

    @testset "Epanechnikov order 2" begin
        for i in 1:len
            w = ws[i]
            ks = [DyadicKDE.kernel(s, w, h, w_min, w_max, "epanechnikov_order_2") for s in ws]
            cs = (ws .- w) ./ h
            @test abs(trapezium_integrate(ks, ws) - 1) <= tol
            @test abs(trapezium_integrate(cs .* ks, ws) - 0) <= tol
        end
    end

    @testset "Epanechnikov order 4" begin
        for i in 1:len
            w = ws[i]
            ks = [DyadicKDE.kernel(s, w, h, w_min, w_max, "epanechnikov_order_4") for s in ws]
            cs = (ws .- w) ./ h
            @test abs(trapezium_integrate(ks, ws) - 1) <= tol
            @test abs(trapezium_integrate(cs .* ks, ws) - 0) <= tol
            @test abs(trapezium_integrate(cs.^2 .* ks, ws) - 0) <= tol
            @test abs(trapezium_integrate(cs.^3 .* ks, ws) - 0) <= tol
        end
    end
end



@testset "ROT bandwidth" begin

    Random.seed!(314159)
    n_data = 1000
    p = [0.25, 0.0, 0.75]

    for rep in 1:5
        data = make_data(n_data, p)
        h_ROT = estimate_ROT_bandwidth(data, "epanechnikov_order_2")
        @test 0.24 <= h_ROT <= 0.25
    end
end



@testset "Estimator" begin

    Random.seed!(314159)
    n_data = 50
    kernel_names = ["epanechnikov_order_2", "epanechnikov_order_4"]
    evals = collect(range(-2.0, stop=2.0, length=10))
    sdp_solver = "cosmo"
    n_resample = 1000
    significance_level = 0.5
    bandwidth = 0.8
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

                RIMSE = get_RIMSE(est)
                @test 0 <= RIMSE <= 0.05

                @suppress display(est)
            end
        end
    end
end



@testset "Errors" begin

    @test_throws "Unknown kernel_name" DyadicKDE.kernel(
        0.0, 0.0, 1.0, -1.0, 1.0, "not_a_valid_kernel_name")

    evals = collect(range(-2.0, stop=2.0, length=10))
    p = [0.5, 0.0, 0.5]
    data = make_data(50, p)
    est = DyadicKernelDensityEstimator(
        "epanechnikov_order_2", 0.8, 0.5,
        10, "not_a_valid_sdp_solver_name", evals, data,
        Dict("p" => p))

    @test_throws "Unknown sdp_solver" fit(est)

    @test_throws "Unknown kernel_name" estimate_ROT_bandwidth(
        data, "epanechnikov_order_4")
end
