using DyadicKDE
using Test
using Random
using Suppressor
using Aqua
using JuliaFormatter

function trapezium_integrate(f::Vector{Float64}, x::Vector{Float64})
    @assert length(f) == length(x)
    n_areas = length(x) - 1
    areas = fill(NaN, n_areas)

    for i in 1:n_areas
        areas[i] = 0.5 * (f[i + 1] + f[i]) * (x[i + 1] - x[i])
    end

    area = sum(areas)
    return area
end

# Aqua tests
Aqua.test_ambiguities(DyadicKDE)
Aqua.test_unbound_args(DyadicKDE)
Aqua.test_undefined_exports(DyadicKDE)
Aqua.test_project_extras(DyadicKDE)
Aqua.test_stale_deps(DyadicKDE, ignore=[:Aqua, :Suppressor, :JuliaFormatter])
Aqua.test_deps_compat(DyadicKDE)

@testset verbose = true "JuliaFormatter" begin
    @test JuliaFormatter.format(DyadicKDE, overwrite=false)
end

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
            @test abs(trapezium_integrate(cs .^ 2 .* ks, ws) - 0) <= tol
            @test abs(trapezium_integrate(cs .^ 3 .* ks, ws) - 0) <= tol
        end
    end
end

@testset "ROT bandwidth" begin
    Random.seed!(314159)
    n_data = 1000
    p = [0.25, 0.0, 0.75]

    for rep in 1:5
        W = make_dyadic_data(n_data, p)
        h_ROT = estimate_ROT_bandwidth(W, "epanechnikov_order_2")
        @test 0.24 <= h_ROT <= 0.25
    end
end

@testset "Estimator" begin
    Random.seed!(314159)
    n_data = 50
    kernel_names = ["epanechnikov_order_2", "epanechnikov_order_4"]
    evals = collect(range(-2.0, stop=2.0, length=3))
    sdp_solver = "cosmo"
    n_resample = 1000
    significance_level = 0.5
    bandwidth = 0.8
    ps = [[0.5, 0.0, 0.5],
          [0.25, 0.0, 0.75],
          [0.2, 0.2, 0.6]]

    for rep in 1:5
        for p in ps
            W = make_dyadic_data(n_data, p)

            for kernel_name in kernel_names
                est = DyadicKernelDensityEstimator(kernel_name, bandwidth, significance_level,
                                                   n_resample, sdp_solver, evals, W, Dict("p" => p))

                fit(est)

                @test all(est.bci[1, :] .< est.ucb[1, :])
                @test all(est.ucb[1, :] .< est.pci[1, :])
                @test all(est.pci[1, :] .< est.fhat)
                @test all(est.fhat .< est.pci[2, :])
                @test all(est.pci[2, :] .< est.ucb[2, :])
                @test all(est.ucb[2, :] .< est.bci[2, :])

                ucb_average_width = get_average_width(est.ucb)
                bci_average_width = get_average_width(est.bci)
                pci_average_width = get_average_width(est.pci)

                @test pci_average_width < ucb_average_width
                @test ucb_average_width < bci_average_width

                ucb_coverage = get_coverage(est.ucb, DyadicKDE.get_f(p, evals))
                bci_coverage = get_coverage(est.bci, DyadicKDE.get_f(p, evals))
                pci_coverage = get_coverage(est.pci, DyadicKDE.get_f(p, evals))

                @test pci_coverage <= ucb_coverage
                @test ucb_coverage <= bci_coverage

                RIMSE = get_RIMSE(est.fhat, DyadicKDE.get_f(p, evals))
                @test 0 <= RIMSE <= 0.1

                @suppress display(est)
            end
        end
    end
end

@testset "Counterfactual" begin
    Random.seed!(314159)
    n_data = 50
    kernel_names = ["epanechnikov_order_2", "epanechnikov_order_4"]
    evals = collect(range(-2.0, stop=2.0, length=10))
    sdp_solver = "cosmo"
    n_resample = 1000
    significance_level = 0.5
    bandwidth = 0.8
    p = [0.2, 0.2, 0.6]

    for rep in 1:5
        W = make_dyadic_data(n_data, p)
        X0 = [Int(1 + round(2 * rand())) for _ in 1:n_data]
        X1 = [Int(1 + round(2 * rand())) for _ in 1:n_data]

        for kernel_name in kernel_names
            est = CounterfactualDyadicKernelDensityEstimator(kernel_name, bandwidth,
                                                             significance_level, n_resample,
                                                             sdp_solver, evals, W, X0, X1,
                                                             Dict("p" => p))

            fit(est)

            @test all(est.bci[1, :] .< est.ucb[1, :])
            @test all(est.ucb[1, :] .< est.pci[1, :])
            @test all(est.pci[1, :] .< est.fhat)
            @test all(est.fhat .< est.pci[2, :])
            @test all(est.pci[2, :] .< est.ucb[2, :])
            @test all(est.ucb[2, :] .< est.bci[2, :])

            ucb_average_width = get_average_width(est.ucb)
            pci_average_width = get_average_width(est.pci)
            bci_average_width = get_average_width(est.bci)

            @test pci_average_width < ucb_average_width
            @test ucb_average_width < bci_average_width

            ucb_coverage = get_coverage(est.ucb, DyadicKDE.get_f(p, evals))
            pci_coverage = get_coverage(est.pci, DyadicKDE.get_f(p, evals))
            bci_coverage = get_coverage(est.bci, DyadicKDE.get_f(p, evals))

            @test pci_coverage <= ucb_coverage
            @test ucb_coverage <= bci_coverage

            RIMSE = get_RIMSE(est.fhat, DyadicKDE.get_f(p, evals))
            @test 0 <= RIMSE <= 0.05

            @suppress display(est)
        end
    end
end

@testset "ParametricCounterfactual" begin
    Random.seed!(314159)
    n_data = 50
    kernel_names = ["epanechnikov_order_2", "epanechnikov_order_4"]
    evals = collect(range(-2.0, stop=2.0, length=10))
    sdp_solver = "cosmo"
    n_resample = 1000
    significance_level = 0.5
    bandwidth = 0.8
    p = [0.2, 0.2, 0.6]

    for rep in 1:5
        W = make_dyadic_data(n_data, p)
        X0 = [Int(1 + round(2 * rand())) for _ in 1:n_data]
        X1 = [Int(1 + round(2 * rand())) for _ in 1:n_data]
        phat0 = [0.4, 0.3, 0.3]
        phat1 = [0.3, 0.4, 0.3]

        for kernel_name in kernel_names
            est = ParametricCounterfactualDyadicKernelDensityEstimator(kernel_name, bandwidth,
                                                                       significance_level,
                                                                       n_resample,
                                                                       sdp_solver, evals, W, X0, X1,
                                                                       phat0, phat1,
                                                                       Dict("p" => p))

            fit(est)

            @test all(est.bci[1, :] .< est.ucb[1, :])
            @test all(est.ucb[1, :] .< est.pci[1, :])
            @test all(est.pci[1, :] .< est.fhat)
            @test all(est.fhat .< est.pci[2, :])
            @test all(est.pci[2, :] .< est.ucb[2, :])
            @test all(est.ucb[2, :] .< est.bci[2, :])

            ucb_average_width = get_average_width(est.ucb)
            pci_average_width = get_average_width(est.pci)
            bci_average_width = get_average_width(est.bci)

            @test pci_average_width < ucb_average_width
            @test ucb_average_width < bci_average_width

            ucb_coverage = get_coverage(est.ucb, DyadicKDE.get_f(p, evals))
            pci_coverage = get_coverage(est.pci, DyadicKDE.get_f(p, evals))
            bci_coverage = get_coverage(est.bci, DyadicKDE.get_f(p, evals))

            @test pci_coverage <= ucb_coverage
            @test ucb_coverage <= bci_coverage

            RIMSE = get_RIMSE(est.fhat, DyadicKDE.get_f(p, evals))
            @test 0 <= RIMSE <= 0.05

            @suppress display(est)
        end
    end
end

@testset "Errors" begin
    Random.seed!(314159)

    @test_throws ErrorException("Unknown kernel_name") DyadicKDE.kernel(0.0, 0.0, 1.0, -1.0, 1.0,
                                                                        "not_a_valid_kernel_name")

    evals = collect(range(-2.0, stop=2.0, length=10))
    p = [0.5, 0.0, 0.5]
    W = make_dyadic_data(50, p)
    est = DyadicKernelDensityEstimator("epanechnikov_order_2", 0.8, 0.5, 10,
                                       "not_a_valid_sdp_solver_name", evals, W, Dict("p" => p))

    @test_throws ErrorException("Unknown sdp_solver") fit(est)

    @test_throws ErrorException("Unknown kernel_name") estimate_ROT_bandwidth(W,
                                                                              "epanechnikov_order_4")
end
