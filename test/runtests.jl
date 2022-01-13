using DyadicKDE
using Test
using LinearAlgebra
using Random



function make_data(n_data::Int, p::Vector{Float64})

    c = cumsum(p)
    r = rand(n_data)
    data_A = -1.0 * (r .<= c[1]) + 1.0 * (r .>= c[2])
    data_V = randn((n_data, n_data))
    data = UpperTriangular(data_A .* data_A' .+ data_V)

    return data
end



function mean(x::Vector)

    return sum(x) / length(x)
end



function phi(t::Float64)

    return (2 * pi)^(-0.5) * exp(-(t^2) / 2)
end



function get_f(est::DyadicKernelDensityEstimator)

    p = est.meta["p"]
    f = (p[1]^2 + p[3]^2) * phi.(est.evals .- 1)
    f += p[2]*(2 - p[2]) * phi.(est.evals)
    f += 2*p[1]*p[3] * phi.(est.evals .+ 1)
    return f
end



function get_RIMSE(est::DyadicKernelDensityEstimator)

    f = get_f(est)
    RMSE = mean((est.fhat .- f).^2)
    return RMSE
end



function get_ucb_coverage(est::DyadicKernelDensityEstimator)

    return all(est.ucb[1,:] .<= get_f(est) .<= est.ucb[2,:])
end



function get_pci_coverage(est::DyadicKernelDensityEstimator)

    return all(est.pci[1,:] .<= get_f(est) .<= est.pci[2,:])
end



function get_ucb_average_width(est::DyadicKernelDensityEstimator)

    return sum(est.ucb[2,:] .- est.ucb[1,:]) / est.n_evals
end



function get_pci_average_width(est::DyadicKernelDensityEstimator)

    return sum(est.pci[2,:] .- est.pci[1,:]) / est.n_evals
end



@testset "paper_replication" begin

    # specify parameters
    n_data = 500
    n_repeats = 2
    degeneracies = ["total", "partial", "none"]
    kernel_names = ["epanechnikov_order_2", "epanechnikov_order_4"]
    evals = collect(range(-2.0, stop=2.0, length=30))
    sdp_solver = "cosmo"
    n_resample = 10000
    significance_level = 0.05
    results = []
    ps = Dict(
        "total" => [0.5, 0.0, 0.5],
        "partial" => [0.25, 0.0, 0.75],
        "none" => [0.2, 0.2, 0.6],
    )

    # run experiments
    Random.seed!(314159)
    for degen in degeneracies
        for kernel_name in kernel_names

            global counter = 0
            Threads.@threads for rep in 1:n_repeats

                global counter += 1
                println("$degen, $kernel_name, $counter / $n_repeats")
                data = make_data(n_data, ps[degen])
                h_ROT = estimate_ROT_bandwidth(data, "epanechnikov_order_2")

                est = DyadicKernelDensityEstimator(
                    kernel_name, h_ROT, significance_level,
                    n_resample, sdp_solver, evals, data,
                    Dict("degen" => degen, "p" => ps[degen]))

                fit(est)
                push!(results, est)

            end
        end
    end


    println("\nResults:")
    for degen in degeneracies
        for kernel_name in kernel_names

            results_reduced = [est for est in results if
                                   est.meta["degen"] == degen
                                   && est.kernel_name == kernel_name]

            p = ps[degen]
            h_ROT = mean([est.bandwidth for est in results_reduced])
            RIMSE = mean([get_RIMSE(est) for est in results_reduced])
            ucb_coverage = mean([get_ucb_coverage(est) for est in results_reduced])
            pci_coverage = mean([get_pci_coverage(est) for est in results_reduced])
            ucb_width = mean([get_ucb_average_width(est) for est in results_reduced])
            pci_width = mean([get_pci_average_width(est) for est in results_reduced])

            println("π = $p")
            println("Degeneracy: $degen")
            println("h_ROT = $h_ROT")
            println("Kernel: $kernel_name")
            println("RIMSE: $RIMSE")
            println("UCB coverage: $ucb_coverage")
            println("UCB width: $ucb_width")
            println("PCI coverage: $pci_coverage")
            println("PCI width: $pci_width")
            println()

        end
    end

end
