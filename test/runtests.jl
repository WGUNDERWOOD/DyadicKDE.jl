using DyadicKDE
using Test
using LinearAlgebra
using Random



function make_data(n_data::Int, p::Vector{Float64})

    c = cumsum(p)
    r = rand(n_data)
    data_A = -1.0 * (r .<= c[1]) + 1.0 * (r .>= c[2])
    data_V = UpperTriangular(randn((n_data, n_data)))
    data = UpperTriangular(data_A .* data_A' .* data_V)

    return data
end



function mean(x::Vector)

    return sum(x) / length(x)
end



@testset "paper_replication" begin

    # specify parameters
    n_data_plot = 100
    n_data_table = 50
    n_repeats = 2
    degeneracies = ["total", "partial", "none"]
    kernel_names = ["epanechnikov_order_2", "epanechnikov_order_4"]
    evals_plot = collect(range(-2.0, stop=2.0, length=100))
    evals_table = collect(range(-2.0, stop=2.0, length=50))
    sdp_solver = "mosek"
    n_resample = 10000
    significance_level = 0.05
    results_plot = []
    results_table = []
    ps = Dict(
        "total" => [0.5, 0.0, 0.5],
        "partial" => [0.25, 0.0, 0.75],
        "none" => [0.2, 0.2, 0.6],
    )

    # run plot experiments
    #Random.seed!(314159)
    #for degen in degeneracies
        #for kernel_name in kernel_names

            #global counter = 0
            #Threads.@threads for rep in 1:n_repeats

                #global counter += 1
                #println("Degeneracy: $degen, Kernel: $kernel_name, Repeat: $counter/$n_repeats")
                #data = make_data(n_data_plot, ps[degen])
                #h_ROT = estimate_ROT_bandwidth(data, "epanechnikov_order_2")

                #est = DyadicKernelDensityEstimator(
                    #kernel_name, h_ROT, significance_level,
                    #n_resample, sdp_solver, evals_plot, data,
                    #Dict("degen" => degen))

                #fit(est)
                #push!(results_plot, est)

            #end
        #end
    #end

    # run table experiments
    Random.seed!(314159)
    for degen in degeneracies
        for kernel_name in kernel_names

            global counter = 0
            Threads.@threads for rep in 1:n_repeats

                global counter += 1
                println("$degen, $kernel_name, $rep / $n_repeats")
                data = make_data(n_data_table, ps[degen])
                h_ROT = estimate_ROT_bandwidth(data, "epanechnikov_order_2")

                est = DyadicKernelDensityEstimator(
                    kernel_name, h_ROT, significance_level,
                    n_resample, sdp_solver, evals_table, data,
                    Dict("degen" => degen))

                fit(est)
                push!(results_table, est)

                println(est.meta)

            end
        end
    end


    println("\nResults:")
    for degen in degeneracies
        for kernel_name in kernel_names

            results_reduced = [est for est in results_table if
                                   est.meta["degen"] == degen
                                   && est.kernel_name == kernel_name]

            p = ps[degen]
            h_ROT = mean([est.bandwidth for est in results_reduced])

            println("π = $p")
            println("Degeneracy: $degen")
            println("h_ROT = $h_ROT")
            println("Kernel: $kernel_name")

            # TODO add RIMSE and coverage rates

            println()

        end
    end


end
