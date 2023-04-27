# generate the table found in
# https://arxiv.org/abs/2201.05967

using DyadicKDE

# specify parameters
n_data = 3000
n_evals = 50
n_repeats = 2000
degeneracies = ["total", "partial", "none"]
kernel_names = ["epanechnikov_order_2", "epanechnikov_order_4"]
evals = collect(range(-2.0, stop=2.0, length=n_evals))
sdp_solver = "cosmo"
n_resample = 10000
significance_level = 0.05
results = []
ps = Dict("total" => [0.5, 0.0, 0.5],
          "partial" => [0.25, 0.0, 0.75],
          "none" => [0.2, 0.2, 0.6])

println("Running experiments")
for degen in degeneracies
    for kernel_name in kernel_names
        global counter = 0
        Threads.@threads for rep in 1:n_repeats
            global counter += 1
            println("Degeneracy: $degen, Kernel: $kernel_name, Repeat $counter/$n_repeats")
            W = make_dyadic_data(n_data, ps[degen])
            h_ROT = estimate_ROT_bandwidth(W, "epanechnikov_order_2")

            est = DyadicKernelDensityEstimator(kernel_name, h_ROT, significance_level,
                                               n_resample, sdp_solver, evals, W,
                                               Dict("degen" => degen, "p" => ps[degen]))

            fit(est)
            push!(results, est)
        end
    end
end

println("\nResults:")
for degen in degeneracies
    for kernel_name in kernel_names
        results_reduced = [est
                           for est in results
                           if
                           est.meta["degen"] == degen
                           &&
                           est.kernel_name == kernel_name]

        p = ps[degen]
        h_ROT = mean([est.bandwidth for est in results_reduced])
        RIMSE = mean([get_RIMSE(est.fhat, DyadicKDE.get_f(est.meta["p"], est.evals))
                      for est in results_reduced])
        ucb_coverage = mean([get_coverage(est.ucb,
                                          DyadicKDE.get_f(est.meta["p"], est.evals))
                             for est in results_reduced])
        pci_coverage = mean([get_coverage(est.pci,
                                          DyadicKDE.get_f(est.meta["p"], est.evals))
                             for est in results_reduced])
        ucb_width = mean([get_average_width(est.ucb) for est in results_reduced])
        pci_width = mean([get_average_width(est.pci) for est in results_reduced])

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
