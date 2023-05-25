using DyadicKDE

REPDIR = @__DIR__() * "/"
PLOTDIR = REPDIR * "plots/"
DATADIR = REPDIR * "data/"

include(REPDIR * "plot_helpers.jl")

# specify parameters
n_data = 100
n_evals = 100
degeneracies = ["total", "partial", "none"]
kernel_name = "epanechnikov_order_2"
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
    println("Degeneracy: $degen")
    W = make_dyadic_data(n_data, ps[degen])
    h_ROT = estimate_ROT_bandwidth(W, "epanechnikov_order_2")

    est = DyadicKernelDensityEstimator(kernel_name, h_ROT, significance_level,
                                       n_resample, sdp_solver, evals, W,
                                       Dict("degen" => degen, "p" => ps[degen]))

    fit(est)
    push!(results, est)
end

println("Making plots")
linewidth = 1.0
y_lim = [-0.01, 0.45]
e = 0.22
pci_marker = PyPlot.matplotlib.path.Path([[-1, 0], [-1, e], [-1, -e], [-1, 0], [1, 0], [1, e],
                                          [1, -e], [1, 0]])
handle_f = PyPlot.matplotlib.lines.Line2D([0], [0], color="k", lw=linewidth, label="\$f_W(w)\$")
handle_fhat = PyPlot.matplotlib.lines.Line2D([0], [0], color="k", lw=linewidth,
                                             linestyle=(0, (1, 1)), label="\$\\widehat f_W(w)\$")
handle_ucb = PyPlot.matplotlib.patches.Patch(facecolor="lightgray", edgecolor="lightgray",
                                             label="UCB")
handle_pci = PyPlot.matplotlib.lines.Line2D([0], [0], color="black", lw=0, markersize=17,
                                            marker=pci_marker, label="PCI")
handles = [handle_f, handle_fhat, handle_ucb, handle_pci]

for est in results
    degen = est.meta["degen"]
    p = est.meta["p"]
    f = generate_f(est.evals, p)

    # plot f
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.plot(est.evals, f, color="black", linewidth=linewidth)

    # plot fhat
    ax.plot(est.evals, est.fhat,
            color="black", linewidth=linewidth, linestyle=(0, (1, 1)))

    # plot ucb
    ax.fill_between(est.evals, est.ucb[1, :], est.ucb[2, :],
                    color="lightgray", linewidth=0.0)

    # plot pci
    plot_confidence_intervals(ax, est.evals, est.pci[1, :],
                              est.pci[2, :], 10, "black", linewidth, "-")

    # save plot
    PyPlot.xlabel("Evaluation point, \$w\$")
    plt.ylim(y_lim)
    x_min = minimum(est.evals)
    x_max = maximum(est.evals)
    plt.xlim((x_min, x_max))
    plt.yticks(range(0.0, stop=0.4, step=0.1))
    legend(handles=handles, loc="upper left")
    plt.ylabel("Density", labelpad=4.0)
    plt.tight_layout()
    PyPlot.savefig(PLOTDIR * "outcome_plot_$degen.pdf")
    close("all")
end
