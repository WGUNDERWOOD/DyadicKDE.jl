# generate the trade plots found in
# https://arxiv.org/abs/2201.05967

using CSV
using DataFrames
using LinearAlgebra
using DyadicKDE
include("./plot_helpers.jl")


# set estimation parameters
n_evals = 100
kernel_name = "epanechnikov_order_4"
evals = collect(range(-10.0, stop=10.0, length=n_evals))
sdp_solver = "cosmo"
n_resample = 10000
significance_level = 0.05


# set plot parameters
linewidth = 1.0
y_lim = [-0.003, 0.08]
e = 0.22
pci_marker = PyPlot.matplotlib.path.Path([[-1,0], [-1,e], [-1,-e], [-1,0], [1,0], [1,e], [1,-e], [1,0]])
handle_fhat = PyPlot.matplotlib.lines.Line2D([0], [0], color="k", lw=linewidth, linestyle=(0, (1,1)), label="\$\\widehat f_W(w)\$")
handle_ucb = PyPlot.matplotlib.patches.Patch(facecolor="lightgray", edgecolor="lightgray", label="UCB")
handle_pci = PyPlot.matplotlib.lines.Line2D([0], [0], color="black", lw=0, markersize=17, marker=pci_marker, label="PCI")
handles = [handle_fhat, handle_ucb, handle_pci]


println("Running experiments")
for year in ["1995", "2000", "2005"]

    # load data
    println("Year: $year")
    data_frame = CSV.read("data$year.csv", DataFrame)
    data_array = Array(data_frame[:,2:end])
    data = UpperTriangular(data_array + data_array') ./ 2
    n = size(data, 1)
    N = 0.5 * n * (n-1)
    nz = sum(data .== -Inf)
    nnz = N - nz
    println("Number of nodes n: $n")
    println("Number of total edges N: $N")
    println("Number of non-zero samples: $nnz")
    println("Proportion of non-zero samples: $(nnz/N)")

    # get bandwidth
    h_ROT = estimate_ROT_bandwidth(data, "epanechnikov_order_2")
    println("ROT bandwidth: $h_ROT")
    data[data .== -Inf] .= -1e10

    # fit dyadic kernel density estimator
    est = DyadicKernelDensityEstimator(
        kernel_name, h_ROT, significance_level,
        n_resample, sdp_solver, evals, data,
        Dict("year" => year))

    fit(est)

    # plot fhat
    fig, ax = plt.subplots(figsize=(4,4))
    ax.plot(est.evals, est.fhat,
        color = "black", linewidth=linewidth, linestyle=(0, (1,1)))

    # plot ucb
    ax.fill_between(est.evals, est.ucb[1,:], est.ucb[2,:],
        color="lightgray", linewidth=0.0)

    # plot pci
    plot_confidence_intervals(ax, est.evals, est.pci[1,:],
        est.pci[2,:], 10, "black", linewidth, "-")

    # save plot
    PyPlot.xlabel("Bilateral trade volume")
    plt.ylim(y_lim)
    x_min = minimum(est.evals)
    x_max = maximum(est.evals)
    plt.xlim((x_min, x_max))
    plt.yticks(range(0.0, stop=0.08, step=0.01))
    legend(handles=handles, loc="upper left")
    plt.ylabel("Density", labelpad=4.0)
    plt.tight_layout()
    PyPlot.savefig("trade_plot_$year.pdf")
    close("all")

end
