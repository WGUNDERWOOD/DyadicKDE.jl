using Revise
using CSV
using DataFrames
using LinearAlgebra
using UnicodePlots
using DyadicKDE

include("./plot_helpers.jl")

# set plot parameters
linewidth = 1.0
y_lim = [-0.003, 0.08]
e = 0.22
pci_marker = PyPlot.matplotlib.path.Path(
    [[-1,0], [-1,e], [-1,-e], [-1,0], [1,0], [1,e], [1,-e], [1,0]])
handle_fhat = PyPlot.matplotlib.lines.Line2D(
    [0], [0], color="k", lw=linewidth, linestyle=(0, (1,1)),
    label="\$\\widehat f_W(w)\$")
handle_ucb = PyPlot.matplotlib.patches.Patch(
    facecolor="lightgray", edgecolor="lightgray", label="UCB")
handle_pci = PyPlot.matplotlib.lines.Line2D(
    [0], [0], color="black", lw=0, markersize=17,
    marker=pci_marker, label="PCI")
handles = [handle_fhat, handle_ucb, handle_pci]

# estimation parameters
n_evals = 50
kernel_name = "epanechnikov_order_4"
evals = collect(range(-10.0, stop=10.0, length=n_evals))
sdp_solver = "mosek"
n_resample = 10000
significance_level = 0.05

years = ["1995", "2000", "2005"]

for year0 in years
    for year1 in years

        println("$year0 ==> $year1")
        println()

        # read data
        data_W = DataFrame(CSV.File("data_W_" * year1 * ".csv"))
        data_W = UpperTriangular(Array(data_W))
        h_ROT = estimate_ROT_bandwidth(data_W, "epanechnikov_order_2")
        data_W[data_W .== -Inf] .= -1e10
        data_X0 = DataFrame(CSV.File("data_X_" * year0 * ".csv"))
        data_X0 = Array(data_X0.GDP_per_capita_bracket)
        data_X1 = DataFrame(CSV.File("data_X_" * year1 * ".csv"))
        data_X1 = Array(data_X1.GDP_per_capita_bracket)

        if year0 == year1

            # print metadata
            n = size(data_W, 1)
            N = 0.5 * n * (n-1)
            nz = sum(data_W .<= -1e9)
            nnz = N - nz
            println("Year: ", year1)
            println("Number of nodes n: $n")
            println("Number of total edges N: $N")
            println("Number of non-zero samples: $nnz")
            println("Proportion of non-zero samples: $(round(nnz/N, digits=3))")
            println()

            # fit estimator
            est = DyadicKernelDensityEstimator(
                kernel_name, h_ROT, significance_level,
                n_resample, sdp_solver, evals,
                data_W, Dict())

            fit(est)

            # plot original
            fig, ax = plt.subplots(figsize=(4,4))
            ax.plot(est.evals, est.fhat,
                    color = "black", linewidth=linewidth, linestyle=(0, (1,1)))
            ax.fill_between(est.evals, est.ucb[1,:], est.ucb[2,:],
                            color="lightgray", linewidth=0.0)
            plot_confidence_intervals(ax, est.evals, est.pci[1,:],
                                      est.pci[2,:], 10, "black", linewidth, "-")
            PyPlot.xlabel("Bilateral trade volume")
            plt.ylim(y_lim)
            x_min = minimum(est.evals)
            x_max = maximum(est.evals)
            plt.xlim((x_min, x_max))
            plt.yticks(range(0.0, stop=0.08, step=0.01))
            legend(handles=handles, loc="upper left")
            plt.ylabel("Density", labelpad=4.0)
            plt.tight_layout()
            PyPlot.savefig("trade_plot_" * year1 * ".pdf")
            PyPlot.savefig("trade_plot_" * year1 * ".png")
            close("all")

        else

            # fit estimator
            est_cf = CounterfactualDyadicKernelDensityEstimator(
                kernel_name, h_ROT, significance_level,
                n_resample, sdp_solver, evals,
                data_W, data_X0, data_X1,
                Dict())

            fit(est_cf)

            # plot counterfactual
            fig, ax = plt.subplots(figsize=(4,4))
            ax.plot(est_cf.evals, est_cf.fhat,
                    color = "black", linewidth=linewidth, linestyle=(0, (1,1)))
            ax.fill_between(est_cf.evals, est_cf.ucb[1,:], est_cf.ucb[2,:],
                            color="lightgray", linewidth=0.0)
            plot_confidence_intervals(ax, est_cf.evals, est_cf.pci[1,:],
                                      est_cf.pci[2,:], 10, "black", linewidth, "-")
            PyPlot.xlabel("Bilateral trade volume")
            plt.ylim(y_lim)
            x_min = minimum(est_cf.evals)
            x_max = maximum(est_cf.evals)
            plt.xlim((x_min, x_max))
            plt.yticks(range(0.0, stop=0.08, step=0.01))
            legend(handles=handles, loc="upper left")
            plt.ylabel("Density", labelpad=4.0)
            plt.tight_layout()
            PyPlot.savefig("trade_plot_" * year0 * "_" * year1 * ".pdf")
            PyPlot.savefig("trade_plot_" * year0 * "_" * year1 * ".png")
            close("all")
        end
    end
end
