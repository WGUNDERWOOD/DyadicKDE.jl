using CSV
using DataFrames
using LinearAlgebra
using DyadicKDE

REPDIR = @__DIR__() * "/"
PLOTDIR = REPDIR * "plots/"
DATADIR = REPDIR * "data/"

include(REPDIR * "plot_helpers.jl")

# set plot parameters
linewidth = 1.0
y_lim = [-0.003, 0.11]
e = 0.22

# estimation parameters
n_evals = 100
kernel_name = "epanechnikov_order_4"
evals = collect(range(-10.0, stop=10.0, length=n_evals))
sdp_solver = "mosek"
n_resample = 10000
significance_level = 0.05

years = ["1995", "2000", "2005"]
year0 = "1995"

for year1 in years
    println("Year: ", year1)

    # get year abbreviations
    yr0 = year0[3:4]
    yr1 = year1[3:4]

    # read data
    W = DataFrame(CSV.File(DATADIR * "data_W_" * year1 * ".csv"))
    W = UpperTriangular(Array(W))
    h_ROT = estimate_ROT_bandwidth(W, "epanechnikov_order_2")
    println("h_ROT: ", h_ROT)
    W[W .== -Inf] .= -1e10
    X0 = DataFrame(CSV.File(DATADIR * "data_X_" * year0 * ".csv"))
    X0 = Array(X0.GDP_bracket)
    X1 = DataFrame(CSV.File(DATADIR * "data_X_" * year1 * ".csv"))
    X1 = Array(X1.GDP_bracket)

    # fit observed estimator
    est = DyadicKernelDensityEstimator(kernel_name, h_ROT, significance_level,
                                       n_resample, sdp_solver, evals,
                                       W, Dict())

    fit(est)

    if year0 == year1

        # legend observed
        handle_fhat = PyPlot.matplotlib.lines.Line2D([0], [0], color="k", lw=linewidth,
                                                     linestyle=(0, (1, 1)),
                                                     label="\$\\widehat f_W^{\\," * yr1 * "}(w)\$")
        handle_ucb = PyPlot.matplotlib.patches.Patch(facecolor="lightgray", edgecolor="lightgray",
                                                     label="UCB for \$f_W^{" * yr1 * "}(w)\$")
        handles = [handle_fhat, handle_ucb]

        # plot observed
        fig, ax = plt.subplots(figsize=(4, 4))
        ax.plot(est.evals, est.fhat,
                color="black", linewidth=linewidth, linestyle=(0, (1, 1)))
        ax.fill_between(est.evals, est.ucb[1, :], est.ucb[2, :],
                        color="lightgray", linewidth=0.0)

        PyPlot.xlabel("Bilateral trade volume")
        plt.ylim(y_lim)
        x_min = minimum(est.evals)
        x_max = maximum(est.evals)
        plt.xlim((x_min, x_max))
        plt.yticks(range(0.0, stop=0.1, step=0.02))
        legend(handles=handles, loc="upper left")
        plt.ylabel("Density", labelpad=4.0)
        plt.tight_layout()
        PyPlot.savefig(PLOTDIR * "trade_plot_" * year1 * ".pdf")
        close("all")

    else

        # fit estimator
        est_cf = CounterfactualDyadicKernelDensityEstimator(kernel_name, h_ROT, significance_level,
                                                            n_resample, sdp_solver, evals,
                                                            W, X0, X1,
                                                            Dict())

        fit(est_cf)

        # legend cf
        handle_fhat = PyPlot.matplotlib.lines.Line2D([0], [0], color="k", lw=linewidth,
                                                     linestyle=(0, (1, 1)),
                                                     label="\$\\widehat f_W^{\\," * yr1 * "}(w)\$")
        handle_ucb = PyPlot.matplotlib.patches.Patch(facecolor="lightgray", edgecolor="lightgray",
                                                     label="UCB for \$f_W^{" * yr1 * "}(w)\$")
        handle_fhat_cf = PyPlot.matplotlib.lines.Line2D([0], [0], color="royalblue", lw=linewidth,
                                                        linestyle=(0, (1, 1)),
                                                        label="\$\\widehat f_W^{\\," * yr1 *
                                                              "\\triangleright" * yr0 * "}(w)\$")
        handle_ucb_cf = PyPlot.matplotlib.patches.Patch(facecolor="cornflowerblue",
                                                        edgecolor="cornflowerblue", alpha=0.4,
                                                        label="UCB for \$f_W^{" * yr1 *
                                                              "\\triangleright" * yr0 * "}(w)\$")
        handles = [handle_fhat, handle_fhat_cf, handle_ucb, handle_ucb_cf]

        # plot counterfactual
        fig, ax = plt.subplots(figsize=(4, 4))

        # est
        ax.plot(est.evals, est.fhat,
                color="black", linewidth=linewidth, linestyle=(0, (1, 1)))
        ax.fill_between(est.evals, est.ucb[1, :], est.ucb[2, :],
                        color="lightgray", linewidth=0.0)

        # est_cf
        ax.plot(est_cf.evals, est_cf.fhat,
                color="royalblue", linewidth=linewidth, linestyle=(0, (1, 1)))
        ax.fill_between(est_cf.evals, est_cf.ucb[1, :], est_cf.ucb[2, :],
                        color="cornflowerblue", linewidth=0.0, alpha=0.3)

        PyPlot.xlabel("Bilateral trade volume")
        plt.ylim(y_lim)
        x_min = minimum(est_cf.evals)
        x_max = maximum(est_cf.evals)
        plt.xlim((x_min, x_max))
        plt.yticks(range(0.0, stop=0.1, step=0.02))
        legend(handles=handles, loc="upper left")
        plt.ylabel("Density", labelpad=4.0)
        plt.tight_layout()
        PyPlot.savefig(PLOTDIR * "trade_plot_" * year0 * "_" * year1 * ".pdf")
        close("all")
    end
end
