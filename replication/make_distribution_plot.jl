using DyadicKDE

REPDIR = @__DIR__() * "/"
PLOTDIR = REPDIR * "plots/"
DATADIR = REPDIR * "data/"

include(REPDIR * "plot_helpers.jl")

# specify parameters
n_evals = 100
degeneracies = ["total", "partial", "none"]
evals = collect(range(-2.0, stop=2.0, length=n_evals))
ps = Dict("total" => [0.5, 0.0, 0.5],
          "partial" => [0.25, 0.0, 0.75],
          "none" => [0.2, 0.2, 0.6])

println("Making plots")
linewidth = 1.0
y_lim = [-0.01, 0.45]
handle_f = PyPlot.matplotlib.lines.Line2D([0], [0], color="k", lw=linewidth, label="\$f_W(w)\$")
label = "\$\\mathrm{Var}[f_{W \\mid A}(w \\mid A_i)]^{1/2}\$"
handle_Var_f_given_A = PyPlot.matplotlib.lines.Line2D([0], [0], color="k", lw=linewidth,
                                                      label=label)
handles = [handle_f, handle_Var_f_given_A]

for degen in degeneracies
    p = ps[degen]
    f = generate_f(evals, p)
    Var_f_given_A = generate_Var_f_given_A(evals, p)

    # plot f
    fig, ax = plt.subplots(figsize=(3, 3))
    ax.plot(evals, f, color="black", linewidth=linewidth)

    # plot Var_f_given_A
    ax.plot(evals, Var_f_given_A .^ (0.5),
            color="black", linewidth=linewidth, linestyle=(0, (1, 1)))

    # save plot
    PyPlot.xlabel("Evaluation point, \$w\$", fontsize=12)
    plt.ylim(y_lim)
    x_min = minimum(evals)
    x_max = maximum(evals)
    plt.xlim((x_min, x_max))
    plt.yticks(range(0.0, stop=0.4, step=0.1), fontsize=11)
    plt.xticks(fontsize=11)
    legend(handles=handles, loc="upper left")
    plt.ylabel("Density", labelpad=4.0, fontsize=12)
    plt.tight_layout()
    PyPlot.savefig(PLOTDIR * "distribution_plot_$degen.pdf")
    close("all")
end
