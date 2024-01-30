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
e = 0.22

# estimation parameters
n_evals = 10
kernel_name = "epanechnikov_order_4"
evals = collect(range(-10.0, stop=10.0, length=n_evals))
sdp_solver = "cosmo"
n_resample = 10000
significance_level = 0.05

years = ["1995", "2000", "2005"]

# histograms and parametric densities of GDP

for year in years
    println("Year: ", year)

    data = DataFrame(CSV.File(DATADIR * "data_X_" * year * ".csv"))
    gdp = Array(data.GDP)
    brackets = Array(data.GDP_bracket)
    n = length(brackets)
    X_levels = maximum(brackets)

    # histogram bins
    breaks = [minimum(gdp[i] for i in 1:n if brackets[i] == j) for j in 1:X_levels]
    push!(breaks, maximum(gdp))

    # parametric MLE fit log(gdp) ~ Normal
    mu = sum(log.(gdp) / n)
    sigma2 = sum(log.(gdp) .^ 2 / n) - mu^2
    x_evals = range(0, 20, step=0.1)
    para = (2 * pi * sigma2)^(-1 / 2) * exp.(-(x_evals .- mu) .^ 2 / (2 * sigma2))

    fig, ax = plt.subplots(figsize=(4, 4))
    ax.hist(log.(gdp), log.(breaks), density=true, color="lightgray", linewidth=0.5,
            edgecolor="black", label="Histogram")
    ax.plot(x_evals, para, color="black", linewidth=1, label="Normal")

    plt.legend()

    plt.ylim(0, 0.25)
    plt.xlim(0, 20)
    plt.xlabel("log GDP")
    plt.ylabel("Density", labelpad=4.0)
    plt.tight_layout()
    PyPlot.savefig(PLOTDIR * "trade_gdp_" * year * ".pdf")
    close("all")
end
