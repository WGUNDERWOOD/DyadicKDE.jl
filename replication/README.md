# Replication files

The results of the empirical studies presented in
[Cattaneo, Feng and Underwood, 2022](https://arxiv.org/abs/2201.05967)
can be replicated by running the Julia scripts found in this directory.
Provided in this document are details of each script, along with a description
of the DOTS data used for the empirical studies.

## Julia scripts

The Julia scripts can be run by executing the shell script `RUNME.sh`.

### `make_distribution_plot.jl`

Generates the dyadic density plots which form Figure 1 in the main paper.

### `make_outcome_plot.jl`

Generates the typical outcomes plots which form Figure 2 in the main paper.

### `make_table.jl`

Generates the numerical results for Table 1 in the main paper.
This script may take several hours to run,
but can be accelerated for example by first starting Julia with
`julia -t 8` to use 8 CPU threads.
Exact results may vary due to Julia's
[pseudorandom number generation](https://docs.julialang.org/en/v1/stdlib/Random/).

### `make_trade_table.jl`

Generates the summary statistics for Table 2 in the main paper.

### `make_trade_plot.jl`

Generates the counterfactual plots for Figure 3 in the main paper
and Figure 1 in the supplemental appendix.

### `make_trade_histogram.jl`

Generates the GDP histogram plots for Figure 2 in the supplemental appendix.

### `make_trade_plot_parametric.jl`

Generates the parametric counterfactual plots for Figure 3 in the supplemental appendix.

## DOTS data

### Covariate data

Covariate data is included in `data_X_1995.csv`
and similarly for the other years.
The columns are as follow:

- `ISO`
  Three-letter code for each country

- `Population`
  The population of each country

- `Population_bracket`
  The decile of each country's population

- `GDP`
  The GDP of each country

- `GDP_bracket`
  The decile of each country's GDP

- `GDP_per_capita`
  The GDP per capita of each country

- `GDP_per_capita_bracket`
  The decile of each country's GDP per capita

- `GATT`
  Whether each country is in the General Agreement on Tariffs and Trade

### Trade data

Trade data is included in `data_W_1995.csv`
and similarly for the other years.
The data is organized as follows.

- The first row is a list of all the ISO country codes.
- The subsequent rows indicate the trade volume between countries,
  forming an upper triangular matrix with missing entries given as -Inf.
