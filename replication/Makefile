all:
	@mkdir -p plots data
	@julia --project=".." --color=yes make_distribution_plot.jl
	@julia --project=".." --color=yes make_outcome_plot.jl
	@julia --project=".." --color=yes make_trade_plot.jl
	@julia --project=".." --color=yes make_trade_histogram.jl
	@julia --project=".." --color=yes make_trade_plot_parametric.jl
	@julia --project=".." --color=yes logo.jl
	@julia --project=".." --color=yes make_table.jl
	@julia --project=".." --color=yes make_trade_table.jl
