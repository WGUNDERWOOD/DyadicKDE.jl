using CSV
using DataFrames
using DyadicKDE

data_W_1995 = DataFrame(CSV.File("data_W_1995.csv"))
data_W_2000 = DataFrame(CSV.File("data_W_2000.csv"))
data_W_2005 = DataFrame(CSV.File("data_W_2005.csv"))

data_X_1995 = DataFrame(CSV.File("data_X_1995.csv"))
data_X_2000 = DataFrame(CSV.File("data_X_2000.csv"))
data_X_2005 = DataFrame(CSV.File("data_X_2005.csv"))

# get bandwidth
h_ROT = estimate_ROT_bandwidth(data, "epanechnikov_order_2")
println("ROT bandwidth: $h_ROT")
data[data .== -Inf] .= -1e10

# fit dyadic kernel density estimator
est = DyadicKernelDensityEstimator(
    kernel_name, h_ROT, significance_level,
    n_resample, sdp_solver, evals, data,
    Dict("year" => year))
