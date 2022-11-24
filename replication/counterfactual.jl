using Revise
using DyadicKDE
using Distributions

n_data = 30
n_evals = 20
pW = [0.25, 0.0, 0.75]
W = make_data(n_data, pW)

pX0 = [0.2, 0.3, 0.5]
pX1 = [0.6, 0.2, 0.2]
X0 = rand(Categorical(pX0), n_data)
X1 = rand(Categorical(pX1), n_data)

kernel_name = "epanechnikov_order_2"
bandwidth = 1.0
significance_level = 0.05
n_resample = 10
sdp_solver = "mosek"
evals = collect(range(-2.0, stop=2.0, length=n_evals))
meta = Dict()

est = CounterfactualDyadicKernelDensityEstimator(
    kernel_name, bandwidth, significance_level,
    n_resample, sdp_solver, evals, W, X0, X1, meta)

DyadicKDE.fit(est)

display(est.fhat_cf)
