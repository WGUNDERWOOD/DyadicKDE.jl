"Generate some example dyadic data using a Gaussian mixture model"
function make_data(n_data::Int, p::Vector{Float64})

    @assert n_data >= 2
    @assert sum(p) == 1.0
    @assert all(p .>= 0)
    @assert length(p) == 3

    c = cumsum(p)
    r = rand(n_data)
    data_A = -1.0 * (r .<= c[1]) + 1.0 * (r .>= c[2])
    data_V = randn((n_data, n_data))
    data = UpperTriangular(data_A .* data_A' .+ data_V)

    return data
end



"Compute the mean of a vector of numbers"
function mean(x::Vector{<:Number})

    return sum(x) / length(x)
end



"Compute the standard normal density function"
function phi(t::Real)

    return (2 * pi)^(-0.5) * exp(-(t^2) / 2)
end



"Get the true density function from example dyadic data"
function get_f(est::DyadicKernelDensityEstimator)

    p = est.meta["p"]
    f = (p[1]^2 + p[3]^2) * phi.(est.evals .- 1)
    f += p[2]*(2 - p[2]) * phi.(est.evals)
    f += 2*p[1]*p[3] * phi.(est.evals .+ 1)
    return f
end



"Compute the root integrated mean squared error"
function get_RIMSE(est::DyadicKernelDensityEstimator)

    f = get_f(est)
    RMSE = sqrt(mean((est.fhat .- f).^2))
    return RMSE
end



"Check if the uniform confidence band covers the true density function"
function get_ucb_coverage(est::DyadicKernelDensityEstimator)

    return all(est.ucb[1,:] .<= get_f(est) .<= est.ucb[2,:])
end



"Check if the pointwise confidence intervals all cover the true density function"
function get_pci_coverage(est::DyadicKernelDensityEstimator)

    return all(est.pci[1,:] .<= get_f(est) .<= est.pci[2,:])
end



"Return the average width of the uniform confidence band"
function get_ucb_average_width(est::DyadicKernelDensityEstimator)

    return sum(est.ucb[2,:] .- est.ucb[1,:]) / est.n_evals
end



"Return the average width of the pointwise confidence intervals"
function get_pci_average_width(est::DyadicKernelDensityEstimator)

    return sum(est.pci[2,:] .- est.pci[1,:]) / est.n_evals
end
