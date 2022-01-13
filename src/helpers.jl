function make_data(n_data::Int, p::Vector{Float64})

    c = cumsum(p)
    r = rand(n_data)
    data_A = -1.0 * (r .<= c[1]) + 1.0 * (r .>= c[2])
    data_V = randn((n_data, n_data))
    data = UpperTriangular(data_A .* data_A' .+ data_V)

    return data
end



function mean(x::Vector)

    return sum(x) / length(x)
end



function phi(t::Float64)

    return (2 * pi)^(-0.5) * exp(-(t^2) / 2)
end



function get_f(est::DyadicKernelDensityEstimator)

    p = est.meta["p"]
    f = (p[1]^2 + p[3]^2) * phi.(est.evals .- 1)
    f += p[2]*(2 - p[2]) * phi.(est.evals)
    f += 2*p[1]*p[3] * phi.(est.evals .+ 1)
    return f
end



function get_RIMSE(est::DyadicKernelDensityEstimator)

    f = get_f(est)
    RMSE = mean((est.fhat .- f).^2)
    return RMSE
end



function get_ucb_coverage(est::DyadicKernelDensityEstimator)

    return all(est.ucb[1,:] .<= get_f(est) .<= est.ucb[2,:])
end



function get_pci_coverage(est::DyadicKernelDensityEstimator)

    return all(est.pci[1,:] .<= get_f(est) .<= est.pci[2,:])
end



function get_ucb_average_width(est::DyadicKernelDensityEstimator)

    return sum(est.ucb[2,:] .- est.ucb[1,:]) / est.n_evals
end



function get_pci_average_width(est::DyadicKernelDensityEstimator)

    return sum(est.pci[2,:] .- est.pci[1,:]) / est.n_evals
end
