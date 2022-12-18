"""
    make_dyadic_data(n_data::Int, p::Vector{Float64})

Generate some example dyadic data using a Gaussian mixture model.
"""
function make_dyadic_data(n_data::Int, p::Vector{Float64})

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



"""
    mean(x::Vector{<:Number})

Compute the mean of a vector of numbers.
"""
function mean(x::Vector{<:Number})

    return sum(x) / length(x)
end



"""
    phi(t::Real)

Compute the standard normal density function.
"""
function phi(t::Real)

    return (2 * pi)^(-0.5) * exp(-(t^2) / 2)
end



"""
    get_f(p::Vector{Float64}, evals::Vector{Float64})

Get the true density function from example dyadic Gaussian mixture data.
"""
function get_f(p::Vector{Float64}, evals::Vector{Float64})

    f = (p[1]^2 + p[3]^2) * phi.(evals .- 1)
    f += p[2]*(2 - p[2]) * phi.(evals)
    f += 2*p[1]*p[3] * phi.(evals .+ 1)
    return f
end



"""
    get_RIMSE(fhat::Vector{Float64}, f::Vector{Float64})

Compute the root integrated mean squared error of an estimate for a function.
"""
function get_RIMSE(fhat::Vector{Float64}, f::Vector{Float64})
    return sqrt(mean((fhat .- f).^2))
end



"""
    get_coverage(cb::Matrix{Float64}, f::Vector{Float64})

Check if a confidence band covers the true density function.
"""
function get_coverage(cb::Matrix{Float64}, f::Vector{Float64})
    return all(cb[1,:] .<= f .<= cb[2,:])
end



"""
    get_average_width(cb::Matrix{Float64})

Return the average width of a confidence band.
"""
function get_average_width(cb::Matrix{Float64})
    return mean(cb[2,:] .- cb[1,:])
end
