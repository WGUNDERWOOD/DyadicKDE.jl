"""
Composite type to represent a dyadic kernel density estimator.
"""
Base.@kwdef mutable struct DyadicKernelDensityEstimator

    # input parameters
    kernel_name::String
    bandwidth::Float64
    significance_level::Float64
    n_resample::Int
    sdp_solver::String
    evals::Vector{Float64}
    W::UpperTriangular{Float64}
    meta::Dict

    # final products
    n_evals::Int
    evals_min::Float64
    evals_max::Float64
    n_data::Int
    N_data::Int
    W_vec::Vector{Float64}
    fhat::Vector{Float64}
    Sigmahat::Symmetric{Float64}
    Sigmahatplus::Symmetric{Float64}
    eigmin_Sigmahat::Float64
    eigmin_Sigmahatplus::Float64
    sdp_error::Float64
    ucb::Array{Float64,2}
    pci::Array{Float64,2}
    bci::Array{Float64,2}
end

"""
    DyadicKernelDensityEstimator(kernel_name, bandwidth, significance_level,
                                 n_resample, sdp_solver, evals, W, meta)

Construct a dyadic kernel density estimator.

# Arguments
- `kernel_name::String`: which kernel to use.
- `bandwidth::Float64`: the bandwidth for the estimator.
- `significance_level::Float64`: for the confidence band/intervals.
- `n_resample::Int`: the number of resamples used to construct the confidence band/intervals.
- `sdp_solver::String`: semi-definite program solver.
- `evals::Vector{Float64}`: points at which to evaluate the density estimator.
- `W::UpperTriangular{Float64}`: array of dyadic data.
- `meta::Dict`: any extra information to pass to the estimator.
"""
function DyadicKernelDensityEstimator(kernel_name::String,
                                      bandwidth::Float64,
                                      significance_level::Float64,
                                      n_resample::Int,
                                      sdp_solver::String,
                                      evals::Vector{Float64},
                                      W::UpperTriangular{Float64},
                                      meta::Dict)
    @assert bandwidth > 0
    @assert 0 <= significance_level <= 1
    @assert n_resample >= 1

    n_evals = length(evals)
    n_data = size(W, 1)
    N_data = Int(n_data * (n_data - 1) // 2)

    est = DyadicKernelDensityEstimator(
                                       # input parameters
                                       kernel_name,
                                       bandwidth,
                                       significance_level,
                                       n_resample,
                                       sdp_solver,
                                       evals,
                                       W,
                                       meta,

                                       # final products
                                       n_evals,
                                       minimum(evals),
                                       maximum(evals),
                                       n_data,
                                       N_data,
                                       fill(NaN, (N_data)),
                                       fill(NaN, (n_evals)),
                                       Symmetric(fill(NaN, (n_evals, n_evals))),
                                       Symmetric(fill(NaN, (n_evals, n_evals))),
                                       NaN,
                                       NaN,
                                       NaN,
                                       fill(NaN, (2, n_evals)),
                                       fill(NaN, (2, n_evals)),
                                       fill(NaN, (2, n_evals)))

    index = 1
    for i in 1:n_data
        for j in 1:n_data
            if i < j
                est.W_vec[index] = est.W[i, j]
                index += 1
            end
        end
    end

    return est
end

function estimate_fhat(est::DyadicKernelDensityEstimator)
    for k in 1:(est.n_evals)
        est.fhat[k] = sum(kernel.(est.W_vec, est.evals[k], est.bandwidth,
                                  est.evals_min, est.evals_max, est.kernel_name))

        est.fhat[k] /= est.N_data
    end
end

function estimate_conditional_expectation(est::DyadicKernelDensityEstimator)
    cond_exp = zeros((est.n_evals, est.n_data))

    for k in 1:(est.n_evals)
        for i in 1:(est.n_data)
            for j in 1:(est.n_data)
                if i < j
                    if abs(est.W[i, j] - est.evals[k]) <= est.bandwidth
                        cond_exp[k, i] += kernel(est.W[i, j], est.evals[k], est.bandwidth,
                                                 est.evals_min, est.evals_max, est.kernel_name)
                    end

                elseif j < i
                    if abs(est.W[j, i] - est.evals[k]) <= est.bandwidth
                        cond_exp[k, i] += kernel(est.W[j, i], est.evals[k], est.bandwidth,
                                                 est.evals_min, est.evals_max, est.kernel_name)
                    end
                end
            end
        end
    end

    cond_exp /= (est.n_data - 1)

    return cond_exp
end

function estimate_Sigmahat(est::DyadicKernelDensityEstimator)
    n = est.n_data

    Sigmahat = zeros((est.n_evals, est.n_evals))
    term1 = zeros((est.n_evals, est.n_evals))
    term2 = zeros((est.n_evals, est.n_evals))
    term3 = zeros((est.n_evals, est.n_evals))

    cond_exp = estimate_conditional_expectation(est)

    for w1 in 1:(est.n_evals)
        for w2 in 1:(est.n_evals)

            # Si Si' term
            for i in 1:n
                term1[w1, w2] += cond_exp[w1, i] * cond_exp[w2, i]
            end
            term1[w1, w2] *= (4 / (n^2))

            # kij kij' term
            w_w1 = est.evals[w1]
            w_w2 = est.evals[w2]

            if abs(w_w1 - w_w2) <= 2 * est.bandwidth
                for i in 1:(est.N_data)
                    s = est.W_vec[i]
                    term2[w1, w2] += kernel(s, w_w1, est.bandwidth, est.evals_min,
                                            est.evals_max, est.kernel_name) *
                                     kernel(s, w_w2, est.bandwidth, est.evals_min,
                                            est.evals_max, est.kernel_name)
                end
            end

            term2[w1, w2] *= (4 / (n^2 * (n - 1)^2))

            # fhat fhat' term
            term3[w1, w2] += est.fhat[w1] * est.fhat[w2]
            term3[w1, w2] *= (4 * n - 6) / (n * (n - 1))

            # combine terms
            Sigmahat[w1, w2] = term1[w1, w2] - term2[w1, w2] - term3[w1, w2]
        end
    end

    est.Sigmahat = Symmetric(Sigmahat)
    return est.eigmin_Sigmahat = eigmin(est.Sigmahat)
end

function matrix_lipschitz_number(mat::Symmetric{Float64}, v::Vector{Float64})
    @assert size(mat, 1) == length(v)
    steps = diff(mat, dims=1) ./ diff(v)
    return maximum(steps)
end

function estimate_Sigmahatplus(Sigmahat::Symmetric{Float64}, sdp_solver::String)
    eigmin_Sigmahat = eigmin(Sigmahat)

    if eigmin_Sigmahat >= 0
        # do nothing if already psd
        Sigmahatplus = Sigmahat

    else
        # formulate optimization problem
        d = diag(Sigmahat)
        sinv = 1 ./ sqrt.(d .+ d')
        C = Semidefinite(size(Sigmahat, 1))
        objective = maximum(abs(sinv .* (C - Sigmahat)))
        problem = minimize(objective)

        # solve optimization problem
        if sdp_solver == "mosek"
            solve!(problem, Mosek.Optimizer, silent_solver=true)
        elseif sdp_solver == "cosmo"
            solve!(problem, COSMO.Optimizer, silent_solver=true)
        else
            error("Unknown sdp_solver")
        end

        # get answer
        Sigmahatplus = Symmetric(evaluate(C))
        mineig_optsol = max(-eigmin(Sigmahatplus), 0)
        Sigmahatplus = Symmetric(Sigmahatplus + 2 * mineig_optsol * I)
    end

    eigmin_Sigmahatplus = eigmin(Sigmahatplus)
    sdp_error = maximum(abs.(Sigmahatplus - Sigmahat))
    @assert eigmin_Sigmahatplus >= 0

    return Dict("Sigmahatplus" => Sigmahatplus,
                "eigmin_Sigmahatplus" => eigmin_Sigmahatplus,
                "sdp_error" => sdp_error)
end

function get_quantile(A::Vector{<:Real}, q::Real)
    sorted_A = sort(A)
    idx = round(Int, length(A) * q)
    quantile = sorted_A[idx]

    return quantile
end

function estimate_uniform_confidence_band(Sigmahatplus::Symmetric{Float64}, n_resample::Int64,
                                          significance_level::Float64, fhat::Vector{Float64})

    # make GP distribution and scaling
    n_evals = length(fhat)
    gp_mu = zeros(n_evals)
    gp = MvNormal(gp_mu, Sigmahatplus)
    Sigmahatplus_diag_sqrt = sqrt.(diag(Sigmahatplus))

    sup_scaled_gps = zeros(n_resample)

    for rep in 1:n_resample

        # sample and scale gp
        Z_hat = rand(gp)
        Z_hat_scaled = Z_hat ./ Sigmahatplus_diag_sqrt

        # get sup absolute value
        sup_scaled_gps[rep] = maximum(abs.(Z_hat_scaled))
    end

    # get quantile
    c = get_quantile(sup_scaled_gps, 1 - significance_level)
    ucb_width = c .* Sigmahatplus_diag_sqrt
    ucb = fill(NaN, (2, n_evals))
    ucb[1, :] .= fhat .- ucb_width
    ucb[2, :] .= fhat .+ ucb_width

    return ucb
end

function estimate_pointwise_confidence_intervals(Sigmahatplus::Symmetric{Float64},
                                                 significance_level::Float64,
                                                 fhat::Vector{Float64})
    n_evals = length(fhat)
    pci = fill(NaN, (2, n_evals))

    for k in 1:n_evals
        variance = Sigmahatplus[k, k]
        normal = Normal(0.0, sqrt(variance))
        q = Statistics.quantile(normal, 1 - (significance_level / 2))

        pci[1, k] = fhat[k] - q
        pci[2, k] = fhat[k] + q
    end

    return pci
end

function estimate_bonferroni_confidence_intervals(Sigmahatplus::Symmetric{Float64},
                                                  significance_level::Float64,
                                                  fhat::Vector{Float64})
    n_evals = length(fhat)
    bci = fill(NaN, (2, n_evals))
    bonferroni_significance_level = significance_level / n_evals

    for k in 1:n_evals
        variance = Sigmahatplus[k, k]
        normal = Normal(0.0, sqrt(variance))
        q = Statistics.quantile(normal, 1 - (bonferroni_significance_level / 2))

        bci[1, k] = fhat[k] - q
        bci[2, k] = fhat[k] + q
    end

    return bci
end

"""
    estimate_ROT_bandwidth(W::UpperTriangular{Float64},
                           kernel_name::String)

Estimate a rule-of-thumb bandwidth from dyadic data.
"""
function estimate_ROT_bandwidth(W::UpperTriangular{Float64}, kernel_name::String)

    # get data_vec
    n_data = size(W, 1)
    N_data = Int(n_data * (n_data - 1) // 2)
    W_vec = Float64[]

    for i in 1:n_data
        for j in 1:n_data
            if i < j
                if -Inf < W[i, j] < Inf
                    push!(W_vec, W[i, j])
                end
            end
        end
    end

    if kernel_name == "epanechnikov_order_2"

        # sigmahat_W
        mean_W = mean(W_vec)
        sum_of_squares_W = sum((W_vec .- mean_W) .^ 2)
        sigmahat_W_squared = sum_of_squares_W / (length(W_vec) - 1)
        sigmahat_W = sqrt(sigmahat_W_squared)

        # IQRhat
        upper_quartile = get_quantile(W_vec, 0.75)
        lower_quartile = get_quantile(W_vec, 0.25)
        IQRhat = upper_quartile - lower_quartile

        # n_rate
        n_rate = N_data^(-1 / 5)

        return 2.435 * min(sigmahat_W, IQRhat / 1.349) * n_rate

    else
        error("Unknown kernel_name")
    end
end

"""
    fit(est::DyadicKernelDensityEstimator)

Fit a dyadic kernel density estimator to data.
"""
function fit(est::DyadicKernelDensityEstimator)
    estimate_fhat(est)
    estimate_Sigmahat(est)
    SDP = estimate_Sigmahatplus(est.Sigmahat, est.sdp_solver)
    est.Sigmahatplus = SDP["Sigmahatplus"]
    est.eigmin_Sigmahatplus = SDP["eigmin_Sigmahatplus"]
    est.sdp_error = SDP["sdp_error"]
    est.ucb = estimate_uniform_confidence_band(est.Sigmahatplus, est.n_resample,
                                               est.significance_level, est.fhat)
    est.pci = estimate_pointwise_confidence_intervals(est.Sigmahatplus, est.significance_level,
                                                      est.fhat)
    return est.bci = estimate_bonferroni_confidence_intervals(est.Sigmahatplus,
                                                              est.significance_level, est.fhat)
end

"""
    display(est::DyadicKernelDensityEstimator)

Display a dyadic kernel density estimator.
"""
function Base.display(est::DyadicKernelDensityEstimator)
    println("DyadicKernelDensityEstimator")
    println("----------------------------")

    println("Kernel: ", est.kernel_name)
    println("Bandwidth: ", est.bandwidth)
    println("Significance level: ", est.significance_level)
    println("Num. GP resamples: ", est.n_resample)
    println("SDP solver: ", est.sdp_solver)
    println("Num. data nodes (n): ", est.n_data)
    println("Num. data points (N): ", est.N_data)
    println("Num. evaluation points: ", est.n_evals)
    println("Min. evaluation point: ", est.evals_min)
    println("Max. evaluation point: ", est.evals_max)
    println("Min. eigval of Sigmahat: ", est.eigmin_Sigmahat)
    println("Min. eigval of Sigmahatplus: ", est.eigmin_Sigmahatplus)
    return println("SDP error: ", est.sdp_error)
end
