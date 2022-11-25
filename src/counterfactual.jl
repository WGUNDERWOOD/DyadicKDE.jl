"""
Composite type to represent a counterfactual dyadic kernel density estimator.
"""
Base.@kwdef mutable struct CounterfactualDyadicKernelDensityEstimator

    # input parameters
    kernel_name::String
    bandwidth::Float64
    significance_level::Float64
    n_resample::Int
    sdp_solver::String
    evals::Vector{Float64}
    W::UpperTriangular{Float64}
    X0::Vector{Int}
    X1::Vector{Int}
    meta::Dict

    # final products
    n_evals::Int
    evals_min::Float64
    evals_max::Float64
    n_data::Int
    N_data::Int
    X_levels::Int
    W_vec::Vector{Float64}
    fhat_cf::Vector{Float64}
    Sigmahat_cf::Symmetric{Float64}
    Sigmahatplus_cf::Symmetric{Float64}
    eigmin_Sigmahat_cf::Float64
    eigmin_Sigmahatplus_cf::Float64
    sdp_error::Float64
    ucb::Array{Float64, 2}
    pci::Array{Float64, 2}
    bci::Array{Float64, 2}

end



"""
    CounterfactualDyadicKernelDensityEstimator(kernel_name, bandwidth, significance_level,
                                               n_resample, sdp_solver, evals,
                                               W, X0, X1, meta)

Construct a counterfactual dyadic kernel density estimator.

# Arguments
- `kernel_name::String`: which kernel to use.
- `bandwidth::Float64`: the bandwidth for the estimator.
- `significance_level::Float64`: for the confidence band/intervals.
- `n_resample::Int`: the number of resamples used to construct the confidence band/intervals.
- `sdp_solver::String`: semi-definite program solver.
- `evals::Vector{Float64}`: points at which to evaluate the density estimator.
- `W::UpperTriangular{Float64}`: array of treated dyadic data.
- `X0::Vector{Int}`: categorical vector of untreated covariates.
- `X1::Vector{Int}`: categorical vector of treated covariates.
- `meta::Dict`: any extra information to pass to the estimator.
"""
function CounterfactualDyadicKernelDensityEstimator(
    kernel_name::String,
    bandwidth::Float64,
    significance_level::Float64,
    n_resample::Int,
    sdp_solver::String,
    evals::Vector{Float64},
    W::UpperTriangular{Float64},
    X0::Vector{Int},
    X1::Vector{Int},
    meta::Dict)

    @assert bandwidth > 0
    @assert 0 <= significance_level <= 1
    @assert n_resample >= 1

    n_evals = length(evals)
    n_data = size(W, 1)
    N_data = Int(n_data * (n_data - 1) // 2)
    X_levels = max(maximum(X0), maximum(X1))

    @assert length(X0) == n_data
    @assert length(X1) == n_data
    @assert sort(unique(X0)) == collect(1:X_levels)
    @assert sort(unique(X1)) == collect(1:X_levels)

    est = CounterfactualDyadicKernelDensityEstimator(

        # input parameters
        kernel_name,
        bandwidth,
        significance_level,
        n_resample,
        sdp_solver,
        evals,
        W,
        X0,
        X1,
        meta,

        # final products
        n_evals,
        minimum(evals),
        maximum(evals),
        n_data,
        N_data,
        X_levels,
        fill(NaN, (N_data)),
        fill(NaN, (n_evals)),
        Symmetric(fill(NaN, (n_evals, n_evals))),
        Symmetric(fill(NaN, (n_evals, n_evals))),
        NaN,
        NaN,
        NaN,
        fill(NaN, (2, n_evals)),
        fill(NaN, (2, n_evals)),
        fill(NaN, (2, n_evals)),
    )

    index = 1
    for i in 1:n_data
        for j in 1:n_data
            if i < j
                est.W_vec[index] = est.W[i,j]
                index += 1
            end
        end
    end

    return est

end



function estimate_phat(est::CounterfactualDyadicKernelDensityEstimator,
                       treatment::Int64)

    if treatment == 0
        X = est.X0
    elseif treatment == 1
        X = est.X1
    else
        error("invalid treatment")
    end

    phat = fill(0.0, est.X_levels)

    for k in 1:est.X_levels
        phat[k] = sum(X .== k) / est.n_data
    end

    return phat
end



function estimate_psihat(est::CounterfactualDyadicKernelDensityEstimator)

    phat0 = estimate_phat(est, 0)
    phat1 = estimate_phat(est, 1)
    psihat_values = phat0 ./ phat1
    psihat = [psihat_values[x] for x in est.X1]

end



function estimate_kappahat(est::CounterfactualDyadicKernelDensityEstimator)

    phat0 = estimate_phat(est, 0)
    phat1 = estimate_phat(est, 1)

    term1 = zeros((est.n_data, est.n_data))
    term2 = zeros((est.n_data, est.n_data))

    for r in 1:est.n_data
        for i in 1:est.n_data
            term1[r,i] = ((est.X0[r] == est.X1[i]) - phat0[est.X1[i]]) / phat1[est.X1[i]]
            term2[r,i] = (phat0[est.X1[i]] / phat1[est.X1[i]]) *
                ((est.X1[r] == est.X1[i]) - phat1[est.X1[i]]) / phat1[est.X1[i]]
        end
    end

    kappahat = term1 - term2

    return kappahat
end



function estimate_fhat_cf(est::CounterfactualDyadicKernelDensityEstimator)

    est.fhat_cf .= 0
    psihat = estimate_psihat(est)

    for k in 1:est.n_evals
        for i in 1:est.n_evals
            for j in 1:est.n_evals
                if i < j
                    est.fhat_cf[k] +=
                        kernel(est.W[i,j], est.evals[k], est.bandwidth,
                               est.evals_min, est.evals_max,
                               est.kernel_name) .*
                        psihat[i] * psihat[i]
                end
            end
        end

        est.fhat_cf[k] /= est.N_data
    end

end


function estimate_conditional_expectation(est::CounterfactualDyadicKernelDensityEstimator)

    S = zeros((est.n_evals, est.n_data))
    Stilde = zeros((est.n_evals, est.n_data))
    cond_exp = zeros((est.n_evals, est.n_data))

    psihat = estimate_psihat(est)
    kappahat = estimate_kappahat(est)

    # make S
    for k in 1:est.n_evals

        for i in 1:est.n_data
            for j in 1:est.n_data

                if i < j
                    if abs(est.W[i,j] - est.evals[k]) <= est.bandwidth
                        S[k,i] += kernel(
                            est.W[i,j], est.evals[k], est.bandwidth,
                            est.evals_min, est.evals_max, est.kernel_name
                        ) * psihat[j]
                    end

                elseif j < i
                    if abs(est.W[j,i] - est.evals[k]) <= est.bandwidth
                        S[k,i] += kernel(
                            est.W[j,i], est.evals[k], est.bandwidth,
                            est.evals_min, est.evals_max, est.kernel_name
                        ) * psihat[j]
                    end
                end
            end
        end
    end

    S /= (est.n_data - 1)

    # make Stilde
    for k in 1:est.n_evals

        for i in 1:est.n_data
            for j in 1:est.n_data
                for r in 1:est.n_data
                    if (j < r) && (i != j) && (i != r)

                        if abs(est.W[r,j] - est.evals[k]) <= est.bandwidth
                            Stilde[k,i] += kernel(
                                est.W[r,j], est.evals[k], est.bandwidth,
                                est.evals_min, est.evals_max, est.kernel_name
                            ) * psihat[r] * kappahat[i,j]
                        end
                    end
                end
            end
        end
    end

    Stilde *= 2 / ((est.n_data - 1) * (est.n_data - 2))

    # make cond_exp
    for k in 1:est.n_evals
        for i in 1:est.n_data
            cond_exp[k,i] = psihat[i] .* S[k,i] + Stilde[k,i]
        end
    end

    return cond_exp
end



function estimate_Sigmahat_cf(est::CounterfactualDyadicKernelDensityEstimator)

    n = est.n_data

    Sigmahat_cf = zeros((est.n_evals, est.n_evals))
    term1 = zeros((est.n_evals, est.n_evals))
    term2 = zeros((est.n_evals, est.n_evals))
    term3 = zeros((est.n_evals, est.n_evals))

    cond_exp = estimate_conditional_expectation(est)
    psihat = estimate_psihat(est)

    for w1 in 1:est.n_evals
        for w2 in 1:est.n_evals

            # conditional expectation term
            for i in 1:n
                term1[w1,w2] += cond_exp[w1,i] * cond_exp[w2,i]
            end
            term1[w1,w2] *= (4 / (n^2))

            # kij kij' psii^2 psij^2 term
            w_w1 = est.evals[w1]
            w_w2 = est.evals[w2]

            if abs(w_w1 - w_w2) <= 2 * est.bandwidth
                for i in 1:est.n_data
                    for j in 1:est.n_data
                        if i < j
                        term2[w1,w2] +=
                            kernel(
                                est.W[i,j], w_w1, est.bandwidth, est.evals_min,
                                est.evals_max, est.kernel_name
                            ) * kernel(
                                est.W[i,j], w_w2, est.bandwidth, est.evals_min,
                                       est.evals_max, est.kernel_name
                            ) * psihat[i]^2 * psihat[j]^2
                        end
                    end
                end
            end

            term2[w1,w2] *= (4 / (n^3 * (n - 1)))

            # fhat_cf fhat_cf' term
            term3[w1,w2] += est.fhat_cf[w1] * est.fhat_cf[w2]
            term3[w1,w2] *= (4 / n)

            # combine terms
            # TODO should this be (-term2)?
            Sigmahat_cf[w1,w2] = term1[w1,w2] + term2[w1,w2] - term3[w1,w2]

        end
    end

    est.Sigmahat_cf = Symmetric(Sigmahat_cf)
    est.eigmin_Sigmahat_cf = eigmin(est.Sigmahat_cf)

end


function estimate_Sigmahatplus_cf(est::CounterfactualDyadicKernelDensityEstimator)

    # TODO reduce code reuse here and later

    if est.eigmin_Sigmahat_cf >= 0
        # do nothing if already psd
        est.Sigmahatplus_cf = est.Sigmahat_cf

    else
        # formulate optimization problem
        d = diag(est.Sigmahat_cf)
        sinv = 1 ./ sqrt.(d .+ d')
        C = Semidefinite(est.n_evals)
        objective = maximum(abs(sinv .* (C - est.Sigmahat_cf)))
        problem = minimize(objective)

        # solve optimization problem
        if est.sdp_solver == "mosek"
            solve!(problem, Mosek.Optimizer, silent_solver=true)
        elseif est.sdp_solver == "cosmo"
            solve!(problem, COSMO.Optimizer, silent_solver=true)
        else
            error("Unknown sdp_solver")
        end

        # get answer
        Sigmahatplus_cf = Symmetric(evaluate(C))
        mineig_optsol = max(-eigmin(Sigmahatplus_cf), 0)
        est.Sigmahatplus_cf = Symmetric(Sigmahatplus_cf + 2 * mineig_optsol * I)

        # check lipschitz property
        n = est.n_data
        h = est.bandwidth
        L = matrix_lipschitz_number(est.Sigmahatplus_cf, est.evals)
        @assert L <= 4 / (n * h^3)

    end

    est.eigmin_Sigmahatplus_cf = eigmin(est.Sigmahatplus_cf)
    est.sdp_error = maximum(abs.(est.Sigmahatplus_cf - est.Sigmahat_cf))
    @assert est.eigmin_Sigmahatplus_cf >= 0

end

#=



function get_quantile(A::Vector{<:Real}, q::Real)

    sorted_A = sort(A)
    idx = round(Int, length(A) * q)
    quantile = sorted_A[idx]

    return quantile
end



function estimate_uniform_confidence_band(est::DyadicKernelDensityEstimator)

    # make GP distribution and scaling
    gp_mu = zeros((est.n_evals))
    gp = MvNormal(gp_mu, est.Sigmahatplus)
    Sigmahatplus_diag_sqrt = sqrt.(diag(est.Sigmahatplus))

    sup_scaled_gps = zeros((est.n_resample))

    for rep in 1:est.n_resample

        # sample and scale gp
        Z_hat = rand(gp)
        Z_hat_scaled = Z_hat ./ Sigmahatplus_diag_sqrt

        # get sup absolute value
        sup_scaled_gps[rep] = maximum(abs.(Z_hat_scaled))

    end

    # get quantile
    c = get_quantile(sup_scaled_gps, 1-est.significance_level)
    ucb_width = c .* Sigmahatplus_diag_sqrt
    est.ucb[1,:] .= est.fhat .- ucb_width
    est.ucb[2,:] .= est.fhat .+ ucb_width

end



function estimate_pointwise_confidence_intervals(est::DyadicKernelDensityEstimator)

    for k in 1:est.n_evals

        w = est.evals[k]
        variance = est.Sigmahatplus[k,k]
        normal = Normal(0.0, sqrt(variance))
        q = Statistics.quantile(normal, 1-(est.significance_level/2))

        est.pci[1,k] = est.fhat[k] - q
        est.pci[2,k] = est.fhat[k] + q

    end

end



function estimate_bonferroni_confidence_intervals(est::DyadicKernelDensityEstimator)

    bonferroni_significance_level = est.significance_level / est.n_evals

    for k in 1:est.n_evals

        w = est.evals[k]
        variance = est.Sigmahatplus[k,k]
        normal = Normal(0.0, sqrt(variance))
        q = Statistics.quantile(normal, 1-(bonferroni_significance_level/2))

        est.bci[1,k] = est.fhat[k] - q
        est.bci[2,k] = est.fhat[k] + q

    end

end



"""
    estimate_ROT_bandwidth(data::UpperTriangular{Float64},
                           kernel_name::String)

Estimate a rule-of-thumb bandwidth from dyadic data.
"""
function estimate_ROT_bandwidth(data::UpperTriangular{Float64},
                                kernel_name::String)

    # get data_vec
    n_data = size(data, 1)
    N_data = Int(n_data * (n_data - 1) // 2)
    data_vec = Float64[]

    for i in 1:n_data
        for j in 1:n_data
            if i < j
                if -Inf < data[i,j] < Inf
                    push!(data_vec, data[i,j])
                end
            end
        end
    end


    if kernel_name == "epanechnikov_order_2"

        # sigmahat_W
        mean_W = sum(data_vec) / length(data_vec)
        sum_of_squares_W = sum((data_vec .- mean_W).^2)
        sigmahat_W_squared = sum_of_squares_W / (length(data_vec) - 1)
        sigmahat_W = sqrt(sigmahat_W_squared)

        # IQRhat
        upper_quartile = get_quantile(data_vec, 0.75)
        lower_quartile = get_quantile(data_vec, 0.25)
        IQRhat = upper_quartile - lower_quartile

        # n_rate
        n_rate = N_data^(-1/5)

        return 2.435 * min(sigmahat_W, IQRhat / 1.349) * n_rate

    else
        error("Unknown kernel_name")
    end

end


=#

"""
    fit(est::CounterfactualDyadicKernelDensityEstimator)

Fit a counterfactual dyadic kernel density estimator to data.
"""
function fit(est::CounterfactualDyadicKernelDensityEstimator)

    estimate_fhat_cf(est)
    estimate_Sigmahat_cf(est)
    estimate_Sigmahatplus_cf(est)
    #estimate_uniform_confidence_band(est)
    #estimate_pointwise_confidence_intervals(est)
    #estimate_bonferroni_confidence_intervals(est)

end




"""
    display(est::CounterfactualDyadicKernelDensityEstimator)

Display a counterfactual dyadic kernel density estimator.
"""
function Base.display(est::CounterfactualDyadicKernelDensityEstimator)

    println("CounterfactualDyadicKernelDensityEstimator")
    println("------------------------------------------")

    println("Kernel: ", est.kernel_name)
    println("Bandwidth: ", est.bandwidth)
    println("Significance level: ", est.significance_level)
    println("Num. GP resamples: ", est.n_resample)
    println("SDP solver: ", est.sdp_solver)
    println("Num. data nodes (n): ", est.n_data)
    println("Num. data points (N): ", est.N_data)
    println("Num. covariate levels: ", est.X_levels)
    println("Num. evaluation points: ", est.n_evals)
    println("Min. evaluation point: ", est.evals_min)
    println("Max. evaluation point: ", est.evals_max)
    println("Min. eigval of Sigmahat: ", est.eigmin_Sigmahat_cf)
    println("Min. eigval of Sigmahatplus: ", est.eigmin_Sigmahatplus_cf)
    println("SDP error: ", est.sdp_error)

end
