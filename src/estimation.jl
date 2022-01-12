#using LinearAlgebra
#using Distributions
#using Statistics
#using Convex
#using Mosek
#using MosekTools

Base.@kwdef mutable struct DyadicKernelDensityEstimator

    # input parameters
    kernel_name::String
    bandwidth::Float64
    significance_level::Float64
    n_resample::Int
    evals::Vector{Float64}
    data::UpperTriangular{Float64}

    # final products
    fhat::Vector{Float64}
    Sigmahat::Symmetric{Float64}
    Sigmahatplus::Symmetric{Float64}
    ucb::Array{Float64, 2}
    pci::Array{Float64, 2}
    bci::Array{Float64, 2}

end



function DyadicKernelDensityEstimator(
    kernel_name::String,
    bandwidth::Float64,
    significance_level::Float64,
    n_resample::Int,
    evals::Vector{Float64},
    data::UpperTriangular{Float64})

    n_evals = length(evals)
    n_data = size(data, 1)

    estimator = DyadicKernelDensityEstimator(

        # input parameters
        kernel_name::String,
        bandwidth::Float64,
        significance_level::Float64,
        n_resample::Int,
        evals::Vector{Float64},
        data::UpperTriangular{Float64},

        # final products
        fill(NaN, (n_evals)),
        Symmetric(fill(NaN, (n_evals, n_evals))),
        Symmetric(fill(NaN, (n_evals, n_evals))),
        fill(NaN, (2, n_evals)),
        fill(NaN, (2, n_evals)),
        fill(NaN, (2, n_evals)),
    )

    return estimator

end


#=



function DyadicKernelDensityEstimator(params::Dict,
    data::DyadicData, evals::DyadicEvaluationPoints)

    estimator = DyadicKernelDensityEstimator(params)
    fit(params, estimator, data, evals)

    return estimator

end



function kernel(s::Float64, w::Float64, bandwidth::Float64, w_min::Float64,
                w_max::Float64, kernel_name::String)

    if abs(w - s) <= bandwidth
        if w_min <= w <= w_max
            if w_min <= s <= w_max

                at = min((w - w_min)/bandwidth, 1)
                bt = min((w_max - w)/bandwidth, 1)
                st = (s - w)/bandwidth

                if kernel_name == "epanechnikov_order_2"

                    if (at == 1.0) && (bt == 1.0)
                        return 0.75 * (1/bandwidth) * (1 - st^2)

                    else
                        c0 = -((48*(3*at^4-3*at^3*bt+at^2*(-5+3*bt^2)+bt^2*(-5+3*bt^2)+
                            at*(5*bt-3*bt^3)))/((at+bt)^3*(60+3*at^4-12*at^3*bt-44*bt^2+3*bt^4-
                            4*at*bt*(-8+3*bt^2)+at^2*(-44+30*bt^2))))
                        c1 = -((180*(at^3-at^2*bt+at*(-2+bt^2)-bt*(-2+bt^2)))/
                            ((at+bt)^3*(60+3*at^4-12*at^3*bt-44*bt^2+3*bt^4-4*at*bt*(-8+3*bt^2)+
                            at^2*(-44+30*bt^2))))
                        c2 = (48*(3*at^4-3*at^3*bt+at^2*(-5+3*bt^2)+bt^2*(-5+3*bt^2)+at*(5*bt-3*bt^3)))/
                            ((at+bt)^3*(60+3*at^4-12*at^3*bt-44*bt^2+3*bt^4-4*at*bt*(-8+3*bt^2)+
                            at^2*(-44+30*bt^2)))
                        c3 = (180*(at^3-at^2*bt+at*(-2+bt^2)-bt*(-2+bt^2)))/
                            ((at+bt)^3*(60+3*at^4-12*at^3*bt-44*bt^2+3*bt^4-4*at*bt*(-8+3*bt^2)+
                            at^2*(-44+30*bt^2)))

                        return (c0 + c1*st + c2*st^2 + c3*st^3) / bandwidth
                    end


                elseif kernel_name == "epanechnikov_order_4"

                    if (at == 1.0) && (bt == 1.0)
                        return (1/bandwidth) * ((45/32) - (75/16)*st^2 + (105/32)*st^4)

                    else
                        c0 = -((240*(25*at^12-225*at^11*bt+45*at^10*(-7+25*bt^2)+
                            at^9*(2835*bt-4125*bt^3)-30*at^7*bt*(252-1120*bt^2+705*bt^4)+
                            15*at^8*(56-945*bt^2+825*bt^4)+45*at^2*bt^4*(-588+840*bt^2-315*bt^4+25*bt^6)-
                            9*at*bt^5*(-588+840*bt^2-315*bt^4+25*bt^6)+bt^6*(-588+840*bt^2-315*bt^4+25*bt^6)-
                            15*at^3*bt^3*(-2548+4480*bt^2-2240*bt^4+275*bt^6)+
                            15*at^4*bt^2*(-1764+5460*bt^2-4270*bt^4+825*bt^6)+
                            14*at^6*(-42+2700*bt^2-4575*bt^4+1775*bt^6)-6*at^5*bt*(-882+11200*bt^2-13125*bt^4+3525*bt^6)))/
                            ((at+bt)^7*(8820+5*at^8-80*at^7*bt-13160*bt^2+5376*bt^4-540*bt^6+5*bt^8-80*at^5*bt*(-33+26*bt^2)+
                            20*at^6*(-27+34*bt^2)-16*at^3*bt*(651-690*bt^2+130*bt^4)+at^4*(5376-8940*bt^2+3130*bt^4)-
                            16*at*bt*(-560+651*bt^2-165*bt^4+5*bt^6)+4*at^2*(-3290+5334*bt^2-2235*bt^4+170*bt^6))))
                        c1 = -((2100*(at-bt)*(15*at^10-120*at^9*bt+at^8*(-238+555*bt^2)+at^7*(1904*bt-1920*bt^3)-
                            16*at^5*bt*(345-679*bt^2+285*bt^4)+at^6*(690-5656*bt^2+3930*bt^4)-
                            8*at*bt^3*(-504+690*bt^2-238*bt^4+15*bt^6)+bt^4*(-504+690*bt^2-238*bt^4+15*bt^6)-
                            16*at^3*bt*(-252+870*bt^2-679*bt^4+120*bt^6)+at^2*bt^2*(-8568+13290*bt^2-5656*bt^4+555*bt^6)+
                            2*at^4*(-252+6645*bt^2-7798*bt^4+1965*bt^6)))/
                            ((at+bt)^7*(8820+5*at^8-80*at^7*bt-13160*bt^2+5376*bt^4-540*bt^6+5*bt^8-80*at^5*bt*(-33+26*bt^2)+
                            20*at^6*(-27+34*bt^2)-16*at^3*bt*(651-690*bt^2+130*bt^4)+at^4*(5376-8940*bt^2+3130*bt^4)-
                            16*at*bt*(-560+651*bt^2-165*bt^4+5*bt^6)+4*at^2*(-3290+5334*bt^2-2235*bt^4+170*bt^6))))
                        c2 = (1200*(5*at^12-45*at^11*bt+15*at^10*(-7+15*bt^2)+at^9*(945*bt-825*bt^3)-
                            6*at^7*bt*(1386-2275*bt^2+705*bt^4)+3*at^8*(308-1575*bt^2+825*bt^4)+70*at^6*(-35+384*bt^2-375*bt^4+71*bt^6)-
                            30*at^5*bt*(-735+1778*bt^2-1071*bt^4+141*bt^6)-9*at*bt^3*(1764-2450*bt^2+924*bt^4-105*bt^6+5*bt^8)+
                            bt^4*(1764-2450*bt^2+924*bt^4-105*bt^6+5*bt^8)+15*at^2*bt^2*(1764-3234*bt^2+1792*bt^4-315*bt^6+15*bt^8)+
                            3*at^4*(588-16170*bt^2+22680*bt^4-8750*bt^6+825*bt^8)+at^3*(-15876*bt+59780*bt^3-53340*bt^5+13650*bt^7-825*bt^9)))/
                            ((at+bt)^7*(8820+5*at^8-80*at^7*bt-13160*bt^2+5376*bt^4-540*bt^6+5*bt^8-
                            80*at^5*bt*(-33+26*bt^2)+20*at^6*(-27+34*bt^2)-16*at^3*bt*(651-690*bt^2+130*bt^4)+
                            at^4*(5376-8940*bt^2+3130*bt^4)-16*at*bt*(-560+651*bt^2-165*bt^4+5*bt^6)+
                            4*at^2*(-3290+5334*bt^2-2235*bt^4+170*bt^6)))
                        c3 = (2100*(at-bt)*(15*at^10-120*at^9*bt-80*at^7*bt*(-25+24*bt^2)+
                            5*at^8*(-50+111*bt^2)-16*at^5*bt*(462-775*bt^2+285*bt^4)+
                            at^6*(924-6100*bt^2+3930*bt^4)-16*at^3*bt*(-630+1302*bt^2-775*bt^4+120*bt^6)+
                            2*at^4*(-630+8274*bt^2-8650*bt^4+1965*bt^6)-8*at*bt*(588-1260*bt^2+924*bt^4-
                            250*bt^6+15*bt^8)+bt^2*(588-1260*bt^2+924*bt^4-250*bt^6+15*bt^8)+
                            at^2*(588-12600*bt^2+16548*bt^4-6100*bt^6+555*bt^8)))/
                            ((at+bt)^7*(8820+5*at^8-80*at^7*bt-13160*bt^2+5376*bt^4-540*bt^6+5*bt^8-80*at^5*bt*(-33+26*bt^2)+
                            20*at^6*(-27+34*bt^2)-16*at^3*bt*(651-690*bt^2+130*bt^4)+at^4*(5376-8940*bt^2+3130*bt^4)-
                            16*at*bt*(-560+651*bt^2-165*bt^4+5*bt^6)+4*at^2*(-3290+5334*bt^2-2235*bt^4+170*bt^6)))
                        c4 = (3360*(15*at^10-135*at^9*bt+135*at^8*(-2+5*bt^2)-45*at^7*bt*(-54+55*bt^2)-3*at^5*bt*(2499-4750*bt^2+
                            1950*bt^4)+at^6*(833-6900*bt^2+4800*bt^4)-9*at*bt^3*(-630+833*bt^2-270*bt^4+15*bt^6)+
                            bt^4*(-630+833*bt^2-270*bt^4+15*bt^6)+15*at^2*bt^2*(-630+1029*bt^2-460*bt^4+45*bt^6)+
                            15*at^4*(-42+1029*bt^2-1230*bt^4+320*bt^6)-5*at^3*bt*(-1134+3724*bt^2-2850*bt^4+495*bt^6)))/
                            ((at+bt)^7*(8820+5*at^8-80*at^7*bt-13160*bt^2+5376*bt^4-540*bt^6+5*bt^8-80*at^5*bt*(-33+26*bt^2)+
                            20*at^6*(-27+34*bt^2)-16*at^3*bt*(651-690*bt^2+130*bt^4)+at^4*(5376-8940*bt^2+3130*bt^4)-
                            16*at*bt*(-560+651*bt^2-165*bt^4+5*bt^6)+4*at^2*(-3290+5334*bt^2-2235*bt^4+170*bt^6)))
                        c5 = (12600*(at-bt)*(2*at^8-16*at^7*bt-8*at^5*bt*(-39+32*bt^2)+at^6*(-39+74*bt^2)-
                            16*at^3*bt*(63-72*bt^2+16*bt^4)+at^4*(126-543*bt^2+284*bt^4)-8*at*bt*(-98+126*bt^2-39*bt^4+2*bt^6)+
                            bt^2*(-98+126*bt^2-39*bt^4+2*bt^6)+at^2*(-98+672*bt^2-543*bt^4+74*bt^6)))/
                            ((at+bt)^7*(8820+5*at^8-80*at^7*bt-13160*bt^2+5376*bt^4-540*bt^6+5*bt^8-80*at^5*bt*(-33+26*bt^2)+
                            20*at^6*(-27+34*bt^2)-16*at^3*bt*(651-690*bt^2+130*bt^4)+at^4*(5376-8940*bt^2+3130*bt^4)-
                            16*at*bt*(-560+651*bt^2-165*bt^4+5*bt^6)+4*at^2*(-3290+5334*bt^2-2235*bt^4+170*bt^6)))

                        return (c0 + c1*st + c2*st^2 + c3*st^3 + c4*st^4 + c5*st^5)/bandwidth
                    end


                else
                    error("unknown kernel")
                end

            end
        end
    end

    return 0.0
end



function estimate_fhat(
    estimator::DyadicKernelDensityEstimator,
    data::DyadicData,
    evals::DyadicEvaluationPoints)

    for k in 1:evals.n_evals

        estimator.fhat[k] = sum(kernel.(data.W_vec, evals.w[k], estimator.bandwidth,
                                        evals.w_min, evals.w_max, estimator.kernel_name))
        estimator.fhat[k] /= data.N_data
    end
end



function estimate_conditional_expectation(
    estimator::DyadicKernelDensityEstimator,
    data::DyadicData,
    evals::DyadicEvaluationPoints)

    cond_exp = zeros((evals.n_evals, data.n_data))

    for k in 1:evals.n_evals

        for i in 1:data.n_data
            for j in 1:data.n_data

                if i < j
                    if abs(data.W[i,j] - evals.w[k]) <= estimator.bandwidth
                        cond_exp[k,i] += kernel(data.W[i,j], evals.w[k],
                                                estimator.bandwidth, evals.w_min, evals.w_max,
                                                estimator.kernel_name)
                    end

                elseif j < i
                    if abs(data.W[j,i] - evals.w[k]) <= estimator.bandwidth
                        cond_exp[k,i] += kernel(data.W[j,i], evals.w[k],
                                                estimator.bandwidth, evals.w_min, evals.w_max,
                                                estimator.kernel_name)
                    end

                end
            end

        end

    end

    cond_exp /= (data.n_data - 1)

    return cond_exp
end



function estimate_Sigmahat(
    estimator::DyadicKernelDensityEstimator,
    data::DyadicData,
    evals::DyadicEvaluationPoints)

    n = data.n_data

    Sigmahat = zeros((evals.n_evals, evals.n_evals))
    term1 = zeros((evals.n_evals, evals.n_evals))
    term2 = zeros((evals.n_evals, evals.n_evals))
    term3 = zeros((evals.n_evals, evals.n_evals))

    cond_exp = estimate_conditional_expectation(estimator, data, evals)

    for w1 in 1:evals.n_evals
        for w2 in 1:evals.n_evals

            # Si Si' term
            for i in 1:n
                term1[w1,w2] += cond_exp[w1,i] * cond_exp[w2,i]
            end
            term1[w1,w2] *= (4 / (n^2))

            # kij kij' term
            w_w1 = evals.w[w1]
            w_w2 = evals.w[w2]

            if abs(w_w1 - w_w2) <= 2 * estimator.bandwidth
                for i in 1:data.N_data
                    term2[w1,w2] +=
                        kernel(data.W_vec[i], w_w1, estimator.bandwidth, evals.w_min,
                               evals.w_max, estimator.kernel_name) *
                        kernel(data.W_vec[i], w_w2, estimator.bandwidth, evals.w_min,
                               evals.w_max, estimator.kernel_name)
                end
            end

            term2[w1,w2] *= (4 / (n^2 * (n - 1)^2))

            # fhat fhat' term
            term3[w1,w2] += estimator.fhat[w1] * estimator.fhat[w2]
            term3[w1,w2] *= (4 * n - 6) / (n * (n - 1))

            # combine terms
            Sigmahat[w1,w2] = term1[w1,w2] - term2[w1,w2] - term3[w1,w2]

        end
    end

    estimator.Sigmahat = Symmetric(Sigmahat)
    estimator.min_eig_Sigmahat = minimum(eigen(estimator.Sigmahat).values)

end



function estimate_Sigmahatplus(
    params::Dict,
    estimator::DyadicKernelDensityEstimator,
    data::DyadicData,
    evals::DyadicEvaluationPoints)

    min_eig_Sigmahat = eigmin(estimator.Sigmahat)

    if min_eig_Sigmahat < 0
        if params["verbose"]
            println("Sigmahat not PSD...")
        end
    end

    if min_eig_Sigmahat >= 0
        estimator.Sigmahatplus = estimator.Sigmahat

    elseif estimator.psd_method == "sdp"
        d = diag(estimator.Sigmahat)
        sinv = 1 ./ sqrt.(d .+ d')
        C = Semidefinite(evals.n_evals)
        objective = maximum(abs(sinv .* (C - estimator.Sigmahat)))
        problem = minimize(objective)
        solve!(problem, Mosek.Optimizer, silent_solver=true)
        Sigmahatplus = Symmetric(evaluate(C))
        mineig_optsol = max(-eigmin(Sigmahatplus), 0)
        estimator.Sigmahatplus = Symmetric(Sigmahatplus + 2 * mineig_optsol * I)

        n = params["n_data"]
        h = params["bandwidth"]
        L = matrix_lipschitz_number(estimator.Sigmahatplus, evals.w)
        @assert L <= 4 / (n * h^3)

    else
        error("unknown psd_method")
    end


    # psd_error
    estimator.psd_error = maximum(abs.(estimator.Sigmahatplus - estimator.Sigmahat))

    # Sigmahatplus eigmin
    estimator.min_eig_Sigmahatplus = eigmin(estimator.Sigmahatplus)

    if params["verbose"]

        println("PSD maximum norm error:")
        println(maximum(abs.(estimator.Sigmahatplus - estimator.Sigmahat)))

        println("Minimum eigenvalue of Sigmahat:")
        println(minimum(eigen(estimator.Sigmahat).values))

        println("Minimum eigenvalue of Sigmahatplus:")
        println(minimum(eigen(estimator.Sigmahatplus).values))

    end


end



function estimate_uniform_confidence_band(
    estimator::DyadicKernelDensityEstimator,
    data::DyadicData,
    evals::DyadicEvaluationPoints)

    # make GP distribution and scaling
    gp_mu = zeros((evals.n_evals))
    gp = MvNormal(gp_mu, estimator.Sigmahatplus)
    Sigmahatplus_diag_sqrt = sqrt.(diag(estimator.Sigmahatplus))

    sup_scaled_gps = zeros((estimator.n_resample))

    for rep in 1:estimator.n_resample

        # sample and scale gp
        Z_hat = rand(gp)
        Z_hat_scaled = Z_hat ./ Sigmahatplus_diag_sqrt

        # get sup absolute value
        sup_scaled_gps[rep] = maximum(abs.(Z_hat_scaled))

    end

    # get quantile
    c = get_quantile(sup_scaled_gps, 1 - estimator.significance_level)
    ucb_width = c .* Sigmahatplus_diag_sqrt
    estimator.ucb[1,:] .= estimator.fhat .- ucb_width
    estimator.ucb[2,:] .= estimator.fhat .+ ucb_width

end



function estimate_pointwise_confidence_band(
    estimator::DyadicKernelDensityEstimator,
    data::DyadicData,
    evals::DyadicEvaluationPoints)

    for k in 1:evals.n_evals

        w = evals.w[k]
        variance = estimator.Sigmahatplus[k,k]
        normal = Normal(0.0, sqrt(variance))
        q = Statistics.quantile(normal, 1-(estimator.significance_level/2))

        estimator.pcb[1,k] = estimator.fhat[k] - q
        estimator.pcb[2,k] = estimator.fhat[k] + q

    end

end



function estimate_bonferroni_confidence_band(
    estimator::DyadicKernelDensityEstimator,
    data::DyadicData,
    evals::DyadicEvaluationPoints)

    bonferroni_significance_level = estimator.significance_level / evals.n_evals

    for k in 1:evals.n_evals

        w = evals.w[k]
        variance = estimator.Sigmahatplus[k,k]
        normal = Normal(0.0, sqrt(variance))
        q = Statistics.quantile(normal, 1-(bonferroni_significance_level/2))

        estimator.bcb[1,k] = estimator.fhat[k] - q
        estimator.bcb[2,k] = estimator.fhat[k] + q

    end

end



function estimate_epanechnikov_ROT_bandwidth(
    estimator::DyadicKernelDensityEstimator,
    data::DyadicData)

    # sigmahat_W
    mean_W = sum(data.W_vec) / data.N_data
    sum_of_squares_W = sum((data.W_vec .- mean_W).^2)
    sigmahat_W_squared = sum_of_squares_W / (data.N_data - 1)
    sigmahat_W = sqrt(sigmahat_W_squared)

    # IQRhat
    upper_quartile = get_quantile(data.W_vec, 0.75)
    lower_quartile = get_quantile(data.W_vec, 0.25)
    IQRhat = upper_quartile - lower_quartile

    # n_rate
    n_rate = data.N_data^(-1/5)

    estimator.bandwidth_ROT = 2.435 * min(sigmahat_W, IQRhat / 1.349) * n_rate

end



function fit(
    params::Dict,
    estimator::DyadicKernelDensityEstimator,
    data::DyadicData,
    evals::DyadicEvaluationPoints)

    if params["kernel_name"] == "epanechnikov_order_2"
        estimate_epanechnikov_ROT_bandwidth(estimator, data)

    elseif params["kernel_name"] == "epanechnikov_order_4"
        # TODO this is not ideal but works
        estimate_epanechnikov_ROT_bandwidth(estimator, data)
    end

    if params["bandwidth"] == "ROT"
        params["bandwidth"] = estimator.bandwidth_ROT
        estimator.bandwidth = estimator.bandwidth_ROT
    end

    estimate_fhat(estimator, data, evals)

    if params["fit_type"] == "partial"
        nothing

    elseif params["fit_type"] == "full"
        estimate_Sigmahat(estimator, data, evals)
        estimate_Sigmahatplus(params, estimator, data, evals)
        estimate_uniform_confidence_band(estimator, data, evals)
        estimate_pointwise_confidence_band(estimator, data, evals)
        estimate_bonferroni_confidence_band(estimator, data, evals)

    else
        error("invalid fit_type")

    end
end

=#
