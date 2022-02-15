using PyPlot

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] = true
rcParams["mathtext.fontset"]  = "cm";
rcParams["font.family"]  = "serif";
rcParams["font.serif"]  = "cm";
plt.ioff()


function phi(t::Real)

    return (2 * pi)^(-0.5) * exp(-(t^2) / 2)
end


function generate_f(evals, p)

    f = (p[1]^2 + p[3]^2) * phi.(evals .- 1)
    f += p[2]*(2 - p[2]) * phi.(evals)
    f += 2*p[1]*p[3] * phi.(evals .+ 1)

    return f
end


function generate_Var_f_given_A(evals, p)

    Var_f_given_A = (p[1] + p[3])^2 * p[2] * phi.(evals).^2
    Var_f_given_A += p[3] * (p[3] * phi.(evals .- 1) + p[1] * phi.(evals .+ 1)).^2
    Var_f_given_A += p[1] * (p[1] * phi.(evals .- 1) + p[3] * phi.(evals .+ 1)).^2
    Var_f_given_A -= ((p[1]^2 + p[3]^2) * phi.(evals .- 1) +
        (p[1] + p[3]) * p[2] * phi.(evals) + 2 * p[1] * p[3] * phi.(evals .+ 1)).^2

    return Var_f_given_A

end


function plot_confidence_intervals(ax, x, y_lower, y_upper, num_intervals,
                                   line_color, line_width, line_style)

    inds = range(1, stop=length(x), length=2*num_intervals+1)
    inds = [inds[i] for i in 1:length(inds) if i%2 == 0]
    inds = Int.(round.(inds))
    @assert size(inds, 1) == num_intervals
    xr = x[inds]
    ylr = y_lower[inds]
    yur = y_upper[inds]

    ax.vlines(
        xr[1], ylr[1], yur[1], linewidth=line_width,
        linestyle=line_style, color=line_color)

    for i in 1:num_intervals

        ax.vlines(
            xr[i], ylr[i], yur[i], linewidth=line_width,
            linestyle=line_style, color=line_color)

        marker = PyPlot.matplotlib.path.Path([[-1,0], [1,0], [-1,0], [1,0]])

        ax.plot(
            xr[i], ylr[i], linewidth=0.0,
            color=line_color, marker=marker,
            markersize=4.0)

        ax.plot(
            xr[i], yur[i], linewidth=0.0,
            color=line_color, marker=marker,
            markersize=4.0)
    end
end
