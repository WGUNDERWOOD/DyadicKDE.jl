using DyadicKDE

n_data = 10
p0 = [0.25, 0.0, 0.75]
p1 = [0.5, 0.0, 0.5]
W0 = make_data(n_data, p0)
W1 = make_data(n_data, p1)

X0 = rand(n_data) .<= 0.5
X1 = rand(n_data) .<= 0.3

# X is {0,1} valued
phat0 = Dict(
    "0" => sum(X0 .== 0) / n_data,
    "1" => sum(X0 .== 1) / n_data,
)

phat1 = Dict(
    "0" => sum(X1 .== 0) / n_data,
    "1" => sum(X1 .== 1) / n_data,
)

function psihat(x)
    return phat0[x] / phat1[x]
end

println(psihat("0"))
println(psihat("1"))
