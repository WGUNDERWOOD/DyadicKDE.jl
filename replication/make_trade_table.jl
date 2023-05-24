# generate the trade network summary plots found in
# https://arxiv.org/abs/2201.05967

using CSV
using LinearAlgebra
using DataFrames

years = ["1995", "2000", "2005"]

for year in years

    W = DataFrame(CSV.File(@__DIR__() * "/data_W_" * year * ".csv"))
    W = UpperTriangular(Array(W))
    n = size(W, 1)
    edges = Int(sum(W + W' .> -Inf) / 2)
    edge_density = 2 * edges / n / (n - 1)
    degrees = [sum(W[i, j] + W[j, i] > -Inf for j in 1:n) for i in 1:n]
    av_degree = sum(degrees) / n
    triangles = 0

    for i in 1:n
        for j in i:n
            if W[i, j] > -Inf
                for k in j:n
                    if W[j, k] > -Inf && W[i, k] > -Inf
                        triangles += 1
                    end
                end
            end
        end
    end

    clustering_coeff = 6 * triangles / sum(degrees .* (degrees .- 1))

    println("Year: ", year)
    println("Nodes: ", n)
    println("Edges: ", edges)
    println("Edge density: ", edge_density)
    println("Average degree: ", av_degree)
    println("Clustering coefficient: ", clustering_coeff)
    println()

end
