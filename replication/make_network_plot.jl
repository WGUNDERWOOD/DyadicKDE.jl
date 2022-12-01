using Revise
using CSV
using DataFrames
using Graphs
using GraphPlot
using Compose
using Cairo

years = ["1995", "2000", "2005"]

for year in years

    # read data
    data_W = DataFrame(CSV.File("data_W_" * year * ".csv"))
    data_W = exp.(Array(data_W))
    data_X = DataFrame(CSV.File("data_X_" * year * ".csv"))

    g = Graph(data_W + data_W')

    nodefillc = "black"
    nodesize = data_X.GDP .^ 0.2
    edgelinewidth = [data_W[src(e),dst(e)] for e in edges(g)] .^ 0.5
    layout = (args...) -> spring_layout(args...; C=40)
    p = gplot(g, layout=layout, nodesize=nodesize,
              edgelinewidth=edgelinewidth, nodefillc=nodefillc)

    n = length(vertices(g))
    N = 0.5 * n * (n-1)
    println("Year: ", year)
    println("Nodes: ", n)
    println("Edges: ", length(edges(g)))
    println("Density: ", density(g))
    println("Clustering coefficient: ", global_clustering_coefficient(g))
    println()

    #draw(PNG("trade_network_1995.png", 40cm, 40cm), p)
    draw(PDF("trade_network_" * year* ".pdf", 40cm, 40cm), p)
end
