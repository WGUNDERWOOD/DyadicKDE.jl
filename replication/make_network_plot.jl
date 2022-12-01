using Revise
using CSV
using DataFrames
using Graphs
using GraphPlot
using Compose
using Cairo

# read data
year = "1995"
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

#draw(PNG("trade_network_1995.png", 40cm, 40cm), p)
draw(PDF("trade_network_1995.pdf", 40cm, 40cm), p)
