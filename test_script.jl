# Include the necessary packages
using Random, Graphs, JuMP, HiGHS, DelimitedFiles, Distributions



include("KEP_readfile.jl")
include("column_generation.jl")


G, edge_weight = read_wmd_file("data_KEP/KEP_151.wmd");
K = 4
max_iter = 10000
transferts = column_generation_ILP(G, K; init_choice="all K=2", max_iter=max_iter)

println("=============[END]=============")
