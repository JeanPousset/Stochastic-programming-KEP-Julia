# Include the necessary packages
using Random, Graphs, JuMP, HiGHS, DelimitedFiles, Distributions



include("KEP_readfile.jl")
include("column_generation.jl")


G, edge_weight = read_wmd_file("data_KEP/KEP_071.wmd");
K = 4
max_iter = 10000
verbosity = 0
SP_method = "Belmann"
SP_order = "random"
C_K, n_transferts_relax = column_generation_ILP(G, K; init_choice="all K=2", SP_method=SP_method, SP_order=SP_order, max_iter=max_iter, verb=verbosity)
selected_cycles = integer_solution(C_K, n_transferts_relax, verb=verbosity)

println("=============[END]=============")
