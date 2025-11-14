include("solver.jl")
include("KEP_readfile.jl")



"""
    KEP_test

Data for a Kidney Exchange Problem and parameters to solve it
"""
struct KEP_test
    G::SimpleDiGraph
    K::Int64
    init_choice::String
    SP_method::String
    SP_order::String
    max_iter::Int64

    # constructor (with default parameters) -> load the graph_file
    function KEP_test(g_file::String; K::Int64=4, init_choice::String="half K=2", SP_method::String="ILP", SP_order::String="random", max_iter::Int64=500)
        G, _ = read_wmd_file(g_file)
        new(G, K, init_choice, SP_method, SP_order, max_iter)
    end
end

"""
    solve_KEP

Solves the relaxation of the KEP Dwantzig-Wolfe formulation with the column generation problem, then solve the integer KEP with the computed cycles
    
# Arguments
* `test::KEP_test` : a Kidney Exchange Problem instance
* `verb::Int64` : verbosity level
"""
function solve_KEP(test::KEP_test; verb::Int64=-1)

    C_K, n_transferts_relax = column_generation_ILP(test.G, test.K, test.init_choice, test.SP_method, test.SP_order, test.max_iter, verb)
    selected_cycles = integer_solution(C_K, n_transferts_relax, verb)

end
