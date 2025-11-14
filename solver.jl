include("column_generation.jl")

struct KEP_test
    G::SimpleDiGraph
    K::Int64
    SP_method::String
    SP_order::String
    max_iter::Int64
end