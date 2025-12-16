"""
    cycle_expectancy

Computes a cycle expectation (without recourse)

# Arguments
* `cycle::Int64[]` : cycle vertices
* `f_rates::Dict{Tuple{Int64,Int64},Float64}` : failure rates for each arc

# Returns 
* `Float64` : cycle expectation
"""
function cycle_expectation(cycle::Int64[], f_rates::Dict{Tuple{Int64,Int64},Float64})::Float64
    E_c = length(cycle) # cycle expectation
    for i in 1:length(cycle)-1
        E_c *= (1 - get(ErrorException, f_rates, (cycle[i], cycle[i+1])))
    end
    E_C *= (1 - get(ErrorException, f_rates, (cycle[end], cycle[0])))
end

