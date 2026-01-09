include("brute_force_approach.jl")

"""
    stoc_ME_brut_force

Solve a stochastic KEP problem with cycles ≤ K. It uses maximum expectation method by listing all cycles of length ≤ K, 
computing their expectations and solving the Dwantzig-Wolfe formulation on this cycle set. 

# Arguments
* `G::SimpleDiGraph` : compatibility graph
* `K::Int64` : max length of searched cycles
* `f_rates::Dict{Tuple{Int64,Int64},Float64}` : failure rates for each arc
* `v_map::Vector{Int64}` : mapping array in case of dealing with a subgraph of G

# Returns
* `Vector{Vector{Int64}}` : selected cycles
* `Float64` : objective value of the integer problem
"""
function stoc_brut_force(G::SimpleDiGraph, K::Int64, f_rates::Dict{Tuple{Int64,Int64},Float64}; v_map::Vector{Int64}=ones(Int64, nv(G)), verb::Int64=-1)
    C_K = enumerate_cycles(G, K)
    p = length(C_K)
    I = vertices(G) # vertex indices

    C_K_i = [Int64[] for _ in 1:length(I)]

    # Compute cycle expectations
    E_c = zeros(length(C_K))
    for k in eachindex(C_K)
        E_c[k] = cycle_expectation(v_map[C_K[k]], f_rates)
    end


    # Build cycle sets that contains the vertex i, ∀ i ∈ I
    #  -> we use sets to don't have duplicates
    for k in eachindex(C_K) # k ∈ ⟦1,p⟧
        for v in C_K[k]
            push!(C_K_i[v], k)
        end
    end

    # indices : 
    # - k ∈ ⟦1,p⟧, where p = |C_K| (k is equivalent to index c in the Dwantzig-Wolfe formulation)
    # - i,v ∈ I (vertices set)

    # Model
    model = Model(HiGHS.Optimizer)
    @variable(model, α[k in 1:p], Bin)
    @objective(model, Max, sum(α[k] * E_c[k] for k in 1:p))
    ## No more than one cycle per donnor
    @constraint(model, [i in I], sum(α[k] for k in C_K_i[i]) <= 1)


    # Resolution
    set_silent(model)
    timer = @timed optimize!(model)

    # Display results
    if (termination_status(model) != MOI.OPTIMAL)
        error("[stoc_brut_force] : non convergence of the ILP")
    end

    if verb >= 1
        println("\n-------------------------------")
        println("Brute force resolution of the KEP problem with K = $K \n")
        println("Resolution status : ", termination_status(model))
        println("computation time : $(timer.time)")
        println("Number of performed exchanges  : ", objective_value(model))
        println("List of performed exchanges :")

        for k in 1:p
            if value(α[k]) == 1 # cycle is choosen
                print("\n*** Cycle chosen: \n\t $(C_K[k][1])")
                for v in C_K[k][2:end]
                    print(" -> $v")
                end
            end
        end
        println("\n-------------------------------\n")
    end
    choosen_cycles = C_K[value.(α).>0]
    return choosen_cycles, sum(length.(choosen_cycles))
end