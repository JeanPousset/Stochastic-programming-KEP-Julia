"""
    enumerate_cycles

Enumerates all cycles of length < K in the given graph

# Arguments
* `G::SimpleDiGraph` : compatibility graph
* `K::Int64` : max length of searched cycles
"""
function enumerate_cycles(G::SimpleDiGraph, K::Int64)::Vector{Vector{Int64}}

    cycles = Vector{Vector{Int64}}()
    c = Vector{Int64}(undef, K) # current cycle so we don't reallocate memory each time we add a vertex
    I = sort(vertices(G))
    @inbounds for v in I
        c[1] = v
        search_cycles!(G, K, cycles, c, 1)
    end
    return cycles
end

# Recursive function to search cycle from a given path
function search_cycles!(G::SimpleDiGraph, K::Int64, cycles::Vector{Vector{Int64}}, c::Vector{Int64}, i::Int64)
    # neighborhood exploration
    for ngb in outneighbors(G, c[i])
        if ngb == c[1]
            push!(cycles, copy(c[1:i]))  # There is a cycle
        elseif i < K && ngb > c[1] && (K < 4 || !detect_subcycles(ngb, c, i))  # No sub-cycles for K < 4 [TO DO]
            c[i+1] = ngb
            search_cycles!(G, K, cycles, c, i + 1)
        end
    end
end

# Returns true if there is a subcycle in the given cycle c
function detect_subcycles(v::Int64, c::Vector{Int64}, i::Int64)::Bool
    for j in 2:i
        if c[j] == v
            return true
        end
    end
    return false
end

# Reshape a cycle so its first vertex is the "smallest"
# function normalize_cycle(cycle::Vector{Int64})
#     n = length(cycle)
#     min_index = argmin(cycle)
#     return vcat(cycle[min_index:end], cycle[1:min_index-1])
# end

"""
    brut_force

Solve a KEP problem with cycles ≤ K by listing all cycles of length ≤ K 
and solving the Dwantzig-Wolfe formulation on this cycle set. 

# Arguments
* `G::SimpleDiGraph` : compatibility graph
* `K::Int64` : max length of searched cycles
"""
function brut_force(G::SimpleDiGraph, K::Int64)
    C_K = enumerate_cycles(G, K)
    p = length(C_K)
    L = length.(C_K) # cycles lengths
    I = vertices(G) # vertex indices

    C_K_i = [Int64[] for _ in 1:length(I)]

    # Build cycle sets that contains the vertex i, ∀ i ∈ I
    #  -> we use sets to don't have duplicates
    for k in eachindex(C_K) # k ∈ ⟦1,p⟧
        print("$k \t")
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
    @objective(model, Max, sum(α[k] * L[k] for k in 1:p))
    ## No more than one cycle per donnor
    @constraint(model, [i in I], sum(α[k] for k in C_K_i[i]) <= 1)


    # Resolution
    set_silent(model)
    timer = @timed optimize!(model)

    # Display results
    println("\n-------------------------------")
    println("Brute force resolution of the KEP problem with K = $K \n")
    println("Resolution status : ", termination_status(model))
    println("computation time : $(timer.time)")
    if (termination_status(model) == MOI.OPTIMAL)
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
    end
    println("\n-------------------------------\n")
end