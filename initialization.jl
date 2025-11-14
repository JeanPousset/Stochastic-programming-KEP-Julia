include("brute_force_approach.jl") # to find the first <= K cycles


"""
    first_cycles_K

Compute the first cycles of length <= K for the column generation depending of the choice


# Arguments
* `init_choice::String` : choice of the kind of first cycles
* `G::SimpleDiGraph` : KEP graph
* `Gs_prime::Vector{SimpleDiGraph}` : outneighbors graphs for each subproblem
* `Gsp_validities::Vector{Bool}` : indicates if there is a path form o to o in G_o_prime
* `Φ::Vector{Vector{Int64}}` : assossiation tables for each outneigbours graph

# Returns
* `cycles_0::Vector{Vector{Int64}}` : first cycles for the column generation 
"""
function first_cycles_K(init_choice::String, G::SimpleDiGraph)::Vector{Vector{Int64}}

    if init_choice == "all K=2"
        cycles_2 = enumerate_cycles(G, 2)   # mask of half the indices
        return cycles_2 # every cycles of length 2
    elseif init_choice == "half K=2"
        cycles_2 = enumerate_cycles(G, 2)   # mask of half the indices
        selected_cycles = rand(1:length(cycles_2), div(length(cycles_2), 2))
        return cycles_2[selected_cycles]   # half of the cycles of length 2
    else
        @error("[first_cycles_K]: unknown choice : \"$init_choice\" for the first cycles")
    end
end

"""
    first_G_op

Gives the cycles from the first paths from o to o in the G_o_prime graphs

# Arguments
* `Gs_prime::Vector{SimpleDiGraph}` : outneighbors graphs for each subproblem
* `Gsp_validities::Vector{Bool}` : indicates if there is a path form o to o in G_o_prime
* `Φ::Vector{Vector{Int64}}` : assossiation tables for each outneigbours graph

# Returns
* `C_K_0::Vector{Vector{Int64}}` a fist set of cycle for the column generation
"""
function first_G_op(Gs_prime::Vector{SimpleDiGraph}, Gsp_validities::Vector{Bool}, Φ::Vector{Vector{Int64}})::Vector{Vector{Int64}}

    C_K_0 = Vector{Vector{Int64}}()

    for o in 1:length(Gs_prime)

        if Gsp_validities[o]

            ## retrieve path
            ind_path = Int64[]
            i = length(Φ[o])
            push!(ind_path, i)
            while i > 1
                i = first(inneighbors(Gs_prime[o], i)) # first inneigbour
                push!(ind_path, i)
            end
            outneighbors

            ## retrives cycles
            path = Φ[o][ind_path[1:end-1]]  # retrieve the cycle vertex indices in G

            ## subcycle segmentation
            i = 1  # vertex to check if it is twice in the path
            while i < length(path)
                v_i = path[i]
                j = i + 2
                while j <= length(path)
                    v = path[j]
                    if v == v_i # v_i is twice in the path
                        sub_c = splice!(path, i+1:j) # cut the subcycle
                        push!(C_K_0, sub_c)
                    end
                    j += 1
                end
                i += 1
            end
            push!(C_K_0, path)  # Save the last sub-cycle
        end
    end

    return C_K_0
end