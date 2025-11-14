using Graphs

function Bellman_search(G::SimpleDiGraph, o::Int, K::Int, Π_dual::Vector{Float64})

    ##Initialisation
    pred = [Dict(v => -1 for v in vertices(G)) for k in 1:K]    #Vector of dictionaries for the predecessors of the vertice v at iteration k

    cost = [Dict(v => -Inf for v in vertices(G)) for k in 1:K]  #Vector of dictionaries for the cost associated with the vertice v at iteration k

    for k in 1:K
        pred[k][o] = o
        cost[k][o] = 0
    end

    is_in_L = falses(nv(G)) #Boolean list of every vertice which cost was updated after k-1 iterations
    is_in_L[o] = true

    ##Cycle search
    for k in 1:(K-1)
        is_in_L_tmp = falses(nv(G)) #Temporary boolean list of every vertice which cost was updated at iteration k
        for v in findall(is_in_L)
            for u in outneighbors(G, v)
                if u < o
                    continue
                end
                if cost[k][v] + (1 - Π_dual[u]) > cost[k+1][u] #Update the cost if it is better
                    cost[k+1][u] = cost[k][v] + (1 - Π_dual[u])
                    pred[k][u] = v
                    is_in_L_tmp[u] = true
                end
            end
        end
        is_in_L = is_in_L_tmp
    end

    ##Cycle reconstruction
    for k in K:-1:2
        #Search for the predecessor of o with the highest associated cost
        for k in K:-1:2
            for v in inneighbors(G, o)
                if cost[k][v] + (1 - Π_dual[o]) > cost[k][o]
                    cost[k][o] = cost[k][v] + (1 - Π_dual[o])
                    pred[k][o] = v
                end
            end
        end
        if pred[k][o] != o
            cycle = [pred[k][o]]
            current = pred[k][o]
            for i = (k-1):-1:1
                u = pred[i][current]
                if u == -1 #No predecessor of current at iteration i
                    continue
                end
                pushfirst!(cycle, u)
                current = u
            end
            #Subcycles, if a subcycle appears in cycle, we just return the subcycle
            if length(cycle) >= 3
                for i = 1:(length(cycle)-1)
                    for j = (i+1):length(cycle)
                        if cycle[i] == cycle[j]
                            cycle = cycle[i:(j-1)]
                            return cycle
                        end
                    end
                end
            end
            return cycle
        end
    end
    return nothing
end

"""
   pricing_Bellmann

Find the cycles of highest cost for every vertices of the graph

# Arguments
* `G::SimpleDiGraph` : KEP graph
* `K::Int64` : maximum length of the transferts
* `Π_dual::Vector{Float64}` : solution of the DW formulation problem restricted to C_K^(k)
* `C_K_k::Vector{Vector{Int}}` : cycles used to solve the master problem

# Returns
* `chosen_cycles::Vector{Int64}` : vector of every valid cycle found
* `true::Boolean` : if no cycle was found and thus we terminate the algorithm
* `false::Boolean` : if at least one cycle was found and thus we continue the algorithm

"""
function pricing_Bellmann(G::SimpleDiGraph, K::Int, Π_dual::Vector{Float64}, C_K_k::Vector{Vector{Int}}, sp_order::Vector{Int64})
    chosen_cycles = Vector{Vector{Int}}()
    for o in sp_order
        cycle = Bellman_search(G, o, K, Π_dual)
        if cycle === nothing || cycle in C_K_k
            continue
        end
        push!(chosen_cycles, cycle)
    end
    if isempty(chosen_cycles)
        return (chosen_cycles, true)
    end
    return (chosen_cycles, false)
end