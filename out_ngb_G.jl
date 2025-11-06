"""
    build_Gs_prime

Builds the |I| sub-graphs of outneighbors G_o_prime and the assocation table φ_o (o ∈ I)

# Arguments
* `G::SimpleDiGraph` : graph of sick/donnor pairs compatibilities
* `K::Int64` : max number of cycles allowed

# Returns
* `Gs_prime` : vector of the outneighbors sub-graphs G_o_prime (o ∈ I)
* `Φ` : vector the the association table φ_o (o ∈ I)
* `Gsp_validities::Vector{Bool}` : contains false if there is no cycles starting from o and passing by vertices > o (infeasible subproblem)
"""
function build_Gs_prime(G::SimpleDiGraph, K::Int64)::Tuple{Vector{Bool},Vector{SimpleDiGraph},Vector{Vector{Int64}}}

    I = vertices(G) # vertex indices (supposed to be ⟦1,length()⟧ without any index missed)
    # results allocation
    Gsp_validities = [false for _ in 1:length(I)] # is_valid[o] == false if there is no cycles starting from o and passing by vertices > o
    Gs_prime = [SimpleDiGraph() for _ in 1:length(I)]
    Φ = [Int64[] for _ in 1:length(I)] # φ ∈ Φ is th table of association between G_o_prime vertices (indices) and G vertices (values)


    # ** used indices: **
    # - o ∈ ⟦1,|I|⟧ : index for sick/donnor pair
    # - j ∈ ⟦2,K-1⟧ : index for the set of outneighbors we are building
    # - i : next free index for the outneighbors graph
    # - on : outneighbor of the current vertex of I_j
    # - l : index of the vertices of the previous outneighbor set I_j-1
    # - m : index of the vertices of the current outneighbor I_j that is being built

    for o in I

        # Initialization
        G_o_prime = SimpleDiGraph()
        add_vertex!(G_o_prime)
        i = 2 # next free index for the outneighbors graph (1 is for o)
        φ_o = Vector{Int64}() # stores corresponding index in G for the j-th outneighbors set
        first_φ_j = 2 # index in φ of the first element of the j-th outneighbors
        push!(φ_o, o) # first vertex is o

        # I_1 : outneighbors of o 
        for on in outneighbors(G, o) # 'on' : outneighbor
            if on > o # o must be the smallest index of the cycle
                add_vertex!(G_o_prime)
                if !add_edge!(G_o_prime, 1, i)
                    @warn "1 : can't add edge (already in)"
                end
                i += 1
                push!(φ_o, on)
            end
        end

        last_φ_j = i - 1 # index in φ of the last element of the j-th outneighbors

        ## I_j : outneighbors of I_j-1 without o
        for j in 2:K-2
            for l in first_φ_j:last_φ_j
                v = φ_o[l]
                for on in outneighbors(G, v)
                    if on >= o # o must be the smallest index of the cycle
                        flag_in = false
                        for m in last_φ_j+1:i-1
                            if on == φ_o[m] # the outneighbor 'on' is already in I_j
                                success = add_edge!(G_o_prime, l, m)
                                if !success
                                    @warn "k : can't add edge (already in)"
                                end
                                flag_in = true
                                break
                            end
                        end
                        if !flag_in # the outneighbor 'on' is not in I_j 
                            add_vertex!(G_o_prime)
                            success = add_edge!(G_o_prime, l, i)
                            if !success
                                @warn "k : can't add edge (not in)"
                            end
                            i += 1
                            push!(φ_o, on)
                        end
                    end
                end
            end
            # update cursors
            first_φ_j = last_φ_j + 1
            last_φ_j = i - 1
        end

        ## I_K-1
        I_K1_empty = true   # flag to know if the K-1 set of outneigbours is empty (i.e. there is no cycle)
        in_ngbrs_o = inneighbors(G, o) # avoid computing it multiple times
        for l in first_φ_j:last_φ_j
            v = φ_o[l]
            # we want the vertex that are outneighbor of I_K-2 and inneighbor of o
            for on in intersect(outneighbors(G, v), in_ngbrs_o)
                if on >= o # o must be the smallest index of the cycle
                    flag_in = false
                    for m in last_φ_j+1:i-1
                        if on == φ_o[m] # on (the outneighbor) is already in I_K-1
                            if !add_edge!(G_o_prime, l, m)
                                @warn "K-1 : can't add edge (already in)"
                            end
                            flag_in = true
                            break
                        end
                    end
                    if !flag_in # the outneighbor on is not already in I_K-1
                        I_K1_empty = false # there is a cycle -> we can save this G_o_prime
                        add_vertex!(G_o_prime)
                        if !add_edge!(G_o_prime, l, i)
                            @warn "K-1 : can't add edge (not already in)"
                        end
                        i += 1
                        push!(φ_o, on)
                    end
                end
            end
        end


        ## I_K : {o}

        # If we don't add any vertex in this outneigbours, there is no cycle from o with only vertices > o 
        if I_K1_empty
            # @warn "[build_GS_prime]: G_$(o)_prime doesn't containt a path ⇒ o-th subproblem is not feasible"
            Gsp_validities[o] = false
            # we don't save G_o_prime and φ[o] (initialized empty at the beginning in the result vectors)
        else
            Gsp_validities[o] = true

            # update cursors
            first_φ_j = last_φ_j + 1
            last_φ_j = i - 1

            # Add last arcs (to I_K = {o}) where o = i 
            add_vertex!(G_o_prime) # for o_out
            for l in first_φ_j:last_φ_j
                !add_edge!(G_o_prime, l, i) && @warn "K : can't add edge"
            end
            push!(φ_o, o) # here at the end i = |φ_o|

            # Save outneighbors graph G_o_prime and association table in result
            Gs_prime[o] = copy(G_o_prime)
            Φ[o] = φ_o
        end
        # print("[$o]")
        # display_G_o_prime(G_o_prime, φ_o)
    end

    return (Gsp_validities, Gs_prime, Φ)
end

"""
    display_G_o_prime

displays a outneigbours graph G_o_prime with the vertex of the initial G

# Arguments
* `G_o_prime::SimpleDiGraph` : outneighbors of the o-th subproblem
* `φ_o::Vector{Int64}` : index association table between G and G_o_prime
"""
function display_G_o_prime(G_o_prime::SimpleDiGraph, φ_o::Vector{Int64})
    println("---------- Display G_o_prime ----------")
    @show φ_o
    for e in edges(G_o_prime)
        println("$(e.src) -> $(e.dst)")
    end
end