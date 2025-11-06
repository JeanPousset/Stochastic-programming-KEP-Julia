"""
    initialize_SP_o

Initializes the o-th ILP fomrulation of the sub-problem (SP_o) without its objective

# Arguments
* `C_K_0::Vector{Vector{Int64}}`: first set of cycle shorter than K 
* `Gop_validity::Bool` : false if there is no cycles starting from o and passing by vertices > o (infeasible subproblem)

## Warning : will create empty problem if it is infeasible (G_o_prime doesn't containt any path)

# Returns
* `model::Model` : ILP formulation of the sub-problem without its objective
"""
function initialize_SP_o(G_o_prime::SimpleDiGraph, Gop_validity::Bool)::Model

    SP_o = Model(HiGHS.Optimizer)

    # build sub problem only if it is feasible
    if Gop_validity

        I_o = vertices(G_o_prime)

        #  Problem : maximum length for a flow of 1
        # → obectif (maximum length) to be defined with the obj_o_SP


        @variable(SP_o, x[i in I_o, j in outneighbors(G_o_prime, i)], Bin)
        # ↓ flow constraint ↓
        @constraint(SP_o, [i in I_o[2:end-1]], sum(x[i, j] for j in outneighbors(G_o_prime, i)) == sum(x[j, i] for j in inneighbors(G_o_prime, i)))
        # ↓ exactly 1 transfert per sick/donor pair ↓

        @constraint(SP_o, sum(x[1, j] for j in outneighbors(G_o_prime, 1)) == 1)
        # @constraint(SP_o, sum(x[j, 1] for j in outneighbors(G_o_prime, I_o[end])) == 1) # useless ?? (?2)
    end

    return SP_o
end

"""
    objective_SP_o!

Update the objective of the o-th sub-problem (SP_o) ILP formulation with the new Π_dual values. 
The weight of each edge is obtnaited with the index association table between G and G_o_prime : φ_o

# Arguments
* `SP_o::Model` : o-th sub-problem (SP_o) ILP formulation with old Π_dual values to be updated
* `G_o_prime::SimpleDiGraph` : outneighbors graph of (SP_o)
* `φ_o::Vector{Int64}` : index association table between G and G_o_prime
* `Π_dual::Vector{Float64}` : solution of the DW formulation problem restricted to C_K^(k)
"""
function objective_SP_o!(SP_o::Model, G_o_prime::SimpleDiGraph, φ_o::Vector{Int64}, Π_dual::Vector{Float64})
    @objective(SP_o, Max, sum(SP_o[:x][e.src, e.dst] * (1 - Π_dual[φ_o[e.dst]]) for e in edges(G_o_prime)))
end


"""
    SP_o_ILP_solver

Solves the ILP formulation of the o-th subproblem 

# Arguments
* `SP_o::Model` :  o-th subproblem
* `G_o_prime::SimpleDiGraph` : outneighbors of the o-th subproblem
* `verb::Int64=0` : verbosity (0: nothing, >= 1: print optimal value z_o, >= 2: solver logs)
"""
function solve_SP_o_ILP(SP_o::Model, G_o_prime::SimpleDiGraph, verb::Int64=0)::Tuple{Vector{Int64},Bool}

    set_optimizer_attribute(SP_o, "output_flag", verb >= 2) # activate or deactivates solver logs
    timer = @timed optimize!(SP_o)

    z_o = objective_value(SP_o)
    I_o = vertices(G_o_prime)

    if (termination_status(SP_o) == MOI.OPTIMAL)
        verb >= 1 && println("z_o := $z_o") # Display optimal value z_o (if verbosity >= 1)
    else
        @show SP_o
        @warn "[sub_pb]: subproblem o (ILP method) is not optimal (status : $(termination_status(SP_o)))"
    end

    ind_cycle = Vector{Int64}()  # indices in φ_o of the cycle vertices that gives z_o > 0 (! not the vertices of the G graph)
    stop_SP = false


    # Stopping condition for all o-th sub-problem
    if z_o > 0.  # we need to save the cycle that gives this z_o > o  
        stop_SP = true
        i = 1
        push!(ind_cycle, i)
        while i != length(I_o)  # stop when it reach the last vertex of G_o_prime : o
            for on in outneighbors(G_o_prime, i)
                if value(SP_o[:x][i, on]) > 0
                    i = on
                    push!(ind_cycle, i)
                    break
                end
            end
        end
    end
    # else : returns an empty cycle : but stopping flag is false

    return (ind_cycle, stop_SP)
end


"""
    princing_ILP

Solves the (SP_o) subproblems until it finds a cycle that gives z_o > 0

# Arguments
* `SP::Vector{Model}` : the o-th subproblems (SP_o) for o ∈ I
* `Gs_prime::Vector{SimpleDiGraph}` : outneighbors graphs for each subproblem
* `Φ::Vector{Vector{Int64}}` : assossiation tables for each outneigbours graph
* `Π_dual::Vector{Float64}` : solution of the DW formulation problem restricted to C_K^(k)
* `Gsp_validities::Vector{Bool}` : contains false there is no cycles starting from o and passing by vertices > o (infeasible subproblem)
"""
function princing_ILP(SP::Vector{Model}, Gs_prime::Vector{SimpleDiGraph}, Φ::Vector{Vector{Int64}}, Gsp_validities::Vector{Bool}, Π_dual::Vector{Float64})::Tuple{Vector{Vector{Int64}},Bool}


    cycles_k = Vector{Vector{Int64}}()

    for o in eachindex(SP) # o ∈ I

        # We solve the sub-problem only if it is valid
        if Gsp_validities[o]
            println("Begin pricing subproblem $o")

            objective_SP_o!(SP[o], Gs_prime[o], Φ[o], Π_dual) # update the o-th subproblem with the new values of Π_dual
            ind_cycle, flag_stop = solve_SP_o_ILP(SP[o], Gs_prime[o])

            if flag_stop # we found a path c so that z_o(c) > o

                path = Φ[o][ind_cycle[1:end-1]]  # retrieve the cycle vertex indices in G
                # subcycle segmentation
                i = 1  # vertex to check if it is twice in the path
                while i < length(path)
                    v_i = path[i]
                    j = i + 2
                    while j <= length(path)
                        v = path[j]
                        if v == v_i # v_i is twice in the path
                            sub_c = splice!(path, i+1:j) # cut the subcycle 
                            @show sub_c
                            if sum((1.0 .- Π_dual[sub_c])) >= 0 # Check that the cost is positive
                                push!(cycles_k, sub_c) # save the sub-cycle into the result
                            end
                        end
                        j += 1
                    end
                    i += 1
                end

                @show path
                push!(cycles_k, path)  # Save the last sub-cycle

                @show cycles_k

                # There is no need to look for other cycles that gives z_o ((1?) Is there really not any ?) 
                # column generation must keep going so we return false 

                println("end flag_stop")
                return (cycles_k, false)
            end
        end
    end

    # It hasn't find any subcycles that produce z_o > 0 ⇒ we are at the optimum
    return (cycles_k, true)
end