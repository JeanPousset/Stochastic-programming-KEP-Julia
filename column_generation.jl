include("out_ngb_G.jl")
include("DW_problem.jl")
include("pricing.jl")
include("initialisation.jl")

using Random


"""
    column_generation_ILP

Column generation method for the relaxed problem formulation

# Arguments
* `G::SimpleDiGraph` : sick/donnor pair compatibility graph
* `K::Int64` : maximum length of cycles
* `init_choice::String` : method for the first cycle choice
* `SP_method::String` : resolution method for the subproblem ("Bellmann" or "ILP")
* `SP_order::String` : order of subproblem resolution ("sequence", "base" or "random")
* `max_iter::Int64` : maximum number of iteration
* `verb::Int64` : verbosity level

# Returns 
* `Vector{Vector{Int64}}` : explored cycles
* `Int64` : objective value of the relaxed problem
"""
function column_generation_ILP(G::SimpleDiGraph, K::Int64; init_choice::String="half K=2", SP_method::String="ILP", SP_order::String="base", max_iter::Int64=500, verb::Int64=-1)

    I = Vector{Int64}(vertices(G)) # vertex indices (supposed to be ⟦1,length()⟧ without any index missed)

    ### Initialization

    # Dantizg-Wolfe formulation
    DW_dual, Π = initialize_DW_dual(I)

    # sub problem solver
    if SP_method == "ILP"

        # outneighbors graphs and association tables for the sub problems
        c_time = @elapsed begin
            Gsp_validities, Gs_prime, Φ = build_Gs_prime(G, K)
        end
        # lists invalid graphs:

        if verb >= 0
            println("-------------------")
            println("outneigbours graph built in $c_time seconds")
            println("\t -> invalids subproblems (G_o_prime without any path): $(I[.!Gsp_validities])")
            println("-------------------")
        end

        # ↓ warning ↓ : will create empty problem if it is infeasible (G_o_prime doesn't containt any path)
        SP = [initialize_SP_o(Gs_prime[o], Gsp_validities[o]) for o in I] # Sub Problems (SP_o) for o ∈ ⟦1,|I|⟧ 
        pricing = (Π_dual, sp_order) -> princing_ILP(SP, Gs_prime, Φ, Gsp_validities, Π_dual, sp_order, verb)

    elseif SP_method == "Bellmann"
        pricing = (Π_dual, sp_order) -> princing_Bellmann(G, K, Π_dual, sp_order, C_K_k, sp_order) # [TO DO]
    else
        @error "[column_generation]: unknown resolution method for the subproblems"
    end

    ## Choice of the first cycles
    if SP_method == "ILP" && init_choice == "path in G_o_prime"
        supp_C_k = first_G_op(Gs_prime, Gsp_validities, Φ)
    else
        supp_C_k = first_cycles_K(init_choice, G)    # C_K^(0)
    end
    C_K_k = copy(supp_C_k) # C_K^(k), initialized by C_K^(0)

    # subproblem resolution order : initializes with the correct ones only
    sp_order = I[Gsp_validities]
    k = 1

    ### Iterations
    while k < max_iter

        verb >= 0 && print("\r\t iter $k / $max_iter")

        ## DW restricted formulation resolution
        add_cycles_DW_dual!(DW_dual, Π, supp_C_k)
        Π_values = solve_DW_dual(DW_dual, verb)

        ## Sub Porblems (SP_o) resolutions

        # order update
        if SP_order == "random"
            shuffle!(sp_order)
        elseif SP_order == "sequence"   # begin with a different one each time
            circshift!(sp_order, -1)
        end

        # resolution
        supp_C_k, stop_algo = pricing(Π_values, sp_order)


        if stop_algo # It hasn't find any subcycles that produce z_o > 0 ⇒ we are at the optimum

            constraints = all_constraints(DW_dual, AffExpr, MOI.GreaterThan{Float64}) # Check if Affine is needed
            α = dual.(constraints)   # retrieve dual solution
            nb_transferts = objective_value(DW_dual)

            if verb >= 0
                println("============================================== ")
                println("| Column generation (linear) : optimum found |")
                println("\n\t - number of transfert perfomed : $nb_transferts")
                println("\t - transfert performed: $(C_K_k[α.>0.0])\n\n")
                verb >= 1 && @show C_K_k[α.>0.0]
            end

            return C_K_k, round(Int64, nb_transferts)
        end

        k += 1
        append!(C_K_k, supp_C_k) # save new cycles

    end

    @warn "non convergence for the column_generation (iteration $k/$max_iter)"

end

"""
    integer_solution

Solves the integer Dwantzig-Wolfe formulation with the given cycles (results of column_generation)

# Arguments
* `C_K::Vector{Vector{Int64}}` : know cycles (result of column_generation)
* `n_transfert_integer::Int64` : objective value of the relaxed problem
* `verb::Int64` : verbosity level

# Returns
* `Vector{Vector{Int64}}` : selected cycles
* `Int64` : objective value of the integer problem
"""
function integer_solution(C_K::Vector{Vector{Int64}}, n_transferts_relax::Int64; verb::Int64=-1)

    p = length(C_K)
    L = length.(C_K)

    C_K_i = [Int64[] for _ in 1:length(vertices(G))]

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
    integer_DW = Model(HiGHS.Optimizer)
    @variable(integer_DW, α[k in 1:p], Bin)
    @objective(integer_DW, Max, sum(α[k] * L[k] for k in 1:p))
    ## No more than one cycle per donnor
    @constraint(integer_DW, [i in vertices(G)], sum(α[k] for k in C_K_i[i]) <= 1)

    # Resolution
    set_silent(integer_DW)
    optimize!(integer_DW)

    if (termination_status(integer_DW) == MOI.OPTIMAL)
        n_transfert_integer = round(Int64, objective_value(integer_DW))

        if verb >= 0
            println("============================================== ")
            if (n_transferts_relax == n_transfert_integer)
                println("| Integer solution : optimum found |")
            else
                println("| Integer solution : lost from relaxation ($(n_transfert_integer-n_transferts_relax))")
            end
            println("\n\t - number of transfert perfomed : $(n_transfert_integer)")
            println("\t - transfert performed: $(C_K[value.(α) .> 0.0])\n\n")
        end

        return C_K[value.(α).>0], n_transfert_integer
    end
end