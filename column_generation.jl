include("brute_force_approach.jl") # to find the first <= K cycles
include("out_ngb_G.jl")
include("DW_problem.jl")
include("pricing.jl")


"""
    first_cycles_K

Compute the first cycles of length <= K for the column generation depending of the choice


# Arguments
* `G::SimpleDiGraph` : KEP graph
* `K::Int64` : maximum length of the transferts
* `init_choice::String` : choice of the kind of first cycles

# Returns
* `cycles_0::Vector{Vector{Int64}}` : first cycles for the column generation 
"""
function first_cycles_K(G::SimpleDiGraph, K::Int64, init_choice::String)::Vector{Vector{Int64}}

    if init_choice == "half K=2"
        cycles_2 = enumerate_cycles(G, 2)   # mask of half the indices
        selected_cycles = rand(1:length(cycles_2), div(length(cycles_2), 1))
        return cycles_2[selected_cycles]   # half of the cycles of length 2
    else
        @error("[first_cycles_K]: unknown choice : \"$init_choice\" for the first cycles")
    end
end


function column_generation_ILP(G::SimpleDiGraph, K::Int64; init_choice::String="half K=2", SP_method::String="ILP", max_iter::Int64=500)

    I = Vector{Int64}(vertices(G)) # vertex indices (supposed to be ⟦1,length()⟧ without any index missed)
    n = length(I) # number of donor/sick pairs
    k = 0 # iteration index

    ### Initialization

    # for o in I
    #     @show Φ[o]
    # end

    # Dantizg-Wolfe formulation
    DW_dual = initialize_DW_dual(I)
    # constraint_ref = ConstraintRef{Model,MOI.ConstraintIndex}[]
    supp_C_k = first_cycles_K(G, K, init_choice)    # C_K^(0)
    C_K_k = copy(supp_C_k) # C_K^(k), initialized by C_K^(0)


    # @show length(Gs_prime)
    # @show length(I)

    # sub problem solver
    if SP_method == "ILP"

        # outneighbors graphs and association tables for the sub problems
        timer = @timed begin
            Gsp_validities, Gs_prime, Φ = build_Gs_prime(G, K)
        end
        # lists invalid graphs:
        println("-------------------")
        println("outneigbours graph built in $timer")
        println("\t -> invalids subproblems (G_o_prime without any path): $(I[.!Gsp_validities])")
        println("-------------------")

        # ↓ warning ↓ : will create empty problem if it is infeasible (G_o_prime doesn't containt any path)
        SP = [initialize_SP_o(Gs_prime[o], Gsp_validities[o]) for o in I] # Sub Problems (SP_o) for o ∈ ⟦1,|I|⟧ 
        pricing = (Π_dual) -> princing_ILP(SP, Gs_prime, Φ, Gsp_validities, Π_dual)

    elseif SP_method == "Bellmann"
        pricing = (Π_dual) -> princing_Bellmann(Gs_prime, Φ, Π_dual) # [TO DO]
    else
        @error "[column_generation]: unknown resolution method for the subproblems"
    end

    k = 1

    ### Iterations
    while k < max_iter

        println("\t iter $k")

        ## DW restricted formulation resolution
        add_cycles_DW_dual!(DW_dual, supp_C_k)
        Π_dual = solve_DW_dual(DW_dual)

        ## Sub Porblems (SP_o) resolutions
        supp_C_k, stop_algo = pricing(Π_dual)


        if stop_algo # It hasn't find any subcycles that produce z_o > 0 ⇒ we are at the optimum

            constraints = all_constraints(DW_dual, AffExpr, MOI.GreaterThan{Float64}) # Check if Affine is needed
            α = dual.(constraints)   # retrieve dual solution

            println("============================================== ")
            println("| Column generation (linear) : optimum found |")
            println("\n\t - number of transfert perfomed : $(objective_value(DW_dual))")
            println("\t - transfert performed: $(C_K_k[α.>0.0])\n\n")
            @show α
            # @show C_K_k[α.>0.0]
            break
            # The dual of Π_dual (x) can be found in DW_dual attributes
            return C_K_k
        end

        append!(C_K_k, supp_C_k) # save new cycles
        k += 1
    end

    @warn "non convergence for the column_generation (iteration $k/$max_iter)"



end