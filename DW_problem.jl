"""
    initialize_DW_dual

Initializes the dual of the Dantzig-Wolfe decomposition without cycles (i.e. without constraint)

# Arguments
* `I::Vector{Int64}` : vertex (pairs) indices

# Returns
* `model::Model` : dual problem of the Dantzig-Wolfe decomposition without it's constraints
* `Π` : dictionary of problem variables
"""
function initialize_DW_dual(I::Vector{Int64})

    DW_dual = Model(HiGHS.Optimizer)
    @variable(DW_dual, Π[i in I] >= 0.)
    @objective(DW_dual, Min, sum(Π[i] for i in I))
    return DW_dual, Π
end

"""
    add_cycles_DW_dual!

Adds constraints for cycles added to the dual problem of the DW formulation

# Arguments
* `DW_dual::Model` : dual problem of the DW formulation restricted to C_K^(k) to be update_DW_dual
* `Π` : dictionary of problem variables
* `supp_C_k::Vector{Vector{Int64}}` : cycles to add
"""
function add_cycles_DW_dual!(DW_dual::Model, Π, supp_C_k::Vector{Vector{Int64}})
    L = length.(supp_C_k) # cycles' length
    # ↓ dual constraint ↓
    for i in eachindex(supp_C_k)
        c = supp_C_k[i]
        @constraint(DW_dual, sum(Π[v] for v in c) >= L[i])
    end
end

"""
    solve_DW_dual

Solves the dual of the DW formulation restricted to C_K^(k)

# Arguments
* `DW_dual::Model` :  dual problem of the DW formulation restricted to C_K^(k)
* `verb::Int64=0` : verbosity (0: nothing, >= 1: print solution Π, >= 2: solver logs)
"""
function solve_DW_dual(DW_dual::Model, verb::Int64=0)::Vector{Float64}

    set_optimizer_attribute(DW_dual, "output_flag", verb >= 2) # activate or deactivates solver logs
    timer = @timed optimize!(DW_dual)


    if (termination_status(DW_dual) == MOI.OPTIMAL)
        verb >= 1 && println("Π := $(value.(Π))") # Display solution (if verbosity >= 1)
    else
        @warn "[solve_DW_dual]: dual of DW formulation is not optimal"
    end

    return value.(DW_dual[:Π])
end