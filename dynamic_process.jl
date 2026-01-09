include("bf_solver_ME.jl")
include("check_failure.jl")
include("expectation.jl")

"""
# Returns
* `Int64` : total number of completed transfert
"""
function stochastic_process(G::SimpleDiGraph, pra::Dict{Int64,Float64}; K::Int64=4, sto_method="Determinist", distribution::String="Binomial", nb_months::Int=12, waiting_list_leaving_percentage::Int=1, waiting_list_added_percentage::Int=2, verb::Int64=0)::Tuple{Int64,Int64}
    #Vector of randomly taken vertices in the waiting list for the KEP
    waiting_list = sort(shuffle(vertices(G))[1:convert(Int, floor(nv(G) / 2))])

    #Vector of the remaining vertices in the graph G (the vertices which aren't in the waiting list and didn't participate in a successful exchange)
    v_valid = setdiff(vertices(G), waiting_list)
    if verb >= 1
        println("Initial waiting list : $waiting_list")
        println("============================================================")
    end

    faillure_rates = get_failure_rates(G, pra, distribution)
    total_nb_choosen = 0
    total_nb_success = 0

    #Main loop for the stochastic process with one update every month
    for i in 1:nb_months
        if verb >= 1
            println("Month $i")
            println("-----------------------------")
        end
        #Stop the process if there are less than 2 vertices left in the waiting list
        if length(waiting_list) < 2
            break
        end

        #Subgraph of G with only the vertices of the waiting list
        subgraph, vmap = induced_subgraph(G, waiting_list)

        #Possible exchanges

        if sto_method == "Determinist"
            test = KEP_test(subgraph; K=K, SP_method="Bellmann")
            exchanges, _, nb_choosen = solve_KEP(test)
        elseif sto_method == "Max-Expectation"
            exchanges, nb_choosen = stoc_brut_force(subgraph, K, faillure_rates; v_map=vmap)
        else
            println("[stochastic_process] : unknow stochastic approach $sto_method")
            @error("[stochastic_process] : unknow stochastic approach $sto_method")
        end

        #Which exchanges are truly successful
        exchanges_tmp = map(v -> vmap[v], exchanges)
        successful_exchanges, nb_exchanges = check_exchange_failure(exchanges_tmp, faillure_rates)

        total_nb_choosen += nb_choosen
        total_nb_success += nb_exchanges

        #The vertices that are leaving the waiting list: successful exchange or other reason
        v_filter = sort(shuffle(waiting_list)[1:convert(Int, floor(waiting_list_leaving_percentage / 100 * length(waiting_list)))]) #a given % of vertices leave the waiting list for some reason
        for c in successful_exchanges
            for v in c
                push!(v_filter, v)
            end
        end
        v_filter = sort(v_filter)
        for v in v_filter
            waiting_list = filter!(e -> e ≠ v, waiting_list)
            v_valid = filter!(e -> e ≠ v, v_valid)
        end

        #The vertices that are added to the waiting list
        v_add = sort(shuffle(v_valid)[1:convert(Int, floor(waiting_list_added_percentage / 100 * length(v_valid)))]) #a given % of vertices is added to the waiting list
        for v in v_add
            push!(waiting_list, v)
            v_valid = filter!(e -> e ≠ v, v_valid)
        end
        waiting_list = sort(waiting_list)

        #Print
        if verb >= 1
            println("New pairs added to the waiting list: $v_add")
            println("-----------------------------")
            println("Successful exchanges: $successful_exchanges")
            println("-----------------------------")
            println("Pairs leaving the waiting list: $v_filter")
            println("-----------------------------")
            println("Waiting list: $waiting_list")
            println("============================================================")
        end
    end

    return total_nb_success, total_nb_choosen
end