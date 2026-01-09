DISTRIBUTIONS = ["Constant", "Binomial", "BinomialUNOS", "BinomialAPD", "NoFailure"]

"""
    get_failure_rates

Generate failure rates on each edge, and add its value as a property to the edge of the kep_graph. 

# Parameters
* `kep_graph::SimpleDiGraph` : graph describing the pairs and compatibilities
* `pra::Dict{Int64,Float64}` : PRA de chaque malade
* `distribution::String` : type of distirbution of uncertainties; to be chosen in the DISTRIBUTIONS vector
"""
function get_failure_rates(kep_graph::SimpleDiGraph, pra::Dict{Int64,Float64}, distribution::String)::Dict{Tuple{Int64,Int64},Float64}

    failure_rate = Dict{Tuple{Int64,Int64},Float64}() # probabilité d'échec d'un transfert d'organe du donneur d'une paire vers le patient d'une autre

    for edge in edges(kep_graph)
        # Failure rates depend on the chosen distribution of uncertainties
        if distribution == "Constant"
            # constant failure rates equal to 70%
            failure_rate[edge.src, edge.dst] = 0.7
        elseif distribution == "Binomial"
            if rand() < 0.25
                # random failure rates equal to 10% on average for 25% edges
                failure_rate[edge.src, edge.dst] = rand() * 0.2
            else
                # random failure rates equal to 90% on average for 75% edges
                failure_rate[edge.src, edge.dst] = 0.8 + rand() * 0.2
            end
        elseif distribution == "BinomialUNOS"
            # %pra denotes the panel reactive antibody level
            # %pra of the patient < 0.8 : UNOS low sensitized patients
            if pra[edge.dst] < 0.8
                # failure rate equal to 10% if the patient is low sensitized
                failure_rate[edge.src, edge.dst] = 0.1
            else
                # failure rate equal to 90% otherwise 
                failure_rate[edge.src, edge.dst] = 0.9
            end
        elseif distribution == "BinomialAPD"
            # %pra denotes the panel reactive antibody level
            # %pra of the patient < 0.75 : APD low sensitized patients
            if pra[edge.dst] < 0.75
                # failure rate equal to 28% if the patient is low sensitized
                failure_rate[edge.src, edge.dst] = 0.28
            else
                # failure rate equal to 58% otherwise 
                failure_rate[edge.src, edge.dst] = 0.58
            end
        elseif distribution == "NoFailure"
            # failure rates equal to 0
            failure_rate[edge.src, edge.dst] = 0.
        end
    end

    return failure_rate
end



function check_exchange_failure(G::SimpleDiGraph, exchanges::Vector{Vector{Int}}, pra::Dict{Int64,Float64}, distribution::String="Binomial")
    #Get the failure rates
    failure_rates = get_failure_rates(G, pra, distribution)

    #Number of successful exchanges
    nb_successful_exchanges = 0

    #Successful exchange cycles
    cycles = Vector{Vector{Int}}()

    #Main loop to check if a cycle is successful: a cycle is considered successful if every exchange in the cycle is successful
    for c in exchanges
        #Boolean to know if a cycle is successful
        successful_bool = true

        if rand(Bernoulli(1 - failure_rates[(c[end], c[1])]), 1) == [0]
            successful_bool = false
            continue
        end
        for i in 2:length(c)
            if rand(Bernoulli(1 - failure_rates[(c[i-1], c[i])]), 1) == [0]
                successful_bool = false
                break
            end
        end

        if successful_bool
            push!(cycles, c)
            nb_successful_exchanges += length(c)
        end
    end
    return cycles, nb_successful_exchanges
end

function check_exchange_failure(exchanges::Vector{Vector{Int}}, failure_rates::Dict{Tuple{Int64,Int64},Float64})

    #Number of successful exchanges
    nb_successful_exchanges = 0

    #Successful exchange cycles
    cycles = Vector{Vector{Int}}()

    #Main loop to check if a cycle is successful: a cycle is considered successful if every exchange in the cycle is successful
    for c in exchanges
        #Boolean to know if a cycle is successful
        successful_bool = true

        if rand(Bernoulli(1 - failure_rates[(c[end], c[1])]), 1) == [0]
            successful_bool = false
            continue
        end
        for i in 2:length(c)
            if rand(Bernoulli(1 - failure_rates[(c[i-1], c[i])]), 1) == [0]
                successful_bool = false
                break
            end
        end

        if successful_bool
            push!(cycles, c)
            nb_successful_exchanges += length(c)
        end
    end
    return cycles, nb_successful_exchanges
end