using Random, Graphs, JuMP, HiGHS, DelimitedFiles, Distributions, CSV, DataFrames

include("solver.jl")
include("brute_force_approach.jl")
include("KEP_readfile.jl")
include("formulations_compactes.jl")

@warn "To use these tests, you need to download necessary the KEP data available on Preflib (Kidney Data 00036) "



# name of data files
f_prefix = Dict{Int64,Vector{String}}()
f_prefix[16] = vcat(["data_KEP/all_0/00036-0000000$(i).wmd" for i in 1:9], ["data_KEP/all_0/00036-00000010.wmd"])
f_prefix[32] = vcat(["data_KEP/all_0/00036-0000003$(i).wmd" for i in 1:9], ["data_KEP/all_0/00036-00000040.wmd"])
f_prefix[64] = vcat(["data_KEP/all_0/00036-0000007$(i).wmd" for i in 1:9], ["data_KEP/all_0/00036-00000080.wmd"])
f_prefix[128] = vcat(["data_KEP/all_0/00036-0000011$(i).wmd" for i in 1:9], ["data_KEP/all_0/00036-00000120.wmd"])
f_prefix[256] = vcat(["data_KEP/all_0/00036-0000015$(i).wmd" for i in 1:9], ["data_KEP/all_0/00036-00000160.wmd"])
f_prefix[512] = vcat(["data_KEP/all_0/00036-0000019$(i).wmd" for i in 1:9], ["data_KEP/all_0/00036-00000200.wmd"])
f_prefix[1024] = vcat(["data_KEP/all_0/00036-0000023$(i).wmd" for i in 1:9], ["data_KEP/all_0/00036-00000240.wmd"])
f_prefix[2048] = vcat(["data_KEP/all_0/00036-0000027$(i).wmd" for i in 1:9], ["data_KEP/all_0/00036-00000280.wmd"])



"""
compares the various choice for the first cycles
"""
function init_choice_comparaison()
    init_choices = ["half K=2", "all K=2", "path in G_o_prime"]

    # pre-solve to load function
    test_init_ILP = KEP_test("data_KEP/all_0/00036-00000001.wmd"; SP_method="ILP")
    solve_KEP(test_init_ILP)

    instances_sizes = [32, 64, 128, 256]

    times = zeros(length(instances_sizes), length(init_choices))
    losses = zeros(length(instances_sizes), length(init_choices))

    for i in eachindex(instances_sizes)
        n = instances_sizes[i]
        @show n
        for j in eachindex(init_choices)
            for k in 1:10
                kep_test = KEP_test(f_prefix[n][k]; init_choice=init_choices[j], SP_method="ILP")
                timer = @timed begin
                    _, n_relax, n_in = solve_KEP(kep_test)
                end
                times[i, j] += timer.time
                losses[i, j] += n_relax - n_in
            end
            times[i, j] /= 10 # avg
            losses[i, j] /= (10 * n) # normalize with instance size
        end
    end

    # Creates dataframe to export
    df_times = DataFrame(times, :auto)
    rename!(df_times, [Symbol(init_choice) for init_choice in init_choices])
    insertcols!(df_times, 1, :Instance_Size => instances_sizes)

    df_losses = DataFrame(losses, :auto)
    rename!(df_losses, [Symbol(init_choice) for init_choice in init_choices])
    insertcols!(df_losses, 1, :Instance_Size => instances_sizes)

    # Exporter into CSV
    CSV.write("Results/init_comparaison_time.csv", df_times)
    CSV.write("Results/init_comparaison_losses.csv", df_losses)

    return times, losses
end

"""
compares the various choice for the subproblems solving order
"""
function sp_order_comparaison()
    order_choice = ["base", "sequence", "random"]

    # pre-solve to load function
    for order in order_choice
        test_init = KEP_test("data_KEP/all_0/00036-00000001.wmd"; SP_method="Bellmann")
        solve_KEP(test_init)
    end

    instances_sizes = [32, 64, 128]

    times = zeros(length(instances_sizes), length(order_choice))
    losses = zeros(length(instances_sizes), length(order_choice))

    for i in eachindex(instances_sizes)
        n = instances_sizes[i]
        @show n
        for j in eachindex(order_choice)
            for k in 1:10
                kep_test = KEP_test(f_prefix[n][k]; SP_order=order_choice[j], SP_method="Bellmann")
                timer = @timed begin
                    _, n_relax, n_int = solve_KEP(kep_test)
                end
                times[i, j] += timer.time
                losses[i, j] += n_relax - n_int
            end
            times[i, j] /= 10 # avg
            losses[i, j] /= (10 * n) # normalize with instance size
        end
    end

    # Creates dataframe to export
    df_times = DataFrame(times, :auto)
    rename!(df_times, [Symbol(order) for order in order_choice])
    insertcols!(df_times, 1, :Instance_Size => instances_sizes)

    df_losses = DataFrame(losses, :auto)
    rename!(df_losses, [Symbol(order) for order in order_choice])
    insertcols!(df_losses, 1, :Instance_Size => instances_sizes)

    # Exporter into CSV
    CSV.write("Results/order_comparaison_time.csv", df_times)
    CSV.write("Results/order_comparaison_losses.csv", df_losses)

    return times, losses
end


""" 
speed comparaison between the given methods 
"""
function methods_speed_comparaison()

    instances_sizes = [32, 64, 128]
    methods = ["brute force", "compact_formulation", "ILP", "Bellmann"]

    times = zeros(length(instances_sizes), length(methods))

    for i in eachindex(instances_sizes)
        n = instances_sizes[i]
        @show n
        for k in 1:10

            @show k
            # brute_force 
            println("\t Brute force")
            G, W = read_wmd_file(f_prefix[n][k])
            timer = @timed brut_force(G, 4; verb=-1) # we take 4 transfert max by default
            times[i, 1] += timer.time

            # compact formulation
            println("\t Compact formulation")
            timer = @timed solve_cycle_K(G, W, 4; verb=-1)
            times[i, 2] += timer.time


            # ILP column generation
            println("\t ILP")
            kep_test = KEP_test(f_prefix[n][k]; init_choice="all K=2", SP_method="ILP")
            timer = @timed solve_KEP(kep_test)
            times[i, 3] += timer.time


            # Bellmann column generation
            println("\t Bellmann")
            kep_test = KEP_test(f_prefix[n][k]; init_choice="all K=2", SP_method="Bellmann")
            timer = @timed solve_KEP(kep_test)
            times[i, 4] += timer.time

        end
        for j in eachindex(methods)
            times[i, j] /= 10 # avg
        end
    end

    # Creates dataframe to export
    df_times = DataFrame(times, :auto)
    rename!(df_times, [Symbol(m) for m in methods])
    insertcols!(df_times, 1, :Instance_Size => instances_sizes)

    # Exporter into CSV
    CSV.write("Results/method_comparaison_time.csv", df_times)

    return times
end



"""
solve KEP problems even for great instance size to determine the time-behaviour, and look at losses and pair means
"""
function evolution_Bellmann()

    # pre-solve to load function
    test_init = KEP_test("data_KEP/all_0/00036-00000001.wmd"; init_choice="all K=2", SP_method="Bellmann")
    solve_KEP(test_init)


    instances_sizes = [32, 64, 128, 256, 512, 1024, 2048]

    times = zeros(length(instances_sizes))
    losses = zeros(length(instances_sizes))
    means = zeros(length(instances_sizes))
    for i in eachindex(instances_sizes)
        n = instances_sizes[i]
        @show n
        for k in 1:10
            kep_test = KEP_test(f_prefix[n][k]; init_choice="all K=2", SP_method="Bellmann")
            timer = @timed begin
                transferts, n_relax, n_int = solve_KEP(kep_test)
            end
            times[i] += timer.time
            losses[i] += n_relax - n_int
            n_int = sum(length.(transferts))
            means[i] += sum(vcat(transferts...)) / length(vcat(transferts...))  # mean of the selected pairs
        end
        times[i] /= 10 # avg
        losses[i] /= (10 * n) # normalize with instance size
        means[i] /= (n * 10) # avg and normalize by size instance

    end

    # Creates dataframe to export in csv
    df = DataFrame(
        Instance_Sizes=instances_sizes,
        Times=times,
        Losses=losses,
        Means=means
    )
    CSV.write("Results/evolution_Bellmann_2.csv", df)
    return times, losses, means
end




""" 
compute optimum value for a range of K to determine if the is a limit above which increasing K is useless 
"""
function max_K_opti()

    instances_sizes = [32, 256, 1024]
    K_range = 2:7

    times = zeros(length(instances_sizes), length(K_range))
    obj_values = zeros(length(instances_sizes), length(K_range))

    for i in eachindex(instances_sizes)
        n = instances_sizes[i]
        @show n
        for j in eachindex(K_range)
            K = K_range[j]
            println("\tK = $K")
            for k in 1:10

                # Bellmann column generation
                kep_test = KEP_test(f_prefix[n][k]; K=K, init_choice="all K=2", SP_method="Bellmann")
                timer = @timed begin
                    _, _, n_int = solve_KEP(kep_test)
                end
                times[i, j] += timer.time
                obj_values[i, j] += n_int
            end
            times[i, j] /= 10 # avg
            obj_values[i, j] /= (10 * n) # avg + normalization by instance size
        end
    end

    # Creates dataframe to export
    df_times = DataFrame(times, :auto)
    rename!(df_times, [Symbol(K) for K in K_range])
    insertcols!(df_times, 1, :Instance_Size => instances_sizes)

    df_obj = DataFrame(obj_values, :auto)
    rename!(df_obj, [Symbol(K) for K in K_range])
    insertcols!(df_obj, 1, :Instance_Size => instances_sizes)

    # Exporter into CSV
    CSV.write("Results/max_K_opti2.csv", df_obj)
    CSV.write("Results/K_time2.csv", df_times)

    return times
end


" All test scripts and algorithms are successfully loaded "