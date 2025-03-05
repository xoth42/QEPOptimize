module Visualizer

# Notes for QEPO
"""
generate_dataframe_without_p2s runs the simulation and generates the dataframe, but this file should not be running the simulation at all. It should only govern visualization, and so these functions should only be building and operating on dataframes.  
"""


"""
Module for Quantum Circuit Analysis
This module contains functions for reading IBM quantum computer calibration data,
plotting performance metrics, and analyzing circuit fidelities.
calling include()

1. Reads calibration data from IBM quantum hardware (sing-qubit/ two-qubit gate error rate, readout error, T1, T2).
2. Provides visualizations of circuit performance from a population.
3. Displays top circuits ranked by performance metrics like fidelity.
4. Recalculates performance metrics (fidelity, success probability) for circuits across different input fidelities.
"""

using CSV
using Plots
using DataFrames
using Statistics
using Quantikz
using Random
using QEPO.Optimizer
using BPGates # efficient representation of purification circuits
export plot_performance_metrics, display_top_circuits

# Overriding Quantikz operation for NoisyBellSwap to represent a SWAP operation
Quantikz.QuantikzOp(op::PauliNoiseBellGate{BellSwap}) = Quantikz.SWAP(op.g.idx1, op.g.idx2)

"""
plot all of the performance of population:
Plots the performance metrics of a population of circuits.
It visualizes the fidelity of purified pairs, logical qubit fidelity, average marginal fidelity, and success probability."""
function plot_performance_metrics(population::Population)
    success_probs           = []
    purified_fidelities     = []
    avg_marginal_fidelities = []
    logical_fidelities      = []
    for indiv in population.individuals
        append!(success_probs, indiv.performance.success_probability)
        append!(purified_fidelities, indiv.performance.purified_pairs_fidelity)
        append!(avg_marginal_fidelities, indiv.performance.average_marginal_fidelity)
        append!(logical_fidelities, indiv.performance.logical_qubit_fidelity)
    end

    p1 = plot(
        purified_fidelities,
        label="Purified Pairs Fidelity",
        xlabel="Individual",
        ylabel="Fidelity",
        # lw=0,
        seriestype=:scatter,
        markershape=:square,
        color=:green,
        title="Purified Pairs Fidelity",
        legendfontsize=8,
        titlefontsize=10,
        guidefontsize=8,
        tickfontsize=8,
        ylim=(-0.2, 1.2)
    )

    p2 = plot(
        logical_fidelities,
        label="Logical Qubit Fidelity",
        xlabel="Individual",
        ylabel="Fidelity",
        # lw=1,
        seriestype=:scatter,
        markershape=:circle,
        color=:blue,
        title="Logical Qubit Fidelity",
        legendfontsize=8,
        titlefontsize=10,
        guidefontsize=8,
        tickfontsize=8,
        ylim=(-0.2, 1.2)
    )

    p3 = plot(
        avg_marginal_fidelities,
        label="Average Marginal Fidelity",
        xlabel="Individual",
        ylabel="Fidelity",
        # lw=0,
        seriestype=:scatter,
        markershape=:diamond,
        color=:red,
        title="Average Marginal Fidelity",
        legendfontsize=8,
        titlefontsize=10,
        guidefontsize=8,
        tickfontsize=8,
        ylim=(-0.2, 1.2)
    )

    p4 = plot(
        success_probs,
        label="Success Probability",
        xlabel="Individual",
        ylabel="Probability",
        # lw=0,
        seriestype=:scatter,
        markershape=:cross,
        color=:purple,
        title="Success Probability",
        legendfontsize=8,
        titlefontsize=10,
        guidefontsize=8,
        tickfontsize=8,
        ylim=(-0.2, 1.2)
    )

    plot(p1, p2, p3, p4, layout=(2, 2), size=(1000, 800))
end

"""
    display_top_circuits(individuals::Vector{Individual}, top_n::Int)

Displays the `top_n` circuits from a population, ranked by purified pairs fidelity.
"""
function display_top_circuits(individuals::Vector{Individual}, top_n::Int)

    # Sort the individuals by their purified_pairs_fidelity in descending order
    sorted_individuals = sort(individuals, by = x -> x.performance.purified_pairs_fidelity, rev = true)

    # Display the top `top_n` circuits
    for (i, indiv) in enumerate(sorted_individuals)
        if i > top_n
            break
        end
        displaycircuit(indiv.ops)
        println(indiv.performance)
    end
end


"""
    generate_dataframe_without_p2s(population, num_individuals, f_ins, num_simulations, noise_model=nothing)

Generates a DataFrame that contains simulation statistics for a population of circuits at different fidelities f_in
    and under noise models (T1 and T2).
"""
function generate_dataframe_without_p2s(population, num_individuals, f_ins, num_simulations, noise_model=nothing)
    dataframe_length = num_individuals * length(f_ins)
    r = zeros(Int, dataframe_length)  # Number of registers
    k = zeros(Int, dataframe_length)  # Number of purified pairs
    n = zeros(Int, dataframe_length)  # Number of raw bell pairs
    circuit_length = zeros(Int, dataframe_length)
    purified_pairs_fidelity = zeros(Float64, dataframe_length)
    logical_qubit_fidelity = zeros(Float64, dataframe_length)
    code_distance = zeros(Int, dataframe_length)
    optimized_for = fill(CostFunction(0), dataframe_length)
    average_marginal_fidelity = zeros(Float64, dataframe_length)
    success_probability = zeros(Float64, dataframe_length)
    error_probabilities = fill([0.0], dataframe_length)
    f_in = zeros(Float64, dataframe_length)
    circuit_hash = zeros(UInt64, dataframe_length)
    individual = fill(population.individuals[1], dataframe_length)

    Threads.@threads for i1 in 1:num_individuals
        indiv = population.individuals[i1]
        representative_indiv = refresh_noise(indiv, 0.9)  # for consistency in hashing of the same circuit across different noise rates
        indiv_hash = hash(representative_indiv)

        for i2 in 1:length(f_ins)
            f = f_ins[i2]
            index = (i1-1) * length(f_ins) + i2

            new_indiv = refresh_noise(indiv, f)

            # check if noise_model is provided and adjust the calculation
            if isnothing(noise_model)
                calculate_performance!(new_indiv, num_simulations)
            else
                calculate_performance!(new_indiv, num_simulations, noise_model )
            end

            r[index] = new_indiv.r
            k[index] = new_indiv.k
            n[index] = total_raw_pairs(new_indiv)
            code_distance[index] = new_indiv.code_distance
            optimized_for[index] = new_indiv.optimize_for
            circuit_length[index] = length(new_indiv.ops)
            purified_pairs_fidelity[index] = new_indiv.performance.purified_pairs_fidelity
            logical_qubit_fidelity[index] = new_indiv.performance.logical_qubit_fidelity
            error_probabilities[index] = new_indiv.performance.error_probabilities
            average_marginal_fidelity[index] = new_indiv.performance.average_marginal_fidelity
            success_probability[index] = new_indiv.performance.success_probability
            f_in[index] = f
            circuit_hash[index] = UInt64(indiv_hash)
            individual[index] = representative_indiv
        end
    end

    df = DataFrame(
        r = r,
        k = k,
        n = n,
        circuit_length = circuit_length,
        purified_pairs_fidelity = purified_pairs_fidelity,
        logical_qubit_fidelity = logical_qubit_fidelity,
        code_distance = code_distance,
        optimized_for = optimized_for,
        average_marginal_fidelity = average_marginal_fidelity,
        success_probability = success_probability,
        error_probabilities = error_probabilities,
        f_in = f_in,
        circuit_hash = circuit_hash,
        individual = individual
    )

    return df
end

# TODO rewrite to not run simulation. Only used to create data frames. 
function generate_dataframe_without_p2s(df::DataFrame, num_individuals, f_ins, num_simulations, noise_model)
    dataframe_length = num_individuals * length(f_ins)
    r = zeros(Int, dataframe_length)  # Number of registers
    k = zeros(Int, dataframe_length)  # Number of purified pairs
    n = zeros(Int, dataframe_length)  # Number of raw bell pairs
    circuit_length = zeros(Int, dataframe_length)
    purified_pairs_fidelity = zeros(Float64, dataframe_length)
    logical_qubit_fidelity = zeros(Float64, dataframe_length)
    code_distance = zeros(Int, dataframe_length)
    optimized_for = fill(CostFunction(0), dataframe_length)
    average_marginal_fidelity = zeros(Float64, dataframe_length)
    success_probability = zeros(Float64, dataframe_length)
    error_probabilities = fill([0.0], dataframe_length)
    f_in = zeros(Float64, dataframe_length)
    circuit_hash = zeros(UInt64, dataframe_length)
    individual = fill(df.individual[1], dataframe_length)

    Threads.@threads for i1 in 1:num_individuals
        indiv = df.individual[i1]
        representative_indiv = refresh_noise(indiv, 0.9)  # for consistency in hashing of the same circuit across different noise rates
        indiv_hash = hash(representative_indiv)

        for i2 in 1:length(f_ins)
            f = f_ins[i2]
            index = (i1-1) * length(f_ins) + i2

            new_indiv = refresh_noise(indiv, f)

            calculate_performance!(new_indiv, num_simulations, noise_model )

            r[index] = new_indiv.r
            k[index] = new_indiv.k
            n[index] = total_raw_pairs(new_indiv)
            code_distance[index] = new_indiv.code_distance
            optimized_for[index] = new_indiv.optimize_for
            circuit_length[index] = length(new_indiv.ops)
            purified_pairs_fidelity[index] = new_indiv.performance.purified_pairs_fidelity
            logical_qubit_fidelity[index] = new_indiv.performance.logical_qubit_fidelity
            error_probabilities[index] = new_indiv.performance.error_probabilities
            average_marginal_fidelity[index] = new_indiv.performance.average_marginal_fidelity
            success_probability[index] = new_indiv.performance.success_probability
            f_in[index] = f
            circuit_hash[index] = UInt64(indiv_hash)
            individual[index] = representative_indiv
        end
    end

    df = DataFrame(
        r = r,
        k = k,
        n = n,
        circuit_length = circuit_length,
        purified_pairs_fidelity = purified_pairs_fidelity,
        logical_qubit_fidelity = logical_qubit_fidelity,
        code_distance = code_distance,
        optimized_for = optimized_for,
        average_marginal_fidelity = average_marginal_fidelity,
        success_probability = success_probability,
        error_probabilities = error_probabilities,
        f_in = f_in,
        circuit_hash = circuit_hash,
        individual = individual
    )

    return df
end

function plot_success_rate(df::DataFrame)
    # Filter out rows with NaN values in the 'purified_pairs_fidelity' column
    filtered_df = filter(row -> !isnan(row.logical_qubit_fidelity), df)

    # Group by f_in and calculate the mean and standard deviation for each group
    grouped_df = combine(groupby(filtered_df, :f_in),
                         :success_probability => mean => :mean_success_probability,
                         :success_probability => std => :std_success_probability,
                         :success_probability => length => :count)

    # Calculate the standard error
    grouped_df.std_error = grouped_df.std_success_probability ./ sqrt.(grouped_df.count)

    # Sort the DataFrame by f_in
    sorted_df = sort(grouped_df, :f_in)

    # Set the default font family for all elements in the plot
    default(fontfamily="Times")

    # Plotting F_L as a function of f_in with shaded error bars
    p = plot(sorted_df.f_in, sorted_df.mean_success_probability,
             ribbon=(sorted_df.std_error),
             marker=:circle,
             label="Success Probability",
             lw=2,
             ms=4,
             lc=:cyan,
             fc=:lightcyan,
             xlabel="F_in (raw bell pair in the bottom register)",
            #  ylabel="F_out",
             ylims=(0, 1.1),
             legend=:topleft,
             grid=:on,
             title="'ibmq_kolkata' (swap, CNOT and measurement noise)",
             titlefont=10)

end

function plot_purified_pairs_fidelity(df::DataFrame, title::String)
    # Filter out rows with NaN values in the 'purified_pairs_fidelity' column
    filtered_df = filter(row -> !isnan(row.purified_pairs_fidelity), df)

    # Group by f_in and calculate the mean and standard deviation for each group
    grouped_df = combine(groupby(filtered_df, :f_in),
                         :purified_pairs_fidelity => mean => :mean_purified_pairs_fidelity,
                         :purified_pairs_fidelity => std => :std_purified_pairs_fidelity,
                         :purified_pairs_fidelity => length => :count)

    # Calculate the standard error
    grouped_df.std_error = grouped_df.std_purified_pairs_fidelity ./ sqrt.(grouped_df.count)

    # Define the upper and lower bounds for the shaded region
    upper_bound = grouped_df.mean_purified_pairs_fidelity .+ grouped_df.std_error
    lower_bound = grouped_df.mean_purified_pairs_fidelity .- grouped_df.std_error

    # Sort the DataFrame by f_in
    sorted_df = sort(grouped_df, :f_in)

    # Set the default font family for all elements in the plot
    default(fontfamily="Times")

    # Plotting F_L as a function of f_in with shaded error bars
    p = plot(sorted_df.f_in, sorted_df.mean_purified_pairs_fidelity,
             ribbon=(sorted_df.std_error),
             marker=:circle,
             label="Purified Pairs Fidelity after distillation",
             lw=2,
             ms=4,
            #  lc=:blue,
            #  fc=:lightblue,
             xlabel="F_in",
             ylabel="F_out",
             ylims=(0.65, 1.01),
             legend=:topleft,
             grid=:on,
             title=title,
             titlefont=10)

    # Plotting y=x line
    plot!(p, sorted_df.f_in, sorted_df.f_in,
          label="y = x",
          lw=2,
          lc=:black,
          linestyle=:dash)

end

function plot_logical_qubit_fidelity(df::DataFrame)
    # Filter out rows with NaN values in the 'purified_pairs_fidelity' column
    filtered_df = filter(row -> !isnan(row.logical_qubit_fidelity), df)

    # Group by f_in and calculate the mean and standard deviation for each group
    grouped_df = combine(groupby(filtered_df, :f_in),
                         :logical_qubit_fidelity => mean => :mean_logical_qubit_fidelity,
                         :logical_qubit_fidelity => std => :std_logical_qubit_fidelity,
                         :logical_qubit_fidelity => length => :count)

    # Calculate the standard error
    grouped_df.std_error = grouped_df.std_logical_qubit_fidelity ./ sqrt.(grouped_df.count)

    # Define the upper and lower bounds for the shaded region
    upper_bound = grouped_df.mean_logical_qubit_fidelity .+ grouped_df.std_error
    lower_bound = grouped_df.mean_logical_qubit_fidelity .- grouped_df.std_error

    # Sort the DataFrame by f_in
    sorted_df = sort(grouped_df, :f_in)

    # Set the default font family for all elements in the plot
    default(fontfamily="Times")

    # Plotting F_L as a function of f_in with shaded error bars
    p = plot(sorted_df.f_in, sorted_df.mean_logical_qubit_fidelity,
             ribbon=(sorted_df.std_error),
             marker=:circle,
             label="Logical Qubit Fidelity",
             lw=2,
             ms=4,
             lc=:blue,
             fc=:lightblue,
             xlabel="F_in",
             ylabel="F_out",
             ylims=(0.7, 1.1),
             legend=:topleft,
             grid=:on,
             title="'ibmq_kolkata' (swap, CNOT and measurement noise)",
             titlefont=10)

    # Plotting y=x line
    plot!(p, sorted_df.f_in, sorted_df.f_in,
          label="y = x",
          lw=2,
          lc=:black,
          linestyle=:dash)

end


####################################################################
# Functions for Analyzing Quantum Circuit Results from a DataFrame #
####################################################################

"""
Convert DataFrame string row to Individual
Converts a string representation of an `Individual` back into an `Individual` object.
This is useful when an `Individual` object has been serialized as a string in a DataFrame."""
function convert_string_to_individual(individual_str::String)
    # Parse the individual string to create an Individual instance
    # Assuming the string is formatted correctly for evaluation
    parsed = eval(Meta.parse(individual_str))
    return parsed
end


"""
Recalculates the performance metrics of each individual in the DataFrame and
updates the relevant columns with the new performance data.
It then plots the purified pairs fidelity.
"""
function performance_of_individual_df(df::DataFrame, title::String)

    # Convert the individual column
    df.individual = map(convert_string_to_individual, df.individual)

    # Initialize the performance column
    df.performance = Vector{Performance}(undef, nrow(df))

    # Calculate performance for each individual and update the DataFrame
    for row in eachrow(df)
        indiv = row.individual
        performance = calculate_performance!(indiv, 1000000)
        row.purified_pairs_fidelity = performance.purified_pairs_fidelity
        row.logical_qubit_fidelity = performance.logical_qubit_fidelity
        row.average_marginal_fidelity = performance.average_marginal_fidelity
        row.success_probability = performance.success_probability
        row.code_distance = indiv.code_distance
        row.optimized_for = string(indiv.optimize_for)
        row.performance = performance
    end

    plot_purified_pairs_fidelity(df, title)

end


"""
Displays the top `top_n` circuits from a DataFrame for a given `f_in` value,
ranked by a selected performance metric.
The performance metric can be fidelity,success probability, or other defined metrics.
"""
function display_top_circuits_df(df::DataFrame, f_in::Float64, top_n::Int, performance_metric::Int)

    # Define the mapping of performance metric
    performance_column, title = if performance_metric == 0
        (:logical_qubit_fidelity, "Logical Qubit Fidelity")
    elseif performance_metric == 1
        (:purified_pairs_fidelity, "Purified Pairs Fidelity")
    elseif performance_metric == 2
        (:average_marginal_fidelity, "Average Marginal Fidelity")
    elseif performance_metric == 3
        (:success_probability, "Success Probability")
    else
        error("Invalid performance metric. Use 0 for logical_qubit_fidelity, 1 for purified_pairs_fidelity, 2 for average_marginal_fidelity, 3 for success_rate.")
    end

    df_filtered = filter(row -> row.f_in == f_in , df)

    # Convert the individual column
    df_filtered.individual = map(convert_string_to_individual, df_filtered.individual)

    # # Initialize the performance column
    df_filtered.performance = Vector{Performance}(undef, nrow(df_filtered))

    # Calculate performance for each individual and update the DataFrame
    for row in eachrow(df_filtered)
        indiv = row.individual
        performance = calculate_performance!(indiv, 10000)
        row.purified_pairs_fidelity = performance.purified_pairs_fidelity
        row.logical_qubit_fidelity = performance.logical_qubit_fidelity
        row.average_marginal_fidelity = performance.average_marginal_fidelity
        row.success_probability = performance.success_probability
        row.code_distance = indiv.code_distance
        row.optimized_for = string(indiv.optimize_for)
        row.performance = performance
    end

    # Sort the DataFrame by the selected performance metric in descending order
    sorted_df = sort(df_filtered, performance_column, rev=true)
    # print(sorted_df)

    println("Top $top_n Circuits based on $title:")
    # Display the top `top_n` circuits
    for i in 1:min(top_n, nrow(sorted_df))
        row = sorted_df[i, :]
        indiv = row.individual
        displaycircuit(indiv.ops)
        perf = row.performance
        println("Circuit: ", indiv.ops)
        println("Performance: $perf")
    end

    # plot_purified_pairs_fidelity(sorted_df)

end

"""
Plots the fidelity metrics of a given individual from the DataFrame.
Useful for visualizing how an individual's fidelity changes across different input fidelities.
"""
function plot_fidelities_of_individuals(df::DataFrame, i::Int)

    df.individual = map(convert_string_to_individual, df.individual)

    calculate_performance(df.individual[i])

    # Set the default font family for all elements in the plot
    default(fontfamily="Times")

    # Plotting F_L as a function of f_in
    plot(df.f_in, df.purified_pairs_fidelity,
             marker=:circle,
             label="F_L",
             xlabel="f_in",
             ylabel="F_L",
             grid=true)

end

end