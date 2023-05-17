"""
This file runs the analyses on a prespecified CSV file outputted from sim.jl
"""

using CSV
using DataFrames
using Statistics

function convert_mutation_data_to_array(mutations::Array{String})

    """
    This function takes in an array of strings, each of which represents an array of floats 
        that are the mutational effects for one population. It converts that to one array of floats
        (so, it combines data from every population together)

    It returns that array of floats
    """
    
    #The array of mutational effects is read as a string by Julia. The below line handles this by 
    #removing the string quotations, spliting up the data by commas, then parses to Float64. This results
    # in an array of arrays.
    mutation_strings = [chop(mutations[i]; head=1, tail=1) for i in eachindex(mutations)] #removes '[' and ']' from strings
    mutation_strings = [split(mutation_strings[i], ',') for i in eachindex(mutation_strings)]
    mutation_strings = [i for i in mutation_strings if i!= ["loat64["]] #eliminates empty arrays (for pops with no fixed mutations)

    mutation_data = [parse.(Float64, mutation_strings[i]) for i in eachindex(mutation_strings)]

    #flattens the above array of arrays to one array
    mutation_data = reduce(vcat,mutation_data)

    return mutation_data
end

function get_mutation_stats(df::DataFrame)

    """
    This function is currently a work in progress, but will eventually return a dataframe with each treatment's 
    mean and standard deviation of their mutational effects. It will also return the mean/standard deviation of the
    top N mutations, where N is the number of mutations from the large populations, in order to compare the
    extremes of both.

    """
    @assert "mutations" in names(df)    
    mutation_data = convert_mutation_data_to_array(df.mutations)

    
    return length(mutation_data), mean(mutation_data), std(mutation_data)
end

function run_analysis()
    """
    This function runs all the analyses.
    """

    csv_file  = "test.csv"
    df = DataFrame(CSV.File(csv_file))
    
    large_pop_data = filter(:population_size => n -> n == 100000000, df)
    a,b,c = get_mutation_stats(large_pop_data)

    small_pop_data = filter(:generation => n -> n == 569, df)
    a,b,c = get_mutation_stats(small_pop_data)

    small_totalmutsupply_pop_data = filter(:generation => n -> n == 56900, df)
    a,b,c = get_mutation_stats(small_totalmutsupply_pop_data)
    
end

run_analysis()