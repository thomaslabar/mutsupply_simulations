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
        (it combines data from every population together)

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

function save_mutation_stats(df::DataFrame, filename::String)

    """
    This function is currently a work in progress, but will eventually return a dataframe with each treatment's 
    mean and standard deviation of their mutational effects. It will also return the mean/standard deviation of the
    top N mutations, where N is the number of mutations from the large populations, in order to compare the
    extremes of both.

    """
    @assert "mutations" in names(df)    
    @assert "treatment" in names(df)  

    mut_stats = DataFrame(treatment = String[],
                          criteria = String[],
                          number = Int[],
                          mean = Float64[],
                          st_dev = Float64[])

    treatment_names = sort(unique(df.treatment)) #I sort this so "Large_repaired is first, which is required in the loop below)

    large_n = 0

    for t in treatment_names
        treatment_data = filter(:treatment => n -> n == t, df)
        mutation_data = convert_mutation_data_to_array(treatment_data.mutations)
        push!(mut_stats,[t,"All",length(mutation_data),mean(mutation_data),std(mutation_data)])
        if t == "Large_repaired"
            large_n = length(mutation_data)
        else
            #This code looks only at the n mutations with the largest effect, where n is the total number of 
            #mutations across all large populations
            top_mutation_data = sort(mutation_data, rev = true)[1:min(length(mutation_data),large_n)]
            push!(mut_stats,[t,"Top_n",length(top_mutation_data),mean(top_mutation_data),std(top_mutation_data)])
        end
    end

    
    CSV.write(filename, mut_stats)    
end

function run_analysis()
    """
    This function runs all the analyses.
    """

    csv_file  = "test.csv"
    df = DataFrame(CSV.File(csv_file))
    
    save_mutation_stats(df, "test_mutation_stats.csv")
end

run_analysis()