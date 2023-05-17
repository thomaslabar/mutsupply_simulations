"""
This file runs the analyses on a prespecified CSV file outputted from sim.jl
"""

using CSV
using DataFrames
using Statistics

function get_mutation_stats(df::DataFrame)

    mutations = df.mutations

    #The array of mutational effects is read as a string by Julia. The below line handles this by 
    #removing the string quotations, spliting up the data by commas, then parses to Float64. This results
    # in an array of arrays. NOTE: This currently can't handle empty arrays of Float64[] text.
    mutation_data = [parse.(Float64, split(chop(mutations[i]; head=1, tail=1), ',')) for i in eachindex(mutations)]

    #flattens the above array of arrays to one array
    mutation_data = reduce(vcat,mutation_data)

    println(mean(mutation_data))
    println(std(mutation_data))
end

function run_analysis(csv_file::String)
    df = DataFrame(CSV.File(csv_file))
    
    large_pop_data = filter(:population_size => n -> n == 100000000, df)
    get_mutation_stats(large_pop_data)

    small_pop_data = filter(:generation => n -> n == 569, df)
    get_mutation_stats(small_pop_data)
    
end

run_analysis("test.csv")