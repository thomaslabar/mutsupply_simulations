"""
This file runs the analyses on a prespecified CSV file outputted from sim.jl
"""

using CSV
using DataFrames
using Statistics

function run_analysis(csv_file::String)
    df = DataFrame(CSV.File(csv_file))

    large_pop_data = filter(:population_size => n -> n == 100000000, df)

    large_mutations = large_pop_data.mutations

    #The array of mutational effects is read as a string by Julia. The below line handles this by 
    #removing the string quotations, spliting up the data by commas, then converting to Float64. This results
    # in an array of arrays.
    mutation_data = [parse.(Float64, split(chop(large_mutations[i]; head=1, tail=1), ',')) for i in eachindex(large_mutations)]
    
    #flattens the above array of arrays to one array
    large_mutation_data = reduce(vcat,mutation_data)
   
    println(mean(large_mutation_data))
    println(std(large_mutation_data))

    
end

run_analysis("test.csv")