"""
This file runs the analyses on a prespecified CSV file outputted from sim.jl
"""

using CSV
using DataFrames
using Statistics
using Plots

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
    
    treatment_names = sort(unique(df.treatment)) #I sort this so "Large_repaired is first, which is required in the loop below)
    @assert "Large_repaired" in treatment_names

    mut_stats = DataFrame(treatment = String[],
                          criteria = String[],
                          number = Int[],
                          mean = Float64[],
                          st_dev = Float64[])

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

function plot_mutation_data(df::DataFrame, filename::String)
    """
    This function will, once completed plot the distribution of fixed mutational effects for each treatment. It 
    will also plot comparisons between the large populations and each small-population treatment.

    It makes the following plots:
        1. Histograms of mutational effects for all 4 treatments (test_all.png)
        2. Histograms of mutational effects for top n mutations for all 4 treatments (test_topn.png)
        3. Histogroms comparing "Large_repaired_all" vs. "Small_equalgen_all" and "Large_repaired_all" vs. "Small_repaired_all"
        4. Histogroms comparing "Large_repaired_topn" vs. "Small_repaired_topn" and "Large_repaired" vs. "Small_totalsupply_topn"

    Note: All histograms are normalized by probability (i.e., the sum of the bar heights equals 1)
    """

    @assert "mutations" in names(df)    
    @assert "treatment" in names(df) 
    
    treatment_names = sort(unique(df.treatment)) #I sort this so "Large_repaired" is first, which is required in the loop below)
    @assert "Large_repaired" in treatment_names

    large_n = 0
    mutation_df = Dict()
    max_s = -1
    min_s = 100
    for t in treatment_names
        treatment_data = filter(:treatment => n -> n == t, df)
        mutation_data = convert_mutation_data_to_array(treatment_data.mutations)
        max_s = max(max_s,maximum(mutation_data))
        min_s = min(min_s,minimum(mutation_data))
        if t == "Large_repaired"
            large_n = length(mutation_data)
        end
        mutation_df[t*"_all"] = mutation_data
        top_mutation_data = sort(mutation_data, rev = true)[1:min(length(mutation_data),large_n)]
        mutation_df[t*"_top_n"] = top_mutation_data
    end

    #Plot all 4 treatments showing all mutations
    p_all = [histogram(mutation_df[t*"_all"], title = t, normalize=:probability, tickfontsize=6) for t in treatment_names]
    plt_all = plot(p_all[1],p_all[2],p_all[3],p_all[4],layout=(2,2), legend = false)
    xlims!(min_s,max_s)
    xlabel!("Fitness Effect")
    ylabel!("Num. Mutations")
    savefig(plt_all,"test_all.png")
    
    #Plot all 4 treatments showing top n mutations, where n is number in large populations
    p_top = [histogram(mutation_df[t*"_top_n"], title = t, normalize=:probability, tickfontsize=6) for t in treatment_names]
    plt_top = plot(p_top[1],p_top[2],p_top[3],p_top[4],layout=(2,2), legend = false)
    xlims!(min_s,max_s)
    xlabel!("Fitness Effect")
    ylabel!("Num. Mutations")
    savefig(plt_top,"test_topn.png")

    #Need different way to plot overlapping histograms
    #Plot comparisons 1) "Large_repaired_all" vs. "Small_equalgen_all" and 2) "Large_repaired_all" vs. "Small_repaired_all"
    y1 = mutation_df["Large_repaired_all"]
    y2 = mutation_df["Small_equalgen_all"]
    y3 = mutation_df["Small_repaired_all"]
    p1 = histogram([y1,y2], title = "Large_repaired_all vs. \nSmall_equalgen_all", titlefontsize = 10, alpha = 0.5, 
                   normalize=:probability, tickfontsize=6)
    p2 = histogram([y1,y3], title = "Large_repaired_all vs. \nSmall_repaired_all", titlefontsize = 10, alpha = 0.5, 
                   normalize=:probability, tickfontsize=6)
    plt_compare_all = plot(p1,p2, layout=(1,2), legend = false)
    xlims!(min_s,max_s)
    xlabel!("Fitness Effect")
    ylabel!("Num. Mutations")
    savefig(plt_compare_all,"test_compare_all.png")

    #Plot comparisons 1) "Large_repaired_topn" vs. "Small_repaired_topn" and 2) "Large_repaired" vs. "Small_totalsupply_topn"
    y4 = mutation_df["Large_repaired_top_n"]
    y5 = mutation_df["Small_repaired_top_n"]
    y6 = mutation_df["Small_totalsupply_top_n"]
    p1 = histogram([y4,y5], title = "Large_repaired_top_n vs. \nSmall_repaired_top_n", titlefontsize = 10, alpha = 0.5, 
                   normalize=:probability, tickfontsize=6)
    p2 = histogram([y4,y6], title = "Large_repaired_top_n vs. \nSmall_totalsupply_top_n", titlefontsize = 10, alpha = 0.5, 
                   normalize=:probability, tickfontsize=6)
    plt_compare_topn = plot(p1,p2, layout=(1,2), legend = false)
    xlims!(min_s,max_s)
    xlabel!("Fitness Effect")
    ylabel!("Num. Mutations")
    savefig(plt_compare_topn,"test_compare_topn.png")
    
end

function run_analysis()
    """
    This function runs all the analyses.
    """

    csv_file  = "test.csv"
    df = DataFrame(CSV.File(csv_file))
    
    save_mutation_stats(df, "test_mutation_stats.csv")
    plot_mutation_data(df, "test_plots.jpg")
end

run_analysis()