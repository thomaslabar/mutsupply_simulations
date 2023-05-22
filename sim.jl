using Distributions
using LinearAlgebra
using CSV

include("sim_functions.jl")

function selection(population::Population)

    clade_fitness = Float64[]
    clade_individuals = Float64[]
    for clade in population.clades
        push!(clade_individuals,clade.individuals)
        push!(clade_fitness,clade.fitness)
    end
    n = population.populationsize

    clade_frequencies = [v/n for v in clade_individuals]

    average_fitness = dot(clade_fitness,clade_frequencies)
    relative_fitness = [v/average_fitness for v in clade_fitness]
    weighted_fitness = [v*(clade_individuals[i]/n) for (i, v) in enumerate(relative_fitness)]

    #For the below assertion, I need something like the sum is close enough to 1 so that it works as an expectation function
    #@assert sum(weighted_fitness) == 1.0
    num_offspring = rand(Multinomial(n,weighted_fitness))
    
    new_clades = Clade[]
    new_max_id = -1
    for (i,clade) in enumerate(population.clades)
        if num_offspring[i] > 0
            new_max_id = max(new_max_id,clade.id)
            push!(new_clades,Clade(clade.id,clade.ancestor,clade.fitness,clade.mutation_rate,clade.s_mean,clade.mutations,num_offspring[i]))
        end
    end

    return Population(new_clades,new_max_id,n,population.generation)

end

function evolve(pop::Population, target_type::String, target::Number)

    @assert target_type == "Generation" || target_type == "Fitness"

    if target_type == "Generation"
        while pop.generation < target
            pop = selection(pop)
            pop = mutation(pop)
        end
    elseif target_type == "Fitness"
        while get_abundant_fitness(pop) < target
            pop = selection(pop)
            pop = mutation(pop)
        end
    end
    return pop
end

function run_sim()

    params = get_parameters("sim.cfg")

    if params.random_seed == 0
        Random.seed!()
    else
        Random.seed!(rs)
    end

    large_populations = [Population([Clade(1,0,params.ancestorfitness,params.mutationrate,params.s_ben,[],
                               params.largepopulationsize)],1,
                               params.largepopulationsize,0) for i in 1:params.replicates]
    small_populations = [Population([Clade(1,0,params.ancestorfitness,params.mutationrate,params.s_ben,[],
                               params.smallpopulationsize)],1,
                               params.smallpopulationsize,0) for i in 1:params.replicates]

    #Evolve large populations to restore fitness
    for i in 1:params.replicates
        large_populations[i] = evolve(large_populations[i],"Fitness",1.0)
    end
    save_populations("test.csv",large_populations,params,"Large_repaired",false)
    println("Large populations repaired and saved")

    #Evolve small populations for same number of generations as mean large population
    gen_target = ceil(get_avg_generation(large_populations))
    for i in 1:params.replicates
        small_populations[i] = evolve(small_populations[i],"Generation",gen_target)
    end
    save_populations("test.csv",small_populations,params,"Small_equalgen",true)
    println("Small populations evolved for same avg. generation as large populations and saved")

    #Evolve small populations to restore fitness
    for i in 1:params.replicates
        small_populations[i] = evolve(small_populations[i],"Fitness",1.0)
    end
    save_populations("test.csv",small_populations,params,"Small_repaired",true)
    println("Small populations repaired and saved")

    #Evolve small populations to same total mutation supply (N*mu*g) as large populations
    for i in 1:params.replicates
        small_populations[i] = evolve(small_populations[i],"Generation",gen_target*(params.largepopulationsize/params.smallpopulationsize))
    end
    save_populations("test.csv",small_populations,params,"Small_totalsupply",true)
    println("Small populations evolved for same total mutation supply as large populations and saved")
end

run_sim()