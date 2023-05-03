using Distributions
using Random
using LinearAlgebra
using CSV
using DataFrames
using ConfParser

struct Clade
    id::Int
    ancestor::Int
    fitness::Float64
    mutation_rate::Float64
    s_mean::Float64
    mutations::Array{Float64}
    individuals::Int
end

struct Population
    clades::Array{Clade}
    max_id::Int
    populationsize::Int64
    generation::Int
end

struct Parameters
    populationsize::Int
    mutationrate::Float64
    s_ben::Float64
    ancestorfitness::Float64
    replicates::Int
    random_seed::Int
    append::Bool
end

mutable struct Experiment
    const params::Parameters
    largereplicates::Array{Population}
    smallreplicates_equalgen::Array{Population}
    smallreplicates_equalfitness::Array{Population}
    smallreplicates_equalsupply::Array{Population}
    generation::Int
    large_repairgen::Int
    small_repairgen::Int
    small_supplygen::Int
end

function mutation(population::Population)
    new_clades = Clade[]
    new_max_id = population.max_id

    for clade in population.clades
        
        #handle mutants
        num_mutants = min(rand(Poisson(clade.individuals*clade.mutation_rate)),clade.individuals)
        for m in 1:num_mutants
            new_max_id += 1
            mutant_s_effect = rand(Exponential(clade.s_mean))
            new_mutations = deepcopy(clade.mutations)
            push!(new_mutations,mutant_s_effect)
            mutant_clade = Clade(new_max_id,clade.id,clade.fitness+mutant_s_effect,clade.mutation_rate,clade.s_mean,new_mutations,1)
            push!(new_clades,mutant_clade)
        end

        #handle non-mutants
        if clade.individuals - num_mutants > 0
            push!(new_clades,Clade(clade.id,clade.ancestor,clade.fitness,clade.mutation_rate,clade.s_mean,clade.mutations,clade.individuals - num_mutants))
        end
    end

    return Population(new_clades,new_max_id,population.populationsize,population.generation+1)     
end

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

function print_abundant_clades(population::Population, num_clades::Int)
    clades = deepcopy(population.clades)
    sort!(clades, by = v -> v.individuals, rev = true)
    for i in 1:num_clades
        println(clades[i])
    end
end

function get_abundant_fitness(population::Population)
    clades = deepcopy(population.clades)
    sort!(clades, by = v -> v.individuals, rev = true)
    return clades[1].fitness
end

function save_populations(filename::String, pops::Array{Population}, params::Parameters)
    finalresults = DataFrame(random_seed = Int[],
                             replicate = Int64[],
                             population_size = Integer[],
                             mutation_rate = Float64[],
                             DFE_mean = Float64[],
                             generation = Integer[],
                             num_individuals = Integer[],
                             num_mutations = Integer[],
                             fitness = Float64[],
                             mutations = Array{Float64}[]
                             )
    
    add_pops_to_dataframe!(finalresults,pops,params)   
    CSV.write(filename, finalresults, append = params.append)    
end

function add_pops_to_dataframe!(df::DataFrame,populations::Array{Population},params::Parameters)
    for (i,pop) in enumerate(populations) 
        clades = deepcopy(pop.clades)
        sort!(clades, by = v -> v.individuals, rev = true)
        clade = clades[1]
        push!(df,[params.random_seed,i,pop.populationsize,params.mutationrate,params.s_ben,pop.generation,
                            clade.individuals,length(clade.mutations),clade.fitness,clade.mutations])
    end
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

    conf = ConfParse("sim.cfg")
    parse_conf!(conf)

    rs = parse(Int,retrieve(conf, "RANDOM_SEED"))
    if rs == 0
        Random.seed!()
    else
        Random.seed!(rs)
    end

    params = Parameters(parse(Int,retrieve(conf, "POPULATION_SIZE")),
                        parse(Float64,retrieve(conf, "MUTATION_RATE")),
                        parse(Float64,retrieve(conf, "BENEFICIAL_EFFECT")),
                        parse(Float64,retrieve(conf, "ANCESTOR_FITNESS")),
                        parse(Int,retrieve(conf, "NUM_REPLICATES")),
                        parse(Int,retrieve(conf, "RANDOM_SEED")),
                        parse(Bool,retrieve(conf, "APPEND_SAVE_FILE")))
    populations = [Population([Clade(1,0,params.ancestorfitness,params.mutationrate,params.s_ben,[],
                               params.populationsize)],1,
                               params.populationsize,0) for i in 1:params.replicates]
   
    for i in 1:params.replicates
        populations[i] = evolve(populations[i],"Fitness",1.0)
        print_abundant_clades(populations[i],1)
    end
    save_populations("test.csv",populations,params)

end

run_sim()