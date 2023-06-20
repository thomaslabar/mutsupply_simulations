using Random
using ConfParser
using DataFrames


"""
The Parameters struct is used to hold all the relevant parameters for the experiment. It is created from
the data in sim.cfg.
"""
struct Parameters
    largepopulationsize::Int
    smallpopulationsize::Int
    mutationrate::Float64
    s_ben::Float64
    ancestorfitness::Float64
    replicates::Int
    random_seed::Int
    savefile::String
end

"""
The Clade struct represents a given genotype (not an individual) and how many individuals are of that genotype.
In addition to the genotype's fitness and the collection of mutations its accumulated over its evolution,
it also contains its beneficial mutation rate and the mean selection coefficient of beneficial mutations 
(the mean of an exponential distribution). These are stored as genotype parameters to allow for their
potential evolution (i.e., global epistasis or mutator strains). However, they are currently kept
constant.

Note that this is not a mutable struct. During every generation, a new Clade object is created for each
genotype
"""
struct Clade
    id::Int
    ancestor::Int
    fitness::Float64
    mutation_rate::Float64
    s_mean::Float64
    mutations::Array{Float64}
    individuals::Int
end

"""
The Population struct is at its heart an array of Clades (see above). It also contains an max_id and 
    population size variables for convienence for some functions (so one does not have to count every
    time). It also keeps track of the (discrete) generations it has evolved for. Like with Clade
    structs, this is immutable and recreated every generation.  
"""
struct Population
    clades::Array{Clade}
    max_id::Int
    populationsize::Int64
    generation::Int
end

function get_parameters(config_file::String)
    """
    This function reads in data from the specified config file (sim.cfg) and the command line
    and then returns a parameter struct to store the data
    """

    conf = ConfParse(config_file)
    parse_conf!(conf)

    #Note: the below code is structured in such a way because a Parameters struct is immutable

    #Read default data in from sim.cfg
    temp_params = Dict()
    temp_params["LARGE_POPULATION_SIZE"] = retrieve(conf, "LARGE_POPULATION_SIZE")
    temp_params["SMALL_POPULATION_SIZE"] = retrieve(conf, "SMALL_POPULATION_SIZE")
    temp_params["MUTATION_RATE"] = retrieve(conf, "MUTATION_RATE")
    temp_params["BENEFICIAL_EFFECT"] = retrieve(conf, "BENEFICIAL_EFFECT")
    temp_params["ANCESTOR_FITNESS"] = retrieve(conf, "ANCESTOR_FITNESS")
    temp_params["NUM_REPLICATES"] = retrieve(conf, "NUM_REPLICATES")
    temp_params["RANDOM_SEED"] = retrieve(conf, "RANDOM_SEED")
    temp_params["SAVE_FILE"] = retrieve(conf, "SAVE_FILE")
    
    #Change default values based on command line inputs
    @assert length(ARGS)%3 == 0 #Should be true because each variable changes requires '-set', '{Parameter}', '{value}'
    true_arg_length = Int(length(ARGS)/3)
    for i in 1:true_arg_length
        @assert ARGS[3*(i-1)+1] == "-set"
        temp_params[ARGS[3*(i-1)+2]] = ARGS[3*i]
    end
    
    #Create Parameters struct
    params = Parameters(parse(Int,temp_params["LARGE_POPULATION_SIZE"]),
                        parse(Int,temp_params["SMALL_POPULATION_SIZE"]),
                        parse(Float64,temp_params["MUTATION_RATE"]),
                        parse(Float64,temp_params["BENEFICIAL_EFFECT"]),
                        parse(Float64,temp_params["ANCESTOR_FITNESS"]),
                        parse(Int,temp_params["NUM_REPLICATES"]),
                        parse(Int,temp_params["RANDOM_SEED"]),
                        temp_params["SAVE_FILE"])

    return params
end

function add_pops_to_dataframe!(df::DataFrame,populations::Array{Population},params::Parameters,treatment::String)
    """
    This function takes as input an array of populations and the simulation parameters and addes them
        to an inputted dataframe. It is used by the function 'save_populations' which writes/appends
        this dataframe to a CSV file. 
        
        Note: it saves the data from the most abundant clade in the population, not the population as a whole.
    """

    for (i,pop) in enumerate(populations) 
        clades = deepcopy(pop.clades)
        sort!(clades, by = v -> v.individuals, rev = true)
        clade = clades[1]
        push!(df,[treatment,i,params.random_seed,pop.populationsize,params.mutationrate,params.s_ben,pop.generation,
                            clade.individuals,length(clade.mutations),clade.fitness,clade.mutations])
    end
end

function save_populations(filename::String, pops::Array{Population}, params::Parameters, treatment::String, append::Bool)
    """
    This function saves data about a population to a dataframe and then to a prespecified CSV file. 
    It can either create a new CSV file or append if the file exists. It saves data about the parameters 
    of the population/experiment and the treatment.
    """

    try
        mkdir("Data")
    catch e
        println("Data directory already exists")
    end

    finalresults = DataFrame(treatment = String[],
                             replicate = Int64[],
                             random_seed = Int[],
                             population_size = Integer[],
                             mutation_rate = Float64[],
                             DFE_mean = Float64[],
                             generation = Integer[],
                             num_individuals = Integer[],
                             num_mutations = Integer[],
                             fitness = Float64[],
                             mutations = Array{Float64}[]
                             )
    
    add_pops_to_dataframe!(finalresults,pops,params,treatment)   
    CSV.write("Data/"*filename, finalresults, append = append)    
end

function get_abundant_fitness(population::Population)
    """
    This function returns the fitness of the most abundant clade (greatest number of individuals)
    in a population.
    """

    clades = deepcopy(population.clades) #I create an extra copy in order to sort without altering the input pop.
    sort!(clades, by = v -> v.individuals, rev = true)
    return clades[1].fitness
end

function get_avg_generation(populations::Array{Population})
    """
    This function calculates the mean generation across an array of replicate populations. it is used
    to calculate the mean number of generations required for the large generations to repair their fitness loss.
    This results is used to determine how long to evolve the small populations in certain treatments.
    """
    return sum([populations[i].generation for i in 1:length(populations)])/length(populations)
end

function print_abundant_clades(population::Population, num_clades::Int)

    """
    This function prints the top "num_clades" most-abundant clades (ranked by most individuals)
    in a given population
    """

    clades = deepcopy(population.clades)
    sort!(clades, by = v -> v.individuals, rev = true)
    for i in 1:num_clades
        println(clades[i])
    end
end

function mutation(population::Population)
    """
    This function takes a population and generates mutations for the next generation. It does not modify
    the current population, but returns a new one. It loops through each clade and calulates the number of 
    mutants by with a Poisson distribution with a mean of the number of individuals in the clade multiplied
    by the clade's mutation rate. Then, it loops through the mutants and samples its fitness effect by 
    generation a random Exponentially-distributed number with mean equal to the clade's s_mean. It creates
    a new clade for each mutant with one individual. Finally, if there are any non-mutant individuals left
    in a clade, it recreates a new clade with the adjusted number (this is necessary because Clades are 
    immutable)
    """
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
            mutant_clade = Clade(new_max_id,clade.id,clade.fitness*(1.0+mutant_s_effect),clade.mutation_rate,clade.s_mean,new_mutations,1)
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
    """
    This function takes a population and performs selection on it (i.e., generates frequencies of each
    clade for the next generation). Like the mutation function, it doesn't modify the input population,
    but returns a new population struct. This function calculates the number of individuals per clade
    in the next generation by first calculating the weighted fitness per clade (where the weights are 
    proportional to the product of the number of individuals in the clade and the clade's relative_fitness
    fitness). These weights are then used for Multinomial sampling to generate the next generation's 
    clade abundances (with total population size kept constant between generations).  
    """

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
    new_clades = [Clade(clade.id,clade.ancestor,clade.fitness,clade.mutation_rate,clade.s_mean,clade.mutations,
                  num_offspring[i]) for (i,clade) in enumerate(population.clades) if num_offspring[i] > 0]
    
    return Population(new_clades,population.max_id,n,population.generation)

end

function evolve(pop::Population, criteria_type::String, criteria::Number)

    """
    This function evolves a population until some criteria is reached. The criteria types can either be
    a number of generations or a target fitness. Evolution works by iterating the population through rounds
    of selection then mutation.
    """

    @assert criteria_type == "Generation" || criteria_type == "Fitness"
    @assert criteria > 0

    if criteria_type == "Generation"
        while pop.generation < criteria
            pop = selection(pop)
            pop = mutation(pop)
        end
    elseif criteria_type == "Fitness"
        while get_abundant_fitness(pop) < criteria
            pop = selection(pop)
            pop = mutation(pop)
        end
    end
    return pop
end