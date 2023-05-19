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
    This function reads in data from the specified config file (sim.cfg) and returns a 
    parameter struct to store the data
    """

    conf = ConfParse(config_file)
    parse_conf!(conf)

    params = Parameters(parse(Int,retrieve(conf, "LARGE_POPULATION_SIZE")),
                        parse(Int,retrieve(conf, "SMALL_POPULATION_SIZE")),
                        parse(Float64,retrieve(conf, "MUTATION_RATE")),
                        parse(Float64,retrieve(conf, "BENEFICIAL_EFFECT")),
                        parse(Float64,retrieve(conf, "ANCESTOR_FITNESS")),
                        parse(Int,retrieve(conf, "NUM_REPLICATES")),
                        parse(Int,retrieve(conf, "RANDOM_SEED")))

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
    CSV.write(filename, finalresults, append = append)    
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