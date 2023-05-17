using Random
using ConfParser


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

function add_pops_to_dataframe!(df::DataFrame,populations::Array{Population},params::Parameters)
    """
    This function takes as input an array of populations and the simulation parameters and addes them
        to an inputted dataframe. It is used by the function 'save_populations' which writes/appends
        this dataframe to a CSV file.
    """

    for (i,pop) in enumerate(populations) 
        clades = deepcopy(pop.clades)
        sort!(clades, by = v -> v.individuals, rev = true)
        clade = clades[1]
        push!(df,[params.random_seed,i,pop.populationsize,params.mutationrate,params.s_ben,pop.generation,
                            clade.individuals,length(clade.mutations),clade.fitness,clade.mutations])
    end
end