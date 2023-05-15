using Random
using ConfParser

"""
This struct is used to hold all the relevant parameters for the experiment. It is created from
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