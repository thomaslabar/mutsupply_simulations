using Distributions
using LinearAlgebra
using CSV

include("sim_functions.jl")

function run_sim()

    params = get_parameters("sim.cfg")

    if params.random_seed == 0
        Random.seed!()
    else
        Random.seed!(params.random_seed)
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
    save_populations(params.savefile,large_populations,params,"Large_repaired",false)
    println("Large populations repaired and saved")

    #Evolve small populations for same number of generations as mean large population
    gen_target = ceil(get_avg_generation(large_populations))
    for i in 1:params.replicates
        small_populations[i] = evolve(small_populations[i],"Generation",gen_target)
    end
    save_populations(params.savefile,small_populations,params,"Small_equalgen",true)
    println("Small populations evolved for same avg. generation as large populations and saved")

    #Evolve small populations to restore fitness
    for i in 1:params.replicates
        small_populations[i] = evolve(small_populations[i],"Fitness",1.0)
    end
    save_populations(params.savefile,small_populations,params,"Small_repaired",true)
    println("Small populations repaired and saved")

    #Evolve small populations to same total mutation supply (N*mu*g) as large populations
    for i in 1:params.replicates
        small_populations[i] = evolve(small_populations[i],"Generation",gen_target*(params.largepopulationsize/params.smallpopulationsize))
    end
    save_populations(params.savefile,small_populations,params,"Small_totalsupply",true)
    println("Small populations evolved for same total mutation supply as large populations and saved")
end

run_sim()