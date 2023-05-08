# Mutation Supply Simulations

This project contains Julia scripts to explore hypotheses concerning the role of population size in altering
the genetic architecture of evolutonary repair in asexual populations.

## Description

This project contains code to simulate the evolution of asexual populations adapting to a 
novel environment with an exponential distribution of beneficial fitness effects (no neutral or deleterious mutations). 
It was designed to simulate an evolutionary repair experiment (cells adapting in response to a genetic perturbation) 
with two population size treatments. A genetic perturbation is simulated by givning the ancestor 
genotype of the experiment a fitness below 1.0 (default ancestral fitness is 0.85). 

The "experiment" works as follows. First, n large populations (default n = 10)
are evolved until the genetic perturbation is repaired (i.e., fitness is restored to 1.0). Then, small
populations evolve until they reach the same number of mean generations required for the large populations
to restore fitness. Next, small populations are evolved until they restore fitness to 1.0. Finally, the same
small populations evolve until they experience the same total experimental mutation supply 
(population size * mutation rate * generations) as the large populations. In other words, if large populations
have 100x more individuals as the small populations, the small populations will evolve for 100x more generations.
At the end of each of these steps, data is saved in CSV format.

The project currently consists of two files:

1. sim.jl, which contains all the structs and functions required to run the simulations in addition to 
actually running the simulation itself. Later, this file will seperate the actual running of the simulation 
from struct and function definitions.

2. sim.cfg, which contains all the parameters that can be changed for the simulation.
    - Large and small population sizes
    - Beneficial mutation rate
    - Mean effect of beneficial mutations (always exponentially-distributed)
    - Ancestor fitness 
    - Number of experimental replicates
    - Random seed
    - Save File Name (all data saved to one file)

## Acknowledgments

* [README template](https://gist.github.com/DomPizzie/7a5ff55ffa9081f2de27c315f5018afc)