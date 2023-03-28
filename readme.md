# Code for "Acquired phototrophy as an evolutionary path to mixotrophy" (2023)
## Alexandra L. Brown*, Grace A. Casarez*, and Holly V. Moeller
*co-first authors
Contact: alexandra_brown (berkeley.edu), hvmoeller (ucsb.edu)

This repository contains the Mathematica code for the simulations and analyses in the manuscript "Acquired phototrophy as an evolutionary path to mixotrophy." 
It models the evolution of acquired phototrophy in a population of phagotrophic plankton (consumers). We study the conditions that cause consumers to evolve to retain the chloroplasts of photosynthetic prey as well as to replicate these chloroplasts so they can be passed on to offspring. In particular, we focus on how consumer evolution is affected by (1) surface light and (2) a trade-off between predation ability and plastid retention replication ability.

To produce the figures in the manuscript, the code was run in Mathematica version 13.

Code was written by Brown and Casarez.

## Files to produce plots
 * `model_and_functions.wls` contains the model equations and functions to analyze ecological and evolutionary dynamics
 * `bifurcation_diagram.wls` contains the code to produce Fig 2
 * `evolution_from_heterotrophy.wls` contains the code to produce Figs 3, 5, S6, and S7
 * `selection_gradient_and_ecological_equilibria.wls` contains the code to produce Figs 4, S3, S4, and S5
  * `alternate_attack_tradeoff.wls` contains the code to produce Figs S1 and S2
 * `phototroph_eats.wls` contains the code to produce Fig S8
 * the folder `phototroph_eating_vs_extinction` contains the code to produce Fig S9
    -  `light_vs_survival_phototroph_doesnt_eat.wls` simulates consumer survival at high light when the phototroph doesn't eat
    - `light_vs_survival_phototroph_eats.wls` simulates consumer survival at high light when the phototroph can eat
    - `plot_high_light_extinction.wls` plots the results of the above simulations
 * `selection_gradient_vs_mortality_parameters.wls` contains the code to produce Fig S10
 * `evolution_vs_photosynthesis.wls` contains the code to produce Fig S11
 * `decay_and_photosynthesis_tradeoffs.wls` contains the code to produce Fig S12

