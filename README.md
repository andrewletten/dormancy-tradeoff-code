## Microbial dormancy supports multi-species coexistence under resource fluctuations

Code and simulation output from:

Letten, Yamamichi, Richardson and Ke. Microbial dormancy supports multi-species coexistence under resource fluctuations.

#### Core scripts

`dormancy-tradeoff-fig1acd.R` - Simulate competition between a gleaner and a dormancy strategist (Fig 1a,c&d in ms).  

`dormancy-tradeoff-fig1b.R` - Simulate competition between a gleaner and a dormancy strategist with different parameter combinations for the dormancy strategist's growth and dormant phase mortality (Fig 1b in ms).  

`dormancy-tradeoff-3sppsims.R` - Simulate competition between all combinations of a gleaner, opportunist and dormancy strategist under different patterns of resource supply.  

`dormancy-tradeoff-post-sim-analysis.R` - Summarize and plot simulation output from `dormancy-tradeoff-3sppsims.R` (Fig 3 in ms).
 
#### Other scripts

`dormancy-tradeoff-spp-params.R` - Parameters for core analyses.

`dormancy-tradeoff-functions.R` - Functions for defining ODE, extracting and wrangling saved sims, plotting etc.

