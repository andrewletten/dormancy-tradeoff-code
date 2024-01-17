####################################################################
## dormancy-tradeoff-spp-params.R:  ################################
## Letten, Yamamichi, Richardson and Ke ############################
## Microbial dormancy supports multi-species coexistence ###########
## under resource fluctuations #####################################
## #################################################################

### Initialise parameters for gleaner (2), opportunist (3) and dormancy (1) strategist
parameters <- c(mu1 = 0.05, # dormancy specialist (active phase) max growth rate
                Ks1 = 1, # dormancy specialist (active phase) half saturation constant
                mu2 = 0.07, # gleaner max growth rate
                Ks2 = 0.1, # gleaner half saturation constant
                mu3 = 0.3, # opportunist max growth rate
                Ks3 = 14, # opportunist half saturation constant
                d = 0, # dilution (resources and cells) rate (only relevant for chemostat dynamics) 
                m = 0.025, # mortality rate of non-dormancy specialists (gleaner or opportunist)
                ma = 0.025, # mortality rate of dormancy specialist (active phase)
                md = 0.0038, # mortality rate of dormancy specialist (dormant phase)
                Qs1 = 0.01, # resource quotas
                Qs2 = 0.01, # resource quotas
                Qs3 = 0.01, # resource quotas
                dorm = 1, # change to zero to inhibit switching to dormancy
                switch_dorm = 0.5, # exponential rate of switching to dormant phase
                switch_active = 0.03, # linear rate of returning to active phase
                R0 = 5) # resource supply concentration under chemostat dynamics (ignored when `d` = 0)
