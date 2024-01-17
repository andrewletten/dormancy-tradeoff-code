####################################################################
## dormancy-tradeoff-3sppsims-chemo.R: continuous resource supply ##
## Letten, Yamamichi, Richardson and Ke ############################
## Microbial dormancy supports multi-species coexistence ###########
## under resource fluctuations #####################################
####################################################################

# Core code for simulating three-way and two-way competition under 
# continuous/chemostat resource supply

### Load required packages and project code
library(deSolve)
library(tidyverse)
source("code/dormancy-tradeoff-functions.R")
source("code/dormancy-tradeoff-spp-params.R")

### Modify outflow/mortality parameters for chemostat dynamics ----------------
parameters["d"] <-  0.025 # outflow/mortality rate for all species
parameters["m"] <-  0 # 0 additional mortality for non-dormancy strategists
parameters["ma"] <-  0 # 0 additional mortality for dormancy specialist (active phase)
parameters["md"] <-  -0.0212 # reduced mortality for dormancy specialist (dormant phase)

### Change starting density of either N1_active, N2 or N3 to zero, to simulate pairwise dynamics 
state <- c(N1_activ = 100,
           N1_dorm = 0,
           N2 = 100,
           N3 = 100,
           R = 10)

### Loop sims over different resource supply concentrations -------------------
ressize <-  seq(5, 20, 1) # Resource concentrations to match pulsed dynamics

dat_list = list()
for (i in 1:length(ressize)){
  totaltime <- 2000
  times <- seq(0, totaltime, by = 1)
  start_time <- Sys.time()
  state["R"] = ressize[i]
  parameters["R0"] = ressize[i]
  out <- ode(y = state, 
           times = times, 
           func = dorm_cr_ode, 
           parms = parameters)
  dat_iter <- out %>% data.frame() 
  dat_list[[i]] <-  dat_iter
  print(i)
  end_time <- Sys.time()
  print(end_time - start_time)
}

### Save output (either as compressed csv files OR as a single .Rdata file) ---
#################################################
## Save individual csv.gz for each sim output
# sppcombo <- "N1N2N3" # Change for pairwise combinations, e.g. "N1N2")
# dir.create(file.path(paste0("data/", sppcombo, "_gz/"))) # same directory as for core sims (pulsed resource supply)
# filenames <- paste(sppcombo, paste0("p", ressize), "i0", sep = "_")
# names(dat_list) <- filenames
# iwalk(dat_list, .progress = TRUE,
#       function(x, y) {write_csv(x, file = paste0("data/", sppcombo, "_gz/", y, ".csv.gz"))})
#################################################

#################################################
## Save as single .Rdata file 
# save(dat_list, file = "data/<FILENAME>.Rdata")
#################################################




