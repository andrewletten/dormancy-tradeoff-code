####################################################################
## dormancy-tradeoff-3sppsims.R: core analysis #####################
## Letten, Yamamichi, Richardson and Ke ############################
## Microbial dormancy supports multi-species coexistence ###########
## under resource fluctuations #####################################
####################################################################

# Core code for simulating three-way and two-way competition under 
# different resource pulsing scenarios

### Load required packages and project code
library(deSolve)
library(tidyverse)
source("code/dormancy-tradeoff-functions.R")
source("code/dormancy-tradeoff-spp-params.R")

### Change starting density of either N1_active, N2 or N3 to zero, to simulate pairwise dynamics 
state <- c(N1_activ = 100,
           N1_dorm = 0,
           N2 = 100,
           N3 = 100,
           R = 10)

### Setup resource parameter space to simulate across
respulsesize <-  seq(5, 20, 1)
respulsefreq <-  seq(24, 336, 24)
respulsecombo <-  expand.grid(respulsesize, respulsefreq)
names(respulsecombo) <- c("size", "freq")

### Loop sims over parameter space (patience required!)
dat_list = list()
for (i in 1:nrow(respulsecombo)){
  pulsesize <- respulsecombo[i,1]
  pulsefreq <- respulsecombo[i,2]
  totaltime <- pulsefreq*ifelse((state[["N1_activ"]] > 0 & state[["N2"]] > 0 & state[["N3"]] > 0), 4000, 2000)
  times <- seq(0, totaltime, by = 1)
  pulseseq <- round(seq(pulsefreq, totaltime, by = pulsefreq), 1)
  start_time <- Sys.time()
  out <- ode(y = state, 
           times = times, 
           func = dorm_cr_ode, 
           parms = parameters,
           events = list(func = respulse, time = pulseseq))
  dat_iter <- out %>% data.frame() 
  dat_list[[i]] <-  dat_iter
  print(i)
  end_time <- Sys.time()
  print(end_time - start_time)
}

### Save output (either as compressed csv files OR as a single .Rdata file)
#################################################
## Save individual csv.gz for each sim output
# sppcombo <- "N1N2N3" # Change for pairwise combinations, e.g. "N1N2")
# dir.create(file.path(paste0("data/", sppcombo, "_gz/")))
# filenames <- pmap_vec(respulsecombo,
#                       function(size, freq) {paste(sppcombo, paste0("p", size), paste0("i", freq), sep = "_")})
# names(dat_list) <- filenames
# 
# iwalk(dat_list, .progress = TRUE,
#       function(x, y) {write_csv(x, file = paste0("data/", sppcombo, "_gz/", y, ".csv.gz"))})
#################################################

#################################################
## Save as single .Rdata file 
# save(dat_list, file = "data/<FILENAME>.Rdata")
#################################################
