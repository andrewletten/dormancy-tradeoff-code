####################################################################
## dormancy-tradeoff-fig1acd.R: reproduce fig 1b ###################
## Letten, Yamamichi, Richardson and Ke ############################
## Microbial dormancy supports multi-species coexistence ###########
## under resource fluctuations #####################################
####################################################################

### Load required packages and project code -----------------------------------
library(deSolve)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot()) 
source("code/dormancy-tradeoff-functions.R")
source("code/dormancy-tradeoff-spp-params.R")
#################################

### Initialise state variables ------------------------------------------------
state <- c(N1_activ = 100, # initial density of active phase of dormant strategist
           N1_dorm = 0, # initial density of dormant phase of dormant strategist
           N2 = 100, # initial density of `gleaner`
           N3 = 0, # initial density of `opportunist`
           R = 10) 

### Resource dynamics ---------------------------------------------------------
totaltime <- 480000 # total simulation time
times <- seq(0, totaltime, by = 1)
pulsefreq <- 240 # resource pulse interval
pulsesize <- 8 # resource pulse magnitude
pulseseq <- round(seq(pulsefreq, totaltime, by = pulsefreq), 1)

### Simulate over parameter space ---------------------------------------------
ks_seq <- 0.05*2^(seq(0, 8, 1)) # Ks1 sequence (active phase)
md_seq <- 0.025/(2^seq(0, 8, by = 1)) # md sequence (dormant phase)
dorm_cost <- expand_grid(ks_seq, md_seq)
names(dorm_cost) <- c("Ks1", "md")

dat_list = list()
for (i in 1:(nrow(dorm_cost))){
  parameters["Ks1"] <- dorm_cost[i,1]
  parameters["md"] <- dorm_cost[i,2]
  start_time <- Sys.time()
  out <- ode(y = state, 
             times = times, 
             func = dorm_cr_ode, 
             parms = parameters,
             events = list(func = respulse, time = pulseseq))
  dat_iter <- out %>% data.frame() 
  dat_list[[i]] <-  dat_iter
  names(dat_list)[i] <- paste0("N1N2_Ks", dorm_cost[i,1], 
                            "_md", dorm_cost[i,2]*1000, 
                            "_p", pulsesize, 
                            "_i", pulsefreq)   
  print(i)
  end_time <- Sys.time()
  print(end_time - start_time)
}

### Summarize coexistence/exclusion outcomes ----------------------------------
results <- summarise_outcomes(dat_list, dorm_cost)

### Make heat map of competitive outcomes -------------------------------------
unicols <- unique(results$colscheme)
fig1b <- ggplot(results %>% filter(md > 0.000001 & Ks1 < 25), 
                    aes(x = as.factor(round(md*1000, 2)), y = as.factor(Ks1))) + 
  geom_tile(aes(fill = spp)) + 
  scale_fill_manual(values = unicols) + 
  theme_cowplot() +
  theme(legend.key.size = unit(0.5, 'cm'), 
        legend.text = element_text(size=10),
        legend.position="bottom",
        legend.justification = "center") + 
  guides(fill = guide_legend(nrow = 1, 
                             title="", 
                             label.position = "bottom")) +
  xlab('Mortality rate (dormant phase)') + ylab('Monod half saturation (active phase)') +
  coord_cartesian(expand = FALSE) +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.position = "none") +
  annotate("text", x=1.7, y=1.75, label="Dormant\nwins",
             color="grey20", size = 4) +
  annotate("text", x=4.25, y=4.5, label="Coexistence",
           color="grey20", size = 4) +
  annotate("text", x=7.5, y=8, label="Gleaner\nwins",
           color="grey20", size = 4) +
  geom_point(x = 3, y = 5.25, shape = 18, size = 3)
fig1b

## Check dynamics are stationary ----------------------------------------------
stat_dynamics(dat_list)










