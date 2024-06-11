####################################################################
## dormancy-tradeoff-fig1acd.R: reproduce fig 2 ###################
## Letten, Yamamichi, Richardson and Ke ############################
## Microbial dormancy supports multi-species coexistence ###########
## under resource fluctuations #####################################
## #################################################################

### Load required packages and project code -----------------------------------
library(deSolve)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(patchwork)
theme_set(theme_cowplot()) 
source("code/dormancy-tradeoff-functions.R")
source("code/dormancy-tradeoff-spp-params.R")

## Initialise state variables -------------------------------------------------
state <- c(N1_activ = 100, # initial density of active phase of dormant strategist
           N1_dorm = 0, # initial density of dormant phase of dormant strategist
           N2 = 100, # initial density of `gleaner`
           N3 = 100, # initial density of `opportunist`
           R = 10)


totaltime <- 300000 # total simulation time
times <- seq(0, totaltime, by = 0.1)
pulsefreq <- 120 # resource pulse interval
pulsesize <- 8 # resource pulse magnitude
pulseseq <- round(seq(pulsefreq, totaltime, by = pulsefreq), 1)


### Simulate competition ------------------------------------------------------
out <- ode(y = state, 
           times = times, 
           func = dorm_cr_ode, 
           parms = parameters,
           events = list(func = respulse, time = pulseseq))


### Wrangle data --------------------------------------------------------------
dat <- out %>% data.frame() %>% 
  mutate(Rinflate = R*10) %>% 
  mutate(Ntot = N1_activ + N1_dorm) %>% 
  select(-c(R, N1_activ, N1_dorm)) %>% 
  pivot_longer(-time, names_to = "state_vars", values_to = "abund")
dat$cons_res <- "consumers"
dat$cons_res[dat$state_vars == "Rinflate"] <- "resource"

### Plot Fig 2B ---------------------------------------------------------------
p2_B <- ggplot(dat, aes(x = time, y = abund)) +
  geom_line(aes(col = state_vars), linewidth = 0.8) +
  coord_cartesian(expand = FALSE) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.position = "none") +
  scale_colour_manual(values = c("#0072B2", "#CC79A7", "#E69F00", "#999999")) + 
  xlim(229070, 229070 + pulsefreq*3) +
  ylim(0,320) +
  ylab("Abundance") + xlab("Time")

p2_B


# save resource distribution over a single cycle ------------------------------
respdf <- out %>% data.frame() %>% filter(time > 19200 & time < (19200 + pulsefreq)) 
respdf$Ntot <- respdf$N1_activ + respdf$N1_dorm

# Realised percapita growth response of dormancy strategist -------------------
# (dormant + active phase)
percapNdorm = c()
for (i in 2:length(respdf$Ntot)){
  percapNdorm[i-1] <- (log(respdf$Ntot[i]) - log(respdf$Ntot[i - 1]))/0.1
}

# Per capita growth response of active phase of dormancy strategist 
percap_active <- monod_func(parameters["mu1"], ks = parameters["Ks1"], respdf$R[-1]) - parameters["m"]

# Per capita growth response of gleaner
percap_reg <- monod_func(parameters["mu2"], ks = parameters["Ks2"], respdf$R[-1]) - parameters["m"]

# Per capita growth response of opportunist
percap_opp <- monod_func(parameters["mu3"], ks = parameters["Ks3"], respdf$R[-1]) - parameters["m"]

percapdf <- data.frame(percap_dorm_integrated = percapNdorm, 
                       percap_reg = percap_reg, 
                       percap_opp = percap_opp,
                       resconc = respdf$R[-1], 
                       percap_active = percap_active)
percapgg <- pivot_longer(percapdf, -resconc)

func_resp <- ggplot(percapgg, aes(x = resconc, y = value)) +
  geom_abline(slope = 0, intercept = 0, col = "grey") +
  geom_line(aes(col = name, alpha = name), linewidth = 1) +
  scale_colour_manual(values = c("#E69F00", "#E69F00", "#CC79A7", "#0072B2")) +
  scale_alpha_manual(values = c(1, 0.5, 1, 1)) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.position = "none") + 
  ylab("per capita growth rate") + xlab("Resource") + 
  xlim(0,8)

(func_resp + p2_B) +
  plot_layout(nrow = 1) +
  plot_annotation(tag_levels = 'A')
