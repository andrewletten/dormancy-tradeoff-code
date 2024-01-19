####################################################################
## dormancy-tradeoff-fig1acd.R: reproduce fig 1a, c & d ############
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
           N3 = 0, # initial density of `opportunist` (ignored for initial analysis [Fig 1])
           R = 10)


totaltime <- 20000 # total simulation time
times <- seq(0, totaltime, by = 0.1)
pulsefreq <- 240 # resource pulse interval
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
  select(-c(R, N3)) %>% 
  pivot_longer(-time, names_to = "state_vars", values_to = "abund")
dat$cons_res <- "consumers"
dat$cons_res[dat$state_vars == "R"] <- "resource"

### Plot Fig 1A ---------------------------------------------------------------
p1_A <- ggplot(dat, aes(x = time, y = abund)) +
  geom_line(aes(col = state_vars, lty = state_vars), linewidth = 0.8) +
  coord_cartesian(expand = FALSE) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.position = "none") +
  scale_colour_manual(values = c("#E69F00", "#E69F00", "#0072B2", "#999999")) + 
  scale_linetype_manual(values = c("solid", "dashed", "solid", "solid")) + 
  xlim(9110, 9760) + ylim(0,350) + 
  ylab("Abundance") + xlab("Time")

p1_A

# save resource distribution over a single cycle ------------------------------
respdf <- out %>% data.frame() %>% filter(time > 19200 & time < (19200 + pulsefreq)) 
respdf$Ntot <- respdf$N1_activ + respdf$N1_dorm

# Realised percapita growth response of dormancy strategist -------------------
# (dormant + active phase)
percapNdorm = c()
for (i in 2:length(respdf$Ntot)){
  percapNdorm[i-1] <- (log(respdf$Ntot[i]) - log(respdf$Ntot[i - 1]))/0.1
}

# Realised per capita growth response of gleaner 
# Note - equivalent to parameterised Monod function
percap_reg = monod_func(parameters["mu2"], ks = parameters["Ks2"], respdf$R[-1]) - parameters["m"]
percap_reg = c()
for (i in 2:length(respdf$N2)){
  percap_reg[i-1] <- (log(respdf$N2[i]) - log(respdf$N2[i - 1]))/0.1
}

# Per capita growth response of active phase of dormancy strategist 
percap_active = monod_func(parameters["mu1"], ks = parameters["Ks1"], respdf$R[-1]) - parameters["m"]


percapdf <- data.frame(percap_dorm_integrated = percapNdorm, 
                       percap_reg = percap_reg, 
                       resconc = respdf$R[-1], 
                       percap_active = percap_active)
percapgg <- pivot_longer(percapdf, -resconc)

func_resp <- ggplot(percapgg, aes(x = resconc, y = value)) +
  geom_abline(slope = 0, intercept = 0, col = "grey") +
  geom_line(aes(col = name, alpha = name), linewidth = 1) +
  scale_colour_manual(values = c("#E69F00", "#E69F00", "#0072B2")) +
  scale_alpha_manual(values = c(1, 0.5, 1)) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.position = "none") + 
  ylab("per capita growth rate") + xlab("Resource") + 
  xlim(0,8)

func_resp_lowres <- ggplot(percapgg, aes(x = resconc, y = value)) +
  geom_abline(slope = 0, intercept = 0, col = "grey") +
  geom_line(aes(col = name, alpha = name), linewidth = 0.7) +
  scale_colour_manual(values = c("#E69F00", "#E69F00", "#0072B2")) +
  scale_alpha_manual(values = c(1, 0.5, 1)) +
  theme(axis.title = element_blank(),
    axis.text = element_text(size = 6),
    legend.position = "none",
    axis.line=element_line(linewidth=0.1),
    axis.ticks=element_line(linewidth=0.1)) + 
  ylab("per capita growth rate") + xlab("Resource") + 
  xlim(0,0.4) + scale_y_continuous(limits = c(-0.025, 0.032), 
                                   breaks = c(-0.025, 0, 0.025))

### Plot Fig 1C
p1_C <- func_resp + inset_element(func_resp_lowres, 
                          0.5, -0.02, 0.95, 0.36,
                          ignore_tag = TRUE)
p1_C

################################
### Inspect resource variability

# resimulate ode model with only gleaner present
state["N1_activ"] <- 0 # initial density of dormancy strategist
state["N2"] <- 100 # initial density of `gleaner`

# simulate competition
out <- ode(y = state, 
           times = times, 
           func = dorm_cr_ode, 
           parms = parameters,
           events = list(func = respulse, time = pulseseq))

# take a subset of dynamics for a single resource pulse cycle
glean_df <- out %>% data.frame() %>% filter(time > 19200 & time < (19200 + pulsefreq)) 
glean_df$Ntot <- glean_df$N1_activ + glean_df$N1_dorm
var(glean_df$R)

# resimulate ode model with only dormancy specialist present
state["N1_activ"] <- 100 # initial density of dormancy strategist
state["N2"] <- 0 # initial density of `gleaner`

out <- ode(y = state, 
           times = times, 
           func = dorm_cr_ode, 
           parms = parameters,
           events = list(func = respulse, time = pulseseq))

# take a subset of dynamics for a single resource pulse cycle
dorm_df <- out %>% data.frame() %>% filter(time > 19200 & time < (19200 + pulsefreq)) 
dorm_df$Ntot <- dorm_df$N1_activ + dorm_df$N1_dorm
var(dorm_df$R)

df_gg <- data.frame(gleaner = glean_df$R, dormant = dorm_df$R) %>% 
  pivot_longer(cols = c("gleaner", "dormant"))

### Plot Fig 1D
p1_D <- ggplot(df_gg) + 
  geom_density(aes(x = value, col = name, fill = name),
               alpha = 0.5, linewidth = 1,
               adjust = 1/2) +
  scale_colour_manual(values = c("#E69F00", "#0072B2"))+ 
  scale_fill_manual(values = c("#E69F00", "#0072B2"))+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.position = "none") +
  coord_cartesian(expand = FALSE) + 
  xlab("Resource") + ylab("Density")

p1_D

###############################
