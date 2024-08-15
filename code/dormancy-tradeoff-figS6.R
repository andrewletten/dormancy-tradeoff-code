####################################################################
## dormancy-tradeoff-figS6.R: plot fig S6 (1c analytically) ########
## Letten, Yamamichi, Richardson and Ke ############################
## Microbial dormancy supports multi-species coexistence ###########
## under resource fluctuations #####################################
####################################################################

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


### Approach 1: analytical solution -------------------------------------------
Analytical <- function(resource, parameters){
  with(as.list(c(R = resource, parameters)), {
    Omega <- -(mu1 * R)/(Ks1 + R) + ma - md
    switchFunc_dorm <- exp(-switch_dorm * R)
    switchFunc_active <- switch_active * R
    first.term <- -0.5 * (Omega + switchFunc_dorm + switchFunc_active + 2 * md)
    second.term <- 0.5 * sqrt((Omega + switchFunc_dorm + switchFunc_active)^2 - 4 * Omega * switchFunc_active)
    eigen.1 <- first.term - second.term
    eigen.2 <- first.term + second.term
    dominant <- ifelse(eigen.1 > eigen.2, eigen.1, eigen.2)  
    return(dominant)
  })
}

resconc <- seq(0.00001, 8, l = 1000)
result.analytical <- rep(0, length = 1000)
for(i in 1:1000){
  result.analytical[i] <- Analytical(resconc[i], parameters)
}
test.analytical <- data.frame(resconc = resconc, value = result.analytical)


# ### Approach 2: numerical timescale separation --------------------------------
# PerCapita <- function(resource, parameters){
#   with(as.list(c(R = resource, parameters)), {
#     A.phi <- -(mu1 * R)/(Ks1 + R) + ma - md
#     B.phi <-  (mu1 * R)/(Ks1 + R) - ma + md - exp(-switch_dorm * R) - switch_active * R
#     C.phi <- switch_active * R
#     phi.star.1 <- (-B.phi - sqrt(B.phi^2 - 4 * A.phi * C.phi)) / (2 * A.phi)
#     phi.star.2 <- (-B.phi + sqrt(B.phi^2 - 4 * A.phi * C.phi)) / (2 * A.phi)
#     criterion.1 <- (phi.star.1 >= 0) && (phi.star.1 <= 1)
#     criterion.2 <- (phi.star.2 >= 0) && (phi.star.2 <= 1)
#     phi.star <- ifelse(criterion.1 == TRUE, phi.star.1, 
#                        ifelse(criterion.2 == TRUE, phi.star.2, NA))  
#     per.capita <- phi.star * (mu1 * R)/(Ks1 + R) - phi.star * (ma - md) - md
#     return(per.capita)
#   })
# }
# 
# resconc <- seq(0.00001, 8, l = 1000)
# result.percapita <- rep(0, length = 1000)
# for(i in 1:1000){
#   result.percapita[i] <- PerCapita(resconc[i], parameters)
# }
# test.percapita <- data.frame(resconc = resconc, value = result.percapita)


# ### Approach 3: numerical dominant eigenvalue ---------------------------------
# A <- function(resource, parameters){
#   with(as.list(c(R = resource, parameters)), {
#     A11 <- (mu1 * R/(Ks1 + R)) - ma - exp(-switch_dorm * R)
#     A12 <- switch_active * R
#     A21 <- exp(-switch_dorm * R)
#     A22 <- -md - switch_active * R
#     mat <- matrix(c(A11, A12, A21, A22), 2, 2, byrow = T)
#     return(mat)
#   })
# }
# resconc <- seq(0.00001, 8, l = 1000)
# result <- rep(0, length = 1000)
# for(i in 1:1000){
#   Temp <- A(resconc[i], parameters)
#   result[i] <- max(eigen(Temp)$value) 
# }
# test.eigen <- data.frame(resconc = resconc, value = result)



### Plot ----------------------------------------------------------------------
func_resp <- 
  ggplot(percapgg, aes(x = resconc, y = value)) +
  geom_abline(slope = 0, intercept = 0, col = "grey") +
  geom_line(aes(col = name, alpha = name), linewidth = 1) +
  scale_colour_manual(values = c("#E69F00", "#E69F00", "#0072B2")) +
  scale_alpha_manual(values = c(1, 0.5, 1)) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.position = "none") + 
  ylab("per capita growth rate") + xlab("Resource") + 
  xlim(0,8) + 
  geom_line(data = test.analytical, aes(x = resconc, y = value), color = "black", size = 0.4) # + 
  # geom_line(data = test.percapita, aes(x = resconc, y = value), color = "red", size = 1.0, linetype = 3)


func_resp_lowres <- 
  ggplot(percapgg, aes(x = resconc, y = value)) +
  geom_abline(slope = 0, intercept = 0, col = "grey") +
  geom_line(aes(col = name, alpha = name), linewidth = 0.7) +
  scale_colour_manual(values = c("#E69F00", "#E69F00", "#0072B2")) +
  scale_alpha_manual(values = c(1, 0.5, 1)) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.position = "none") + 
  ylab("per capita growth rate") + xlab("Resource") + 
  xlim(0,0.4) + scale_y_continuous(limits = c(-0.025, 0.032), 
                                   breaks = c(-0.025, 0, 0.025)) + 
  geom_line(data = test.analytical, aes(x = resconc, y = value), color = "black", size = 0.4) # + 
  # geom_line(data = test.percapita, aes(x = resconc, y = value), color = "red", size = 1.0, linetype = 3)


func_resp + func_resp_lowres


