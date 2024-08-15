####################################################################
## dormancy-tradeoff-figS1.R: plot fig S1 ##########################
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


### New simulation function that accommodates all functional forms ------------------------------------------------------
dorm_cr_ode_allfunc <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dN1_activ <- N1_activ*((mu1*R)/(Ks1 + R) - d - ma) - 
      dorm_exp * dorm * (exp(-switch_dorm*R)) * N1_activ - 
      (1 - dorm_exp) * dorm * (1 - 0.25*switch_dorm*R) * N1_activ + 
      active_exp * dorm * (exp(switch_active*R) - 1) * N1_dorm + 
      (1 - active_exp) * dorm * (switch_active*R) * N1_dorm
    
    dN1_dorm <- dorm*(
      dorm_exp * (exp(-switch_dorm*R)) * N1_activ + 
      (1 - dorm_exp) * (1 - 0.25*switch_dorm*R) * N1_activ - 
      active_exp * (exp(switch_active*R) - 1) * N1_dorm - 
      (1 - active_exp) * (switch_active*R) * N1_dorm - 
        d*N1_dorm - 
        md*N1_dorm)
    
    dN2 <- N2*((mu2*R)/(Ks2 + R) - d - m) 
    dN3 <- N3*((mu3*R)/(Ks3 + R) - d - m) 
    dR <-  d*(R0 - R) - 
      N1_activ*(mu1*R*Qs1)/(Ks1 + R) - 
      N2*(mu2*R*Qs2)/(Ks2 + R) - 
      N3*(mu3*R*Qs3)/(Ks3 + R)
    return(list(c(dN1_activ, dN1_dorm, dN2, dN3, dR)))
  })
}


### Analytical calculation of per capita growth rate -------------------------------------------
Analytical <- function(resource, parameters){
  with(as.list(c(R = resource, parameters)), {
    Omega <- -(mu1 * R)/(Ks1 + R) + ma - md
    switchFunc_dorm <- dorm_exp * exp(-switch_dorm * R) + (1 - dorm_exp) * (1 - switch_dorm/4 * R)
    switchFunc_active <- active_exp * (exp(switch_active * R) - 1) + (1 - active_exp) * (switch_active * R)
    first.term <- -0.5 * (Omega + switchFunc_dorm + switchFunc_active + 2 * md)
    second.term <- 0.5 * sqrt((Omega + switchFunc_dorm + switchFunc_active)^2 - 4 * Omega * switchFunc_active)
    eigen.1 <- first.term - second.term
    eigen.2 <- first.term + second.term
    dominant <- ifelse(eigen.1 > eigen.2, eigen.1, eigen.2)  
    return(dominant)
  })
}


### Function to simulate competition and plot ------------------------------------------------------
GeneratePerCapita <- function(param, dorm_exponential_form, active_exponential_form){
  
  # Create used parameters 
  param.used <- c(param, dorm_exp = dorm_exponential_form, active_exp = active_exponential_form)
  
  # Simulate competition
  out <- ode(y = state, 
             times = times, 
             func = dorm_cr_ode_allfunc, 
             parms = param.used,
             events = list(func = respulse, time = pulseseq))
  
  # Wrangle data 
  dat <- out %>% data.frame() %>% 
    mutate(Rinflate = R*10) %>% 
    select(-c(R, N3)) %>% 
    pivot_longer(-time, names_to = "state_vars", values_to = "abund")
  dat$cons_res <- "consumers"
  dat$cons_res[dat$state_vars == "R"] <- "resource"
  
  # save resource distribution over a single cycle 
  respdf <- out %>% data.frame() %>% filter(time > 19200 & time < (19200 + pulsefreq)) 
  respdf$Ntot <- respdf$N1_activ + respdf$N1_dorm
  
  # Realised percapita growth response of dormancy strategist (dormant + active phase)
  percapNdorm = c()
  for (i in 2:length(respdf$Ntot)){
    percapNdorm[i-1] <- (log(respdf$Ntot[i]) - log(respdf$Ntot[i - 1]))/0.1
  }
  
  # Realised per capita growth response of gleaner 
  percap_reg = monod_func(param.used["mu2"], ks = param.used["Ks2"], respdf$R[-1]) - param.used["m"]
  percap_reg = c()
  for (i in 2:length(respdf$N2)){
    percap_reg[i-1] <- (log(respdf$N2[i]) - log(respdf$N2[i - 1]))/0.1
  }
  
  # Per capita growth response of active phase of dormancy strategist 
  percap_active = monod_func(param.used["mu1"], ks = param.used["Ks1"], respdf$R[-1]) - param.used["m"]
  percapdf <- data.frame(percap_dorm_integrated = percapNdorm, 
                         percap_reg = percap_reg, 
                         resconc = respdf$R[-1], 
                         percap_active = percap_active)
  percapgg <- pivot_longer(percapdf, -resconc)
  
  # Analytical calculation of dormant per capita growth rate
  resconc <- seq(0.00001, 8, l = 1000)
  result.analytical <- rep(0, length = 1000)
  for(i in 1:1000){
    result.analytical[i] <- Analytical(resconc[i], param.used)
  }
  test.analytical <- data.frame(resconc = resconc, value = result.analytical)
  
  # Plot 
  time_series <-
    ggplot(dat, aes(x = time, y = abund)) +
    geom_line(aes(col = state_vars, lty = state_vars), linewidth = 0.8) +
    coord_cartesian(expand = FALSE) +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          legend.position = "none") +
    scale_colour_manual(values = c("#E69F00", "#E69F00", "#0072B2", "#999999")) + 
    scale_linetype_manual(values = c("solid", "dotted", "solid", "solid")) + 
    xlim(9110, 9760) + 
    ylim(0, 500) + 
    ylab("Abundance") + xlab("Time")
  
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
    geom_line(data = test.analytical, aes(x = resconc, y = value), color = "black", size = 0.4)
  
  # func_resp_lowres <- 
  #   ggplot(percapgg, aes(x = resconc, y = value)) +
  #   geom_abline(slope = 0, intercept = 0, col = "grey") +
  #   geom_line(aes(col = name, alpha = name), linewidth = 0.7) +
  #   scale_colour_manual(values = c("#E69F00", "#E69F00", "#0072B2")) +
  #   scale_alpha_manual(values = c(1, 0.5, 1)) +
  #   theme(axis.text = element_text(size = 10),
  #         axis.title = element_text(size = 12),
  #         legend.position = "none") + 
  #   ylab("per capita growth rate") + xlab("Resource") + 
  #   xlim(0,0.4) + scale_y_continuous(limits = c(-0.025, 0.032), 
  #                                    breaks = c(-0.025, 0, 0.025)) + 
  #   geom_line(data = test.analytical, aes(x = resconc, y = value), color = "black", size = 0.4)
  
  # return(time_series + func_resp + func_resp_lowres)
  return(time_series + func_resp)
  
}


### Create plots -----------------------------------------------------------------------------------
DEAL <- GeneratePerCapita(parameters, dorm_exponential_form = 1, active_exponential_form = 0)
DEAE <- GeneratePerCapita(parameters, dorm_exponential_form = 1, active_exponential_form = 1)
DLAE <- GeneratePerCapita(parameters, dorm_exponential_form = 0, active_exponential_form = 1)
DLAL <- GeneratePerCapita(parameters, dorm_exponential_form = 0, active_exponential_form = 0)
Final_figure <- (DEAL / DEAE / DLAE / DLAL) 
Final_figure


