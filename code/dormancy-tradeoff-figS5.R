############################################################################
## dormancy-tradeoff-figS5.R: plot fig S5 (fig 1b with floquet theory) #####
## Letten, Yamamichi, Richardson and Ke ####################################
## Microbial dormancy supports multi-species coexistence ###################
## under resource fluctuations #############################################
############################################################################

### Load required packages and project code -----------------------------------
library(deSolve)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(patchwork)
theme_set(theme_cowplot()) 
source("code/dormancy-tradeoff-functions.R")
source("code/dormancy-tradeoff-spp-params.R")


### Common parameter setting --------------------------------------------------
totaltime <- 100000 # total simulation time
interval <- 0.1
times <- seq(0, totaltime, by = interval)
pulsefreq <- 240 # resource pulse interval
pulsesize <- 8 # resource pulse magnitude
pulseseq <- round(seq(pulsefreq, totaltime, by = pulsefreq), 1)


### Monoculture initial conditions --------------------------------------------
state.gleaner.mono <- c(N1_activ = 0, # initial density of active phase of dormant strategist
                        N1_dorm = 0, # initial density of dormant phase of dormant strategist
                        N2 = 100, # initial density of `gleaner`
                        N3 = 0, # initial density of `opportunist` (ignored for initial analysis [Fig 1])
                        R = 10)

state.dorm.mono <- c(N1_activ = 100, # initial density of active phase of dormant strategist
                     N1_dorm = 0, # initial density of dormant phase of dormant strategist
                     N2 = 0, # initial density of `gleaner`
                     N3 = 0, # initial density of `opportunist` (ignored for initial analysis [Fig 1])
                     R = 10)


### Common function setting ---------------------------------------------------
# Invasion of the dormant species under gleaner monoculture 
only_dorm <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dN1_activ <- N1_activ*((mu1*R)/(Ks1 + R) - d - ma) - 
      dorm*(exp(-switch_dorm*R)*N1_activ - switch_active*N1_dorm*R)
    dN1_dorm <- dorm*(
      exp(-switch_dorm*R)*N1_activ - 
        switch_active*N1_dorm*R - 
        d*N1_dorm - 
        md*N1_dorm)
    dR <- 0
    return(list(c(dN1_activ, dN1_dorm, dR)))
  })
}

# Pulse function using gleaner monoculture resource cycles 
# Note here that "Resource.gleaner.mono" is sourced from global environment 
respulse_only_dorm <- function(Time, State, Pars) {
  with(as.list(State), {
    N1_activ <- N1_activ
    N1_dorm <- N1_dorm
    R <- Resource.gleaner.mono[Time/interval+1]
    return(c(N1_activ, N1_dorm, R))
  })
}

# Invasion of the gleaner species under dormant monoculture 
only_gleaner <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dN2 <- N2*((mu2*R)/(Ks2 + R) - d - m) #population growth integral (Klausmeier, 2008)
    # dN2 <- ((mu2*R)/(Ks2 + R) - d - m) #per capita growth integral (Xiao & Fussmann, 2013)
    dR <- 0
    return(list(c(dN2, dR)))
  })
}

# Pulse function using dormant monoculture resource cycles 
# Note here that "Resource.dorm.mono" is sourced from global environment 
respulse_only_gleaner <- function(Time, State, Pars) {
  with(as.list(State), {
    N2 <- N2
    R <- Resource.dorm.mono[Time/interval+1]
    return(c(N2, R))
  })
}


### Floquet modifier for dormant ----------------------------------------------
### "x" is the numerical input of the target parameter (for uniroot)
### "parms" is a vector of all parameters used for the simulation
### "target" is a character indicating the name of the target parameter 
Floquet_dorm <- function(x, parms, target){
  
  # Set parameters
  parms[target] <- x
  parameters_only_dorm <- parms
  
  # # Recalculate gleaner monoculture with new parameters
  # # This IS NOT necessary because Fig.1B are parameters for the dormant strategy
  # out <- ode(y = state.gleaner.mono, 
  #            times = times, 
  #            func = dorm_cr_ode, 
  #            parms = parameters_only_dorm,
  #            events = list(func = respulse, time = pulseseq))
  # 
  # # Export new gleaner monoculture out of function memory
  # # There must be better ways (e.g., assign()) but for now let's just to it this way
  # # This is necessary because later functions ("respulse_only_dorm") require the new "Resource.gleaner.mono" from the global environment
  # Resource.gleaner.mono <<- 
  #   out %>% 
  #   data.frame() %>% 
  #   filter(time > 9120 & time <= (9120 + pulsefreq)) %>%
  #   select(R) %>%
  #   unlist()
  
  # Set initial conditions
  # Note that "Resource.gleaner.mono" is sourced from global environment so this needs to be calculated outside first
  state_only_dorm <- cbind(c(N1_activ = 1, N1_dorm = 0, R = unname(Resource.gleaner.mono[1])),
                           c(N1_activ = 0, N1_dorm = 1, R = unname(Resource.gleaner.mono[1])))
  times_only_dorm <- seq(0.1, pulsefreq, by = interval)
  
  # Perform one-cycle invasion of dormant to construct fundamental matrix
  Fundamental <- matrix(0, 2, 2)
  for(i in 1:ncol(Fundamental)){
    temp <- ode(y = state_only_dorm[, i],
                times = times_only_dorm,
                func = only_dorm,
                parms = parameters_only_dorm,
                events = list(func = respulse_only_dorm, time = times_only_dorm))
    Fundamental[, i] <- as.vector(temp[nrow(temp), c("N1_activ", "N1_dorm")])
  }
  
  # Calculate eigenvalue for fundamental matrix
  multiplier <- max(Re(eigen(Fundamental)$value))
  
  # Return value to root.solve boundary = 1
  return(multiplier - 1)
  
}


### Floquet modifier for gleaner ----------------------------------------------
### "x" is the numerical input of the target parameter (for uniroot)
### "parms" is a vector of all parameters used for the simulation
### "target" is a character indicating the name of the target parameter 
Floquet_gleaner <- function(x, parms, target){
  
  # Set parameters
  parms[target] <- x
  parameters_only_gleaner <- parms
  
  # Recalculate dormant monoculture with new parameters
  # This IS necessary because Fig.1B are parameters for the dormant strategy
  out <- ode(y = state.dorm.mono, 
             times = times, 
             func = dorm_cr_ode, 
             parms = parameters_only_gleaner,
             events = list(func = respulse, time = pulseseq))
  
  # Export new dormant monoculture out of function memory
  # There must be better ways (e.g., assign()) but for now let's just to it this way
  # This is necessary because later functions ("respulse_only_gleaner") require the new "Resource.dorm.mono" from the global environment
  Resource.dorm.mono <<- 
    out %>% 
    data.frame() %>% 
    filter(time > 99120 & time <= (99120 + pulsefreq)) %>%
    select(R) %>%
    unlist()
  
  # Set initial conditions
  state_only_gleaner <- cbind(c(N2 = 1, R = unname(Resource.dorm.mono[1])))
  times_only_gleaner <- seq(0.1, pulsefreq, by = interval)
  
  # Perform one-cycle invasion of dormant to construct fundamental matrix
  Fundamental <- matrix(0, 1, 1)
  for(i in 1:ncol(Fundamental)){
    temp <- ode(y = state_only_gleaner[, i], 
                times = times_only_gleaner, 
                func = only_gleaner, 
                parms = parameters_only_gleaner,     
                events = list(func = respulse_only_gleaner, time = times_only_gleaner))
    Fundamental[, i] <- as.vector(temp[nrow(temp), c("N2")])
  }
  
  # Calculate eigenvalue for fundamental matrix
  multiplier <- max(Re(eigen(Fundamental)$value))
  
  # Return value to root.solve boundary = 1
  return(multiplier - 1) #population growth baseline = 1 (Klausmeier, 2008)
  # return(multiplier - 0) #per capita growth baseline = 0 (Xiao & Fussmann, 2013)
  
}


### ---------------------------------------------------------------------------
### Floquet multiplier version of Fig.1B --------------------------------------
# Note parameter change for Fig.1B 
source("code/dormancy-tradeoff-spp-params.R")
# parameters["mu1"] <- 0.07

# Simulation setting: target parameter (y-axis) "ks1" across sequence parameter (x-axis) "md"
ks_seq <- 0.05*2^(seq(0, 8, 0.05)) # Ks1 sequence (active phase)
md_seq <- 0.025/(2^seq(0, 9, by = 0.05)) # md sequence (dormant phase)
target.parameter <- "Ks1"
target.seq <- ks_seq
sequence.parameter <- "md"
sequence.seq <- md_seq

# Data storage space
result.dormant <- rep(0, length(md_seq))
result.gleaner <- rep(0, length(md_seq))

# Get gleaner monoculture (to calculate Floquet multiplier for dormant)
out.gleaner.mono <- ode(y = state.gleaner.mono, 
                        times = times, 
                        func = dorm_cr_ode, 
                        parms = parameters,
                        events = list(func = respulse, time = pulseseq))

Resource.gleaner.mono <- 
  out.gleaner.mono %>% 
  data.frame() %>% 
  filter(time > 99120 & time <= (99120 + pulsefreq)) %>%
  select(R) %>%
  unlist()

# Calculate Floquet multiplier for dormant
for (i in 1:length(md_seq)){
  
  tryCatch({
    temp.parameter <- parameters
    temp.parameter[sequence.parameter] <- sequence.seq[i]
    result.dormant[i] <- uniroot(Floquet_dorm, 
                                 interval = c(target.seq[1], 2 * target.seq[length(target.seq)]), 
                                 parms = temp.parameter, target = target.parameter)$root
    }, error=function(e){cat("ERROR :",conditionMessage(e), ", same invasibility across parameter range\n")})

}

# Aggregate data for dormant
Data.dormant <- data.frame(md = md_seq, Ks1 = result.dormant) 

# # Quickly plot invasion boundary for dormant
# ggplot(Data.dormant, aes(x = md, y = Ks1)) +
#   geom_line() +
#   scale_x_continuous(limits = c(min(sequence.seq), max(sequence.seq)), expand = c(0, 0), trans='log2') +
#   scale_y_continuous(limits = c(min(target.seq), max(target.seq)), expand = c(0, 0), trans='log2')

# Calculate Floquet multiplier for gleaner
for (i in 1:length(md_seq)){
  
  tryCatch({
    temp.parameter <- parameters
    temp.parameter[sequence.parameter] <- sequence.seq[i]
    result.gleaner[i] <- uniroot(Floquet_gleaner, 
                                 interval = c(target.seq[1], 2 * target.seq[length(target.seq)]), 
                                 parms = temp.parameter, target = target.parameter)$root
  }, error=function(e){cat("ERROR :",conditionMessage(e), ", same invasibility across parameter range\n")})
  
}

# Aggregate data for gleaner
Data.gleaner <- data.frame(md = md_seq, Ks1 = result.gleaner) 

# # Quickly plot invasion boundary for gleaner
# ggplot(Data.gleaner, aes(x = md, y = Ks1)) +
#   geom_line() +
#   scale_x_continuous(limits = c(min(sequence.seq), max(sequence.seq)), expand = c(0, 0), trans='log2') +
#   scale_y_continuous(limits = c(min(target.seq), max(target.seq)), expand = c(0, 0), trans='log2')

### Plot the two invasion boundaries (adjust boundary for plotting to match main Fig.1B) 
### geom_ribbon version of plotting 
x.ticks.1B <- 0.025/(2^seq(0, 8, by = 1))
y.ticks.1B <- 0.05*2^(seq(0, 8, 1))
# saveRDS(Data.dormant, "Data_dormant_Fig1B.rds") # Commented out after running simulations
# saveRDS(Data.gleaner, "Data_gleaner_Fig1B.rds") # Commented out after running simulations
Data_dormant_Fig1B <- readRDS("Data_dormant_Fig1B.rds")
Data_gleaner_Fig1B <- readRDS("Data_gleaner_Fig1B.rds")
names(Data_dormant_Fig1B) <- c("Sequence", "Boundary")
names(Data_gleaner_Fig1B) <- c("Sequence", "Boundary")
target.seq.1B <- 0.05*2^(seq(0, 8, 0.05)) # Ks1 sequence (active phase)
sequence.seq.1B <- 0.025/(2^seq(0, 8, by = 0.05)) # md sequence (dormant phase)

Fig1B <- 
  ggplot() + 
  geom_ribbon(data = Data_dormant_Fig1B[Data_dormant_Fig1B$Boundary > 0.0, ], # Remove region where root is out of bound 
              aes(x = Sequence, 
                  ymin = min(target.seq.1B), 
                  ymax = Boundary), 
              fill = "#56B4E9") +  
  geom_ribbon(data = Data_gleaner_Fig1B[Data_gleaner_Fig1B$Boundary > 0.0, ], # Remove region where root is out of bound
              aes(x = Sequence, 
                  ymin = min(target.seq.1B), 
                  ymax = Boundary), 
              fill = "#E69F00") +
  scale_x_continuous(limits = c(min(sequence.seq.1B), max(sequence.seq.1B)), 
                     expand = c(0, 0), 
                     trans='log2', 
                     breaks = x.ticks.1B, 
                     labels = as.character(round(x.ticks.1B * 1000, 2))) + 
  scale_y_continuous(limits = c(min(target.seq.1B), max(target.seq.1B)), 
                     expand = c(0, 0), 
                     trans='log2', 
                     breaks = y.ticks.1B, 
                     labels = as.character(y.ticks.1B)) + 
  theme_cowplot() +
  theme(legend.key.size = unit(0.5, 'cm'), 
        legend.text = element_text(size=10),
        legend.justification = "center") + 
  xlab('Mortality rate (dormant phase)') + 
  ylab('Monod half saturation') +
  coord_cartesian(expand = FALSE) +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.position = "none",
        panel.background = element_rect(fill = "#0072B2")) 
# ggsave("Fig1B_Floquet.pdf", width = 6, height = 6)


### ---------------------------------------------------------------------------
### Floquet multiplier version of Fig.S2 --------------------------------------
# Note parameter change for Fig.S2 
source("code/dormancy-tradeoff-spp-params.R")
# parameters["mu1"] <- 0.07
# parameters["Ks1"] <- 0.1
# parameters["md"] <- 0.000781250

# Simulation setting: target parameter (y-axis) "lambda" across sequence parameter (x-axis) "md"
lambda_seq <- 0.0125*2^(seq(0, 7, by = 0.05)) # lambda sequence (active phase)
md_seq <- 0.00625/(2^seq(0, 8, by = 0.05)) # md sequence (dormant phase)
target.parameter <- "switch_dorm"
target.seq <- lambda_seq
sequence.parameter <- "md"
sequence.seq <- md_seq

# Data storage space
result.dormant <- rep(0, length(md_seq))
result.gleaner <- rep(0, length(md_seq))

# Get gleaner monoculture (to calculate Floquet multiplier for dormant)
out.gleaner.mono <- ode(y = state.gleaner.mono, 
                        times = times, 
                        func = dorm_cr_ode, 
                        parms = parameters,
                        events = list(func = respulse, time = pulseseq))

Resource.gleaner.mono <- 
  out.gleaner.mono %>% 
  data.frame() %>% 
  filter(time > 99120 & time <= (99120 + pulsefreq)) %>%
  select(R) %>%
  unlist()

# Calculate Floquet multiplier for dormant
for (i in 1:length(md_seq)){
  
  tryCatch({
    temp.parameter <- parameters
    temp.parameter[sequence.parameter] <- sequence.seq[i]
    result.dormant[i] <- uniroot(Floquet_dorm, 
                                 interval = c(target.seq[1], 16 * target.seq[length(target.seq)]), 
                                 parms = temp.parameter, target = target.parameter)$root
  }, error=function(e){cat("ERROR :",conditionMessage(e), ", same invasibility across parameter range\n")})
  
}

# Aggregate data for dormant
Data.dormant <- data.frame(md = md_seq, lambda = result.dormant) 

# # Quickly plot invasion boundary for dormant
# ggplot(Data.dormant, aes(x = md, y = lambda)) +
#   geom_line() +
#   scale_x_continuous(limits = c(min(sequence.seq), max(sequence.seq)), expand = c(0, 0), trans='log2') +
#   scale_y_continuous(limits = c(min(target.seq), max(target.seq)), expand = c(0, 0), trans='log2')

# Calculate Floquet multiplier for gleaner
for (i in 1:length(md_seq)){
  
  tryCatch({
    temp.parameter <- parameters
    temp.parameter[sequence.parameter] <- sequence.seq[i]
    result.gleaner[i] <- uniroot(Floquet_gleaner, 
                                 interval = c(target.seq[1], 16 * target.seq[length(target.seq)]), 
                                 parms = temp.parameter, target = target.parameter)$root
  }, error=function(e){cat("ERROR :",conditionMessage(e), ", same invasibility across parameter range\n")})
  
}

# Aggregate data for gleaner
Data.gleaner <- data.frame(md = md_seq, lambda = result.gleaner) 

# # Quickly plot invasion boundary for gleaner
# ggplot(Data.gleaner, aes(x = md, y = lambda)) +
#   geom_line() +
#   scale_x_continuous(limits = c(min(sequence.seq), max(sequence.seq)), expand = c(0, 0), trans='log2') +
#   scale_y_continuous(limits = c(min(target.seq), max(target.seq)), expand = c(0, 0), trans='log2')

### Plot the two invasion boundaries (adjust boundary for plotting to match main Fig.S2) 
### geom_ribbon version of plotting 
x.ticks.S2 <- 0.00625/(2^seq(0, 8, by = 1))
y.ticks.S2 <- 0.0125*2^(seq(0, 7, by = 1))
# saveRDS(Data.dormant, "Data_dormant_FigS2.rds") # Commented out after running simulations 
# saveRDS(Data.gleaner, "Data_gleaner_FigS2.rds") # Commented out after running simulations
Data_dormant_FigS2 <- readRDS("Data_dormant_FigS2.rds")
Data_gleaner_FigS2 <- readRDS("Data_gleaner_FigS2.rds")
names(Data_dormant_FigS2) <- c("Sequence", "Boundary")
names(Data_gleaner_FigS2) <- c("Sequence", "Boundary")
target.seq.S2 <- 0.0125*2^(seq(0, 7, by = 0.05)) # lambda sequence (active phase)
sequence.seq.S2 <- 0.00625/(2^seq(0, 8, by = 0.05)) # md sequence (dormant phase)

FigS2 <- 
  ggplot() + 
  geom_ribbon(data = Data_dormant_FigS2, 
              aes(x = Sequence, 
                  ymin = Boundary, 
                  ymax = max(target.seq.S2)), 
              fill = "#56B4E9") +  
  geom_ribbon(data = Data_gleaner_FigS2[Data_gleaner_FigS2$Boundary > 0.0, ], # Remove region where root is out of bound
              aes(x = Sequence,
                  ymin = Boundary,
                  ymax = max(target.seq.S2)),
              fill = "#E69F00") +
  scale_x_continuous(limits = c(min(sequence.seq.S2), max(sequence.seq.S2)), 
                     expand = c(0, 0), 
                     trans='log2', 
                     breaks = x.ticks.S2, 
                     labels = as.character(round(x.ticks.S2 * 1000, 2))) + 
  scale_y_continuous(limits = c(min(target.seq.S2), max(target.seq.S2)),
                     expand = c(0, 0), 
                     trans='log2', 
                     breaks = y.ticks.S2, 
                     labels = as.character(y.ticks.S2)) + 
  theme_cowplot() +
  theme(legend.key.size = unit(0.5, 'cm'), 
        legend.text = element_text(size=10),
        legend.justification = "center") + 
  xlab('Mortality rate (dormant phase)') + 
  ylab(expression(paste('Switch dormant (', lambda, ')'))) + 
  coord_cartesian(expand = FALSE) +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.position = "none",
        panel.background = element_rect(fill = "#0072B2")) 
# ggsave("FigS2_Floquet.pdf", width = 6, height = 6)


### ---------------------------------------------------------------------------
### Floquet multiplier version of Fig.S3 --------------------------------------
# Note parameter change for Fig.S3 
source("code/dormancy-tradeoff-spp-params.R")
# parameters["mu1"] <- 0.07
parameters["Ks1"] <- 0.1

# Simulation setting: target parameter (y-axis) "lambda" across sequence parameter (x-axis) "gamma"
lambda_seq <- 0.0125*2^(seq(0, 7, by = 0.05)) # lambda sequence (active phase)
gamma_seq <- 0.01*(2^seq(0, 7, by = 0.05)) # gamma sequence (dormant phase)
target.parameter <- "switch_dorm"
target.seq <- lambda_seq
sequence.parameter <- "switch_active"
sequence.seq <- gamma_seq

# Data storage space
result.dormant <- rep(0, length(gamma_seq))
result.gleaner <- rep(0, length(gamma_seq))

# Get gleaner monoculture (to calculate Floquet multiplier for dormant)
out.gleaner.mono <- ode(y = state.gleaner.mono, 
                        times = times, 
                        func = dorm_cr_ode, 
                        parms = parameters,
                        events = list(func = respulse, time = pulseseq))

Resource.gleaner.mono <- 
  out.gleaner.mono %>% 
  data.frame() %>% 
  filter(time > 99120 & time <= (99120 + pulsefreq)) %>%
  select(R) %>%
  unlist()

# Calculate Floquet multiplier for dormant
for (i in 1:length(gamma_seq)){
  
  tryCatch({
    temp.parameter <- parameters
    temp.parameter[sequence.parameter] <- gamma_seq[i]
    result.dormant[i] <- uniroot(Floquet_dorm, 
                                 interval = c(target.seq[1], 2 * target.seq[length(target.seq)]), 
                                 parms = temp.parameter, target = target.parameter)$root
  }, error=function(e){cat("ERROR :",conditionMessage(e), ", same invasibility across parameter range\n")})
  
}

# Aggregate data for dormant
Data.dormant <- data.frame(gamma = gamma_seq, lambda = result.dormant) 

# # Quickly plot invasion boundary for dormant
# ggplot(Data.dormant, aes(x = gamma, y = lambda)) +
#   geom_line() +
#   scale_x_continuous(limits = c(min(sequence.seq), max(sequence.seq)), expand = c(0, 0), trans='log2') +
#   scale_y_continuous(limits = c(min(target.seq), max(target.seq)), expand = c(0, 0), trans='log2')

# Calculate Floquet multiplier for gleaner
for (i in 1:length(gamma_seq)){
  
  tryCatch({
    temp.parameter <- parameters
    temp.parameter[sequence.parameter] <- gamma_seq[i]
    result.gleaner[i] <- uniroot(Floquet_gleaner, 
                                 interval = c(target.seq[1], 2 * target.seq[length(target.seq)]), 
                                 parms = temp.parameter, target = target.parameter)$root
  }, error=function(e){cat("ERROR :",conditionMessage(e), ", same invasibility across parameter range\n")})
  
}

# Aggregate data for gleaner
Data.gleaner <- data.frame(gamma = gamma_seq, lambda = result.gleaner) 

# # Quickly plot invasion boundary for gleaner
# ggplot(Data.gleaner, aes(x = gamma, y = lambda)) +
#   geom_line() +
#   scale_x_continuous(limits = c(min(sequence.seq), max(sequence.seq)), expand = c(0, 0), trans='log2') +
#   scale_y_continuous(limits = c(min(target.seq), max(target.seq)), expand = c(0, 0), trans='log2')

### Plot the two invasion boundaries (adjust boundary for plotting to match main Fig.S3) 
### geom_ribbon version of plotting 
x.ticks.S3 <- 0.01*(2^seq(0, 7, by = 1))
y.ticks.S3 <- 0.0125*2^(seq(0, 7, by = 1))
# saveRDS(Data.dormant, "Data_dormant_FigS3.rds") # Commented out after running simulations
# saveRDS(Data.gleaner, "Data_gleaner_FigS3.rds") # Commented out after running simulations
Data_dormant_FigS3 <- readRDS("Data_dormant_FigS3.rds")
Data_gleaner_FigS3 <- readRDS("Data_gleaner_FigS3.rds")
names(Data_dormant_FigS3) <- c("Sequence", "Boundary")
names(Data_gleaner_FigS3) <- c("Sequence", "Boundary")
# Simulation setting: target parameter (y-axis) "lambda" across sequence parameter (x-axis) "gamma"
target.seq.S3 <- 0.0125*2^(seq(0, 7, by = 0.05)) # lambda sequence (active phase)
sequence.seq.S3 <- 0.01*(2^seq(0, 7, by = 0.05)) # gamma sequence (dormant phase)

FigS3 <- 
  ggplot() + 
  geom_ribbon(data = Data_dormant_FigS3, 
              aes(x = Sequence, 
                  ymin = Boundary, 
                  ymax = max(target.seq.S3)), 
              fill = "#56B4E9") +  
  geom_ribbon(data = Data_gleaner_FigS3[Data_gleaner_FigS3$Sequence > 0.088, ], # Remove region where root is out of bound
              aes(x = Sequence,
                  ymin = Boundary,
                  ymax = max(target.seq.S3)),
              fill = "#E69F00") +
  scale_x_continuous(limits = c(min(sequence.seq.S3), max(sequence.seq.S3)), 
                     expand = c(0, 0), 
                     trans='log2', 
                     breaks = x.ticks.S3, 
                     labels = as.character(round(x.ticks.S3, 2))) + 
  scale_y_continuous(limits = c(min(target.seq.S3), max(target.seq.S3)), 
                     expand = c(0, 0), 
                     trans='log2', 
                     breaks = y.ticks.S3, 
                     labels = as.character(y.ticks.S3)) + 
  theme_cowplot() +
  theme(legend.key.size = unit(0.5, 'cm'), 
        legend.text = element_text(size=10),
        legend.justification = "center") + 
  xlab(expression(paste('Switch active (', gamma, ')'))) + 
  ylab(expression(paste('Switch dormant (', lambda, ')'))) + 
  coord_cartesian(expand = FALSE) +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.position = "none",
        panel.background = element_rect(fill = "#0072B2")) 
# ggsave("FigS3_Floquet.pdf", width = 6, height = 6)


### ---------------------------------------------------------------------------
### Combine all figures -------------------------------------------------------
# Fig1B / FigS2 / FigS3 + plot_annotation(tag_levels = 'A')
# ggsave("FigS6_Floquet_Vertical.pdf", width = 4.0, height = 10)
# Fig1B + FigS2 + FigS3 + plot_annotation(tag_levels = 'A')
# ggsave("FigS6_Floquet_Horizontal.pdf", width = 10, height = 3.5)
