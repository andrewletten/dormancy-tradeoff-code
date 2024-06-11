## General functions

#' ODE model to simulate resource competition between a gleaner, opportunist 
#' dormancy specialists
dorm_cr_ode <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dN1_activ <- N1_activ*((mu1*R)/(Ks1 + R) - d - ma) - 
      dorm*(exp(-switch_dorm*R)*N1_activ - switch_active*N1_dorm*R)
    dN1_dorm <- dorm*(
      exp(-switch_dorm*R)*N1_activ - 
        switch_active*N1_dorm*R - 
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

#' Event function for resource pulsing
respulse <- function(Time, State, Pars) {
  with(as.list(State), {
    R <- R + pulsesize
    N1_activ <- N1_activ
    N1_dorm <- N1_dorm
    N2 <- N2
    N3 <- N3
    return(c(N1_activ, N1_dorm, N2, N3, R))
  })
}

#' Generate Monod function from Monod params and a vector of resource concs.
monod_func <- function(mu_max, ks, r_vec){
  (mu_max*r_vec)/(ks + r_vec)
}

#' Extract zipped simulations files and combine into a named list
extract_simfiles <- function(sppcombo, dirpath = "data/"){
  file_names <- list.files(path=paste0(dirpath, sppcombo, "_gz"), full.names=TRUE) # Change for pairwise combinations, e.g. "N1N2_gz")
  dat_list <- file_names %>% 
    map(~read_csv(.x, progress = FALSE, show_col_types = FALSE), 
        .progress = TRUE)
  names(dat_list) <- map_vec(file_names, 
                             function(x) {gsub("\\w+/\\w+/(\\w+).csv.gz", "\\1", x[1])[[1]][1]})
  return(dat_list)
  }

#' Extract resource size and supply frequency from named list
extract_sizefreq <- function(dat_list){
  size <- map_vec(names(dat_list), 
                  function(x) {gsub(".+p(\\d+)_.+", replacement = "\\1", x)}) %>% 
    as.numeric()
  freq <- map_vec(names(dat_list), 
                  function(x) {gsub(".+i(\\d+)", replacement = "\\1", x)}) %>% 
    as.numeric()
  list(size, freq)
}

#' Classify coexistence/exclusion outcomes from final resource pulse interval
coexist_or_not <- function(dat_list, respulsecombo){
  lastsamp <- list()
  for (i in 1:length(dat_list)){
    pulsefreq <- ifelse(respulsecombo$size[i] > 0, respulsecombo$freq[i], 1)
    totaltime <- nrow(dat_list[[i]])
    dat_list[[i]]$N1_tot <- dat_list[[i]]$N1_activ + dat_list[[i]]$N1_dorm
    lastsamp[[i]] <- data.frame(
      dat_list[[i]][((totaltime) - (pulsefreq)):(totaltime),] > 0.0000001
    ) %>% 
      summarize(across(c(N1_tot, N2, N3), all)) 
    print(paste("sim", i))
  }
  return(lastsamp)
}

#' Prepare summary data for plotting
summarise_outcomes <- function(dat_list, dorm_cost = NULL){
  sizefreq <- extract_sizefreq(dat_list)
  respulsecombo <- data.frame(size = sizefreq[[1]], freq = sizefreq[[2]])
  lastsamp <- coexist_or_not(dat_list, respulsecombo)
  lastsamp.df <- do.call("rbind", lastsamp)
  if(is.null(dorm_cost)){
    lastsamp.df <- cbind(lastsamp.df, respulsecombo)
  } else {
    lastsamp.df <- cbind(lastsamp.df, dorm_cost)
  }
  lastsamp.df$combos <- do.call(paste, c(lastsamp.df[,1:3], sep = "-"))
  lookup <- data.frame(combos = c("FALSE-TRUE-FALSE",
                                  "FALSE-TRUE-TRUE",
                                  "FALSE-FALSE-TRUE",
                                  "TRUE-TRUE-TRUE", 
                                  "TRUE-FALSE-TRUE",
                                  "TRUE-TRUE-FALSE",
                                  "TRUE-FALSE-FALSE",
                                  "FALSE-FALSE-FALSE"),
                       spp = as.factor(c("Gl", 
                                         "Gl-Op", 
                                         "Op",
                                         "Gl-Op-Do", 
                                         "Op-Do", 
                                         "Gl-Do", 
                                         "Do",
                                         "none"
                       )),
                       colscheme = c("#0072B2", 
                                     "#D55E00", 
                                     "#CC79A7",
                                     "#F0E442",
                                     "#009E73", 
                                     "#56B4E9",
                                     "#E69F00",
                                     "#999999"))
  
  lastsamp_comb <- lastsamp.df %>% 
    merge(lookup, by.y = "combos") %>% 
    mutate(spp = fct_relevel(spp, 
                             "Gl", "Op", "Do", 
                             "Gl-Op", "Gl-Do", "Op-Do", 
                             "Gl-Op-Do", "none")) %>% 
    arrange(spp) 
  if(is.null(dorm_cost)){lastsamp_comb <- lastsamp_comb %>% 
    mutate(freq = replace(freq, freq == 0, "cont.")) %>% 
    mutate(freq = fct_relevel(freq, 
                              "cont.", "24", "48", 
                              "72", "96", "120", 
                              "144", "168", "192",
                              "216", "240", "264", 
                              "288", "312", "336"))}
  return(lastsamp_comb)
}


#' Plot heat map of competition outcomes under different patterns of resource pulsing
plot_comp_resource <- function(dat_gg, outline = NULL){
  unicols <- unique(dat_gg$colscheme)
  if(!is.null(outline)){
    outline2 <- merge(outline, dat_gg, by = c("size", "freq"))
  } else {}
  ggplot(dat_gg, aes(x = size, y = as.factor(freq))) + 
    geom_tile(aes(fill = spp)) +
    scale_fill_manual(values = unicols) + 
    theme_cowplot() +
    theme(legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.text = element_text(size=10),
          legend.position="bottom",
          legend.justification = "center") +
    guides(fill = guide_legend(nrow = 1, 
                               title="", 
                               label.position = "bottom")) +
    xlab('Resource pulse size') + ylab('Resource pulse interval') +
    coord_cartesian(expand = FALSE) +
    {if(!is.null(outline))
      geom_tile(data = outline2, colour = "grey", alpha = 0, linewidth = 1, linejoin = "round")} +
    {if(!is.null(outline))
      geom_tile(data = outline2, aes(fill = spp.y), alpha = 1)} +
    theme(axis.text = element_text(size = 8),
          axis.title= element_text(size = 12))
}

#' Check dynamics are stationary
stat_dynamics <- function(dat_list){
  datstat <- list()
  sizefreq <- extract_sizefreq(dat_list)
  respulsecombo <- data.frame(size = sizefreq[[1]], freq = sizefreq[[2]])
  for (i in 1:length(dat_list)){
    pulsefreq <- ifelse(respulsecombo[i,2] > 0, respulsecombo[i,2], 1)
    totaltime <- nrow(dat_list[[i]])
    datstat[[i]] <- tail(dat_list[[i]][seq(100, totaltime, pulsefreq),])
  }
  names(datstat) <- names(dat_list)
  return(datstat)
}
