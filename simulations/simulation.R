#----------------------------------------------
#Simulation study to assess the GxE interactions and correlations
#Scenarios 1-6
#Feb 12, 2023
#----------------------------------------------


#----------------------------------------------
#Scenario 1: without strata, int1=0, int2=0, assess type 1 error
#simulation work for stratified model
#Nov 7, 2022
#----------------------------------------------


library(caseControlGE)
path <- "/users/zwang4/GXE/simulation_results_11072022/"
source("/users/zwang4/GXE/sim_data_function.R")
source("/users/zwang4/GXE/RetroGE.R")


#look at the population disease prevalence
dat = simFit(
  ncontrol = 500,
  ncase = 500,
  beta0 = -4.8,
  betaG_normPRS = 0.450,
  betaE_bin = 0.15,
  betaE_norm = -0.07,
  betaGE_normPRS_bin = 0,
  betaGE_normPRS_norm = 0,
  E_bin_freq = 0.745,
  strata = F
)
#Disease prevalence:  0.00813


#Start simulation

sim_N = 1000

res = list()
res_spmle = list()
res_caseonly = list()


for (i in 1:sim_N) {
  set.seed(i)
  dat = simFit(
    ncontrol = 500,
    ncase = 500,
    beta0 = -4.8,
    betaG_normPRS = 0.450,
    betaE_bin = 0.15,
    betaE_norm = -0.07,
    betaGE_normPRS_bin = 0,
    betaGE_normPRS_norm = 0,
    E_bin_freq = 0.745,
    strata = F
  )
  dat0 = cbind(dat$D, dat$G, dat$E)
  colnames(dat0) = c("D", "prs", "envir1", "envir2")
  dat0 = data.frame(dat0)
  
  res[[i]] = prs_e_function_gr(
    data = dat0,
    formula = D ~ prs + envir1 + envir2 + envir1:prs + envir2:prs,
    formula_prs = prs ~ 1,
    facVar = NULL,
    numDeriv = F,
    side0 = 2
  )
  
  #case only
  res_caseonly[[i]] = function_caseonly(data_sim = dat, strata = F)
  
  #nonparametric
  spmleFull = spmle(
    D = dat$D,
    G = dat$G,
    E = dat$E,
    pi1 = 0
  )
  res_spmle[[i]] = summary(spmleFull)$coefficients
}

save(res, res_caseonly, res_spmle, file = paste0(path, "sim_s1.rda"))






#-------------------------------------------------
#Scenario 2: with strata c(0.1,0.3,0.5), int1=0, int2=-0.4, no interaction between E and S
#Rerun with a lower rare disease assumption
#simulation work for stratified model
#Jan 6, 2023
#-------------------------------------------------

library(caseControlGE)
path <- "/users/zwang4/GXE/simulation_results_manuscript_edit/rerun/"


#look at the population disease prevalence
dat = simFit(
  ncontrol = 500,
  ncase = 500,
  beta0 = -7,
  betaG_normPRS = 0.450,
  betaE_bin = 0.15,
  betaE_norm = -0.07,
  betaGE_normPRS_bin = 0,
  betaGE_normPRS_norm = -0.4,
  E_bin_freq = 0.745,
  strata = T,
  interact_strata_prs = F,
  beta_strata = c(0.1, 0.3, 0.5),
  mu_strata = c(1, 0.2, 0.3, 0.2)
)
#Disease prevalence: 0.006494





#Start simulation

sim_N = 1000

res = list()
res_spmle = list()
res_caseonly = list()


for (i in 1:sim_N) {
  set.seed(i)
  dat = simFit(
    ncontrol = 500,
    ncase = 500,
    beta0 = -7,
    betaG_normPRS = 0.450,
    betaE_bin = 0.15,
    betaE_norm = -0.07,
    betaGE_normPRS_bin = 0,
    betaGE_normPRS_norm = -0.4,
    E_bin_freq = 0.745,
    strata = T,
    interact_strata_prs = F,
    beta_strata = c(0.1, 0.3, 0.5),
    mu_strata = c(1, 0.2, 0.3, 0.2)
  )
  
  dat0 = cbind(dat$D, dat$G, dat$E, dat$S)
  
  colnames(dat0) = c("D", "prs", "envir1", "envir2", "s1", "s2")
  dat0 = data.frame(dat0)
  dat0$s1 = factor(dat0$s1)
  dat0$envir1 = factor(dat0$envir1)
  
  res[[i]] = prs_e_function_gr(
    data = dat0,
    formula = D ~ prs + envir1 + envir2 + factor(s1) + s2 + envir1:prs + envir2:prs,
    formula_prs = prs ~ factor(s1) + s2,
    facVar = c("s1"),
    numDeriv = F,
    initial_empirical = F,
    initial_eta_sigma = c(1, 0.4, 0.8, 0.9, 0.25, 0.5, 1),
    side0 = 2
  )
  
  
  #case only
  res_caseonly[[i]] = function_caseonly(data_sim = dat, strata = T)
  
  #nonparametric
  s = dat$S
  race = dat$S[, 1]
  #Make sure only include two dummy variables for a factor of 3
  s <- cbind(s, ifelse(race == '2', 1, 0))
  s <- cbind(s, ifelse(race == '3', 1, 0))
  s = as.matrix(s)
  s <- apply(s, 2, as.numeric)
  
  spmleFull = spmle(
    D = dat$D,
    G = dat$G,
    E = cbind(dat$E, s[, c(3, 4, 2)]),
    pi1 = 0
  )
  res_spmle[[i]] = summary(spmleFull)$coefficients
  
  
}

save(res, res_caseonly, res_spmle, file = paste0(path, "sim_s2_rerun.rda"))





#-------------------------------------------------
#Scenario 3: with strata c(0.02,0.2,0.5), int1=0, int2=-0.4, interaction between prs and S 0,0.1,0.2
#simulation work for stratified model
#Jan 6, 2023
#-------------------------------------------------

library(caseControlGE)
path <- "/users/zwang4/GXE/simulation_results_manuscript_edit/"



#look at the population disease prevalence
dat = simFit(
  ncontrol = 500,
  ncase = 500,
  beta0 = -6,
  betaG_normPRS = 0.45,
  betaE_bin = 0.15,
  betaE_norm = -0.07,
  betaGE_normPRS_bin = 0,
  betaGE_normPRS_norm = -0.4,
  E_bin_freq = 0.745,
  strata = T,
  interact_strata_prs = T,
  beta_strata = c(0.1, 0.3, 0.5),
  beta_strata_interact = c(0, 0.1, 0.2),
  mu_strata = c(1, 0.2, 0.3, 0.2)
)
#Disease prevalence:  0.011257




#Start simulation

sim_N = 1000

res = list()
res_spmle = list()
res_caseonly = list()


for (i in 1:sim_N) {
  dat = simFit(
    ncontrol = 500,
    ncase = 500,
    beta0 = -6,
    betaG_normPRS = 0.45,
    betaE_bin = 0.15,
    betaE_norm = -0.07,
    betaGE_normPRS_bin = 0,
    betaGE_normPRS_norm = -0.4,
    E_bin_freq = 0.745,
    strata = T,
    interact_strata_prs = T,
    beta_strata = c(0.1, 0.3, 0.5),
    beta_strata_interact = c(0, 0.1, 0.2),
    mu_strata = c(1, 0.2, 0.3, 0.2)
  )
  
  dat0 = cbind(dat$D, dat$G, dat$E, dat$S)
  
  colnames(dat0) = c("D", "prs", "envir1", "envir2", "s1", "s2")
  dat0 = data.frame(dat0)
  dat0$s1 = factor(dat0$s1)
  dat0$envir1 = factor(dat0$envir1)
  
  res[[i]] = prs_e_function_gr(
    data = dat0,
    formula = D ~ prs + envir1 + envir2 + factor(s1) + s2 + envir1:prs + envir2:prs +
      factor(s1):prs + s2:prs,
    formula_prs = prs ~ factor(s1) + s2,
    facVar = c("s1"),
    initial_empirical = F,
    numDeriv = F,
    initial_eta_sigma = c(1, 0.4, 0.8, 0.9, 0.25, 0.5, 1),
    side0 = 2
  )
  
  #case only
  res_caseonly[[i]] = function_caseonly(data_sim = dat, strata = T)
  
  #nonparametric
  s = dat$S
  race = dat$S[, 1]
  #Make sure only include two dummy variables for a factor of 3
  s <- cbind(s, ifelse(race == '2', 1, 0))
  s <- cbind(s, ifelse(race == '3', 1, 0))
  s = as.matrix(s)
  s <- apply(s, 2, as.numeric)
  
  spmleFull = spmle(
    D = dat$D,
    G = dat$G,
    E = cbind(dat$E, s[, c(3, 4, 2)]),
    pi1 = 0
  )
  res_spmle[[i]] = summary(spmleFull)$coefficients
  
}

save(res, res_caseonly, res_spmle, file = paste0(path, "sim_s3_rerun.rda"))





#-------------------------------------------------
#Scenario 4:  int1=0, int2=-0.4, only strata for prs, no strata in disease model, correlated E2 and S2
#simulation work for stratified model
#Jan 6, 2023
#-------------------------------------------------

library(caseControlGE)
path <- "/users/zwang4/GXE/simulation_results_manuscript_edit/rerun/"

#Modify the simulation function a bit to incorporate correlations between E and S
simFit <- function(ncontrol = 500,
                   ncase = 500,
                   beta0 = -5,
                   betaG_normPRS = 0.450,
                   betaE_bin = 0.143,
                   betaE_norm = -0.019,
                   betaGE_normPRS_bin = 0,
                   betaGE_normPRS_norm = -0.40,
                   E_bin_freq = 0.745,
                   strata = TRUE,
                   interact_strata_prs = FALSE,
                   beta_strata = c(0.1, 0.1, 0.1),
                   beta_strata_interact = c(0.4, 0.4, 0.4),
                   sigma_strata = c(0.25, 0.5, 1),
                   mu_strata = c(1, 0.2, 0.3, 0.2)) {
  n_pop = 1000000
  library(MASS)
  if (strata == TRUE) {
    my_mu1 <-
      c(0, 0)                                   # Specify the means of the variables
    my_Sigma1 <-
      matrix(c(1, 0.4, 0.4, 1),                 # Specify the covariance matrix of the variables
             ncol = 2)
    
    test = mvrnorm(n = n_pop, mu = my_mu1, Sigma = my_Sigma1)
    E_sim_norm <- test[, 1]
    strata_norm = test[, 2]
    
    race = sample(c(1, 2, 3), size = n_pop, replace = T)
    sigma_sim = race
    sigma_sim[race == 1] = sigma_strata[1]
    sigma_sim[race == 2] = sigma_strata[2]
    sigma_sim[race == 3] = sigma_strata[3]
    
    
    race_dummy = array(0, c(n_pop, 2))
    race_dummy[, 1] <- ifelse(race == '2', 1, 0)
    race_dummy[, 2] <- ifelse(race == '3', 1, 0)
    S_sim.tmp = cbind(rep(1, n_pop), race_dummy, strata_norm)
    tmp_mu_sim = as.numeric(as.vector(S_sim.tmp %*% mu_strata))
    
    prs <- rnorm(n = n_pop, mean = tmp_mu_sim, sd = sigma_sim)
    E_sim_bin <- rbinom(n = n_pop, size = 1, prob = E_bin_freq)
    #  E_sim_norm <- rnorm(n=n_pop)
  } else {
    prs <- rnorm(n = n_pop)
    E_sim_bin <- rbinom(n = n_pop, size = 1, prob = E_bin_freq)
    E_sim_norm <- rnorm(n = n_pop)
    
  }
  
  
  
  prs <- as.numeric(as.vector(scale(
    prs, center = TRUE, scale = F
  )))
  
  if (strata == TRUE) {
    if (interact_strata_prs == TRUE) {
      xb <-
        beta0 + betaG_normPRS * prs + betaE_bin * E_sim_bin + betaE_norm * E_sim_norm + betaGE_normPRS_bin *
        E_sim_bin * prs + betaGE_normPRS_norm * E_sim_norm * prs + S_sim.tmp[, -1] %*% beta_strata + beta_strata_interact[1] *
        prs * S_sim.tmp[, 2] + beta_strata_interact[2] * prs * S_sim.tmp[, 3] +
        beta_strata_interact[3] * prs * S_sim.tmp[, 4]
    } else {
      xb <-
        beta0 + betaG_normPRS * prs + betaE_bin * E_sim_bin + betaE_norm * E_sim_norm + betaGE_normPRS_bin *
        E_sim_bin * prs + betaGE_normPRS_norm * E_sim_norm * prs
    }
  } else{
    xb <-
      beta0 + betaG_normPRS * prs + betaE_bin * E_sim_bin + betaE_norm * E_sim_norm + betaGE_normPRS_bin *
      E_sim_bin * prs + betaGE_normPRS_norm * E_sim_norm * prs
  }
  prob <- 1 / (1 + exp(-xb))
  
  D_sim <-
    rbinom(n = n_pop, size = 1, prob = prob)  # 1=case, 0=control
  cat("Disease prevalence: ", mean(D_sim))
  
  
  
  E_sim = cbind(E_sim_bin, E_sim_norm)
  
  mean_prs = mean(prs)
  sd_prs = sd(prs)
  
  # randomly select the ncontrol controls and ncase cases
  id1 = sample(which(D_sim == 0), size = ncontrol, replace = F)
  id2 = sample(which(D_sim == 1), size = ncase, replace = F)
  D = c(D_sim[id1], D_sim[id2])
  G = c(prs[id1], prs[id2])
  E = rbind(as.matrix(E_sim[id1, ]), as.matrix(E_sim[id2, ]))
  if (strata == TRUE) {
    S_sim = cbind(race, strata_norm)
    S = rbind(as.matrix(S_sim[id1, ]), as.matrix(S_sim[id2, ]))
  } else{
    S = NULL
  }
  
  #### Bind D, G, and E into a list
  dat = list(
    D = D,
    G = G,
    E = E,
    S = S,
    prevelance = mean(D_sim),
    mean_prs = mean_prs,
    sd_prs = sd_prs
  )
  return(dat)
  
}




#look at the population disease prevalence
dat = simFit(
  ncontrol = 500,
  ncase = 500,
  beta0 = -5,
  betaG_normPRS = 0.45,
  betaE_bin = 0.143,
  betaE_norm = -0.019,
  betaGE_normPRS_bin = 0,
  betaGE_normPRS_norm = -0.4,
  E_bin_freq = 0.745,
  strata = T,
  interact_strata_prs = F,
  beta_strata = c(0.1, 0.3, 0.5),
  beta_strata_interact = c(0, 0.1, 0.2)
)
#Disease prevalence:  0.00999





#Start simulation

sim_N = 1000

res = list()
res_spmle = list()
res_caseonly = list()


for (i in 1:sim_N) {
  set.seed(i)
  dat = simFit(
    ncontrol = 500,
    ncase = 500,
    beta0 = -5,
    betaG_normPRS = 0.45,
    betaE_bin = 0.15,
    betaE_norm = -0.07,
    betaGE_normPRS_bin = 0,
    betaGE_normPRS_norm = -0.4,
    E_bin_freq = 0.745,
    strata = T,
    interact_strata_prs = F,
    beta_strata = c(0.1, 0.3, 0.5),
    beta_strata_interact = c(0, 0.1, 0.2)
  )
  dat0 = cbind(dat$D, dat$G, dat$E, dat$S)
  colnames(dat0) = c("D", "prs", "envir1", "envir2", "s1", "s2")
  dat0 = data.frame(dat0)
  dat0$s1 = factor(dat0$s1)
  dat0$envir1 = factor(dat0$envir1)
  
  res[[i]] = prs_e_function_gr(
    data = dat0,
    formula = D ~ prs + envir1 + envir2 + envir1:prs + envir2:prs,
    formula_prs = prs ~ factor(s1) + s2,
    facVar = c("s1"),
    initial_empirical = F,
    numDeriv = F,
    initial_eta_sigma = c(1, 0.4, 0.8, 0.9, 0.25, 0.5, 1),
    side0 = 2
  )
  
  #case only
  res_caseonly[[i]] = function_caseonly(data_sim = dat, strata = F)
  
  #nonparametric
  spmleFull = spmle(
    D = dat$D,
    G = dat$G,
    E = dat$E,
    pi1 = dat$prevelance
  )
  res_spmle[[i]] = summary(spmleFull)$coefficients
  
}

save(res, res_caseonly, res_spmle, file = paste0(path, "sim_s4.rda"))






#-------------------------------------------------
#Scenario 5:  int1=0, int2=-0.4, only strata for prs, add strata in disease model but no interaction terms, correlated E2 and S2
#simulation work for stratified model
#Jan 6, 2023
#-------------------------------------------------

library(caseControlGE)
path <- "/users/zwang4/GXE/simulation_results_manuscript_edit/rerun/"


#Modify the simulation function a bit to incorporate correlations between E and S
simFit <- function(ncontrol = 500,
                   ncase = 500,
                   beta0 = -5,
                   betaG_normPRS = 0.450,
                   betaE_bin = 0.143,
                   betaE_norm = -0.019,
                   betaGE_normPRS_bin = 0,
                   betaGE_normPRS_norm = -0.40,
                   E_bin_freq = 0.745,
                   strata = TRUE,
                   interact_strata_prs = FALSE,
                   beta_strata = c(0.1, 0.1, 0.1),
                   beta_strata_interact = c(0.4, 0.4, 0.4),
                   sigma_strata = c(0.25, 0.5, 1),
                   mu_strata = c(1, 0.2, 0.3, 0.2)) {
  n_pop = 1000000
  library(MASS)
  if (strata == TRUE) {
    my_mu1 <-
      c(0, 0)                                   # Specify the means of the variables
    my_Sigma1 <-
      matrix(c(1, 0.4, 0.4, 1),                 # Specify the covariance matrix of the variables
             ncol = 2)
    
    test = mvrnorm(n = n_pop, mu = my_mu1, Sigma = my_Sigma1)
    E_sim_norm <- test[, 1]
    strata_norm = test[, 2]
    
    race = sample(c(1, 2, 3), size = n_pop, replace = T)
    sigma_sim = race
    sigma_sim[race == 1] = sigma_strata[1]
    sigma_sim[race == 2] = sigma_strata[2]
    sigma_sim[race == 3] = sigma_strata[3]
    
    
    race_dummy = array(0, c(n_pop, 2))
    race_dummy[, 1] <- ifelse(race == '2', 1, 0)
    race_dummy[, 2] <- ifelse(race == '3', 1, 0)
    S_sim.tmp = cbind(rep(1, n_pop), race_dummy, strata_norm)
    tmp_mu_sim = as.numeric(as.vector(S_sim.tmp %*% mu_strata))
    
    prs <- rnorm(n = n_pop, mean = tmp_mu_sim, sd = sigma_sim)
    E_sim_bin <- rbinom(n = n_pop, size = 1, prob = E_bin_freq)
    #  E_sim_norm <- rnorm(n=n_pop)
  } else {
    prs <- rnorm(n = n_pop)
    E_sim_bin <- rbinom(n = n_pop, size = 1, prob = E_bin_freq)
    E_sim_norm <- rnorm(n = n_pop)
    
  }
  
  
  
  prs <- as.numeric(as.vector(scale(
    prs, center = TRUE, scale = F
  )))
  
  if (strata == TRUE) {
    if (interact_strata_prs == TRUE) {
      xb <-
        beta0 + betaG_normPRS * prs + betaE_bin * E_sim_bin + betaE_norm * E_sim_norm + betaGE_normPRS_bin *
        E_sim_bin * prs + betaGE_normPRS_norm * E_sim_norm * prs + S_sim.tmp[, -1] %*% beta_strata + beta_strata_interact[1] *
        prs * S_sim.tmp[, 2] + beta_strata_interact[2] * prs * S_sim.tmp[, 3] +
        beta_strata_interact[3] * prs * S_sim.tmp[, 4]
    } else {
      xb <-
        beta0 + betaG_normPRS * prs + betaE_bin * E_sim_bin + betaE_norm * E_sim_norm + betaGE_normPRS_bin *
        E_sim_bin * prs + betaGE_normPRS_norm * E_sim_norm * prs + S_sim.tmp[, -1] %*% beta_strata
    }
  } else{
    xb <-
      beta0 + betaG_normPRS * prs + betaE_bin * E_sim_bin + betaE_norm * E_sim_norm + betaGE_normPRS_bin *
      E_sim_bin * prs + betaGE_normPRS_norm * E_sim_norm * prs
  }
  prob <- 1 / (1 + exp(-xb))
  
  D_sim <-
    rbinom(n = n_pop, size = 1, prob = prob)  # 1=case, 0=control
  cat("Disease prevalence: ", mean(D_sim))
  
  
  
  E_sim = cbind(E_sim_bin, E_sim_norm)
  
  mean_prs = mean(prs)
  sd_prs = sd(prs)
  
  # randomly select the ncontrol controls and ncase cases
  id1 = sample(which(D_sim == 0), size = ncontrol, replace = F)
  id2 = sample(which(D_sim == 1), size = ncase, replace = F)
  D = c(D_sim[id1], D_sim[id2])
  G = c(prs[id1], prs[id2])
  E = rbind(as.matrix(E_sim[id1, ]), as.matrix(E_sim[id2, ]))
  if (strata == TRUE) {
    S_sim = cbind(race, strata_norm)
    S = rbind(as.matrix(S_sim[id1, ]), as.matrix(S_sim[id2, ]))
  } else{
    S = NULL
  }
  
  #### Bind D, G, and E into a list
  dat = list(
    D = D,
    G = G,
    E = E,
    S = S,
    prevelance = mean(D_sim),
    mean_prs = mean_prs,
    sd_prs = sd_prs
  )
  return(dat)
  
}




#look at the population disease prevalence
dat = simFit(
  ncontrol = 500,
  ncase = 500,
  beta0 = -5,
  betaG_normPRS = 0.45,
  betaE_bin = 0.143,
  betaE_norm = -0.019,
  betaGE_normPRS_bin = 0,
  betaGE_normPRS_norm = -0.4,
  E_bin_freq = 0.745,
  strata = T,
  interact_strata_prs = F,
  beta_strata = c(0, 0.3, 0.5),
  beta_strata_interact = c(0, 0.1, 0.2)
)
#Disease prevalence:  0.00999





#Start simulation

sim_N = 1000

res = list()
res_spmle = list()
res_caseonly = list()


for (i in 1:sim_N) {
  set.seed(i)
  dat = simFit(
    ncontrol = 500,
    ncase = 500,
    beta0 = -5,
    betaG_normPRS = 0.45,
    betaE_bin = 0.15,
    betaE_norm = -0.07,
    betaGE_normPRS_bin = 0,
    betaGE_normPRS_norm = -0.4,
    E_bin_freq = 0.745,
    strata = T,
    interact_strata_prs = F,
    beta_strata = c(0.1, 0.3, 0.5),
    beta_strata_interact = c(0, 0.1, 0.2)
  )
  dat0 = cbind(dat$D, dat$G, dat$E, dat$S)
  colnames(dat0) = c("D", "prs", "envir1", "envir2", "s1", "s2")
  dat0 = data.frame(dat0)
  dat0$s1 = factor(dat0$s1)
  dat0$envir1 = factor(dat0$envir1)
  
  res[[i]] = prs_e_function_gr(
    data = dat0,
    formula = D ~ prs + envir1 + envir2 + factor(s1) + s2 + envir1:prs + envir2:prs,
    formula_prs = prs ~ factor(s1) + s2,
    facVar = c("s1"),
    initial_empirical = F,
    numDeriv = F,
    initial_eta_sigma = c(1, 0.4, 0.8, 0.9, 0.25, 0.5, 1),
    side0 = 2
  )
  
  #case only
  res_caseonly[[i]] = function_caseonly(data_sim = dat, strata = T)
  
  #nonparametric
  s = dat$S
  race = dat$S[, 1]
  #Make sure only include two dummy variables for a factor of 3
  s <- cbind(s, ifelse(race == '2', 1, 0))
  s <- cbind(s, ifelse(race == '3', 1, 0))
  s = as.matrix(s)
  s <- apply(s, 2, as.numeric)
  
  spmleFull = spmle(
    D = dat$D,
    G = dat$G,
    E = cbind(dat$E, s[, c(3, 4, 2)]),
    pi1 = 0
  )
  res_spmle[[i]] = summary(spmleFull)$coefficients
  
  
}

save(res, res_caseonly, res_spmle, file = paste0(path, "sim_s5.rda"))





#-------------------------------------------------
#Scenario 6: with strata c(0.02,0.2,0.5), int1=0, int2=-0.4, interaction between prs and S 0,0.1,0.2
#simulation work for stratified model
#Jan 6, 2023
#-------------------------------------------------



library(caseControlGE)
path <- "/users/zwang4/GXE/simulation_results_manuscript_edit/rerun/"
source("/users/zwang4/GXE/sim_data_function.R")



#look at the population disease prevalence
dat = simFit(
  ncontrol = 500,
  ncase = 500,
  beta0 = -6,
  betaG_normPRS = 0.45,
  betaE_bin = 0,
  betaE_norm = 0,
  betaGE_normPRS_bin = 0,
  betaGE_normPRS_norm = 0,
  E_bin_freq = 0.745,
  strata = T,
  interact_strata_prs = T,
  beta_strata = c(0.1, 0.3, 0.5),
  beta_strata_interact = c(0, 0.1, 0.2),
  mu_strata = c(1, 0.2, 0.3, 0.2)
)
#Disease prevalence:  0.011257


#Modify the case-only method function for scenario 6
function_caseonly1 <- function(data_sim, strata = FALSE) {
  D = data_sim$D
  prs = data_sim$G
  
  mean_prs = data_sim$mean_prs
  sd_prs = data_sim$sd_prs
  
  fit_caseonly <-
    lm(prs[D == 1] ~ factor(data_sim$S[D == 1, 1]) + data_sim$S[D == 1, 2])
  
  
  
  beta = fit_caseonly$coefficients
  beta_int = fit_caseonly$coefficients[-1] / (sd_prs) ^ 2
  sd_int = summary(fit_caseonly)$coef[-1, 2] / (sd_prs) ^ 2
  
  beta_prs = (fit_caseonly$coefficients[1] - mean_prs) / (sd_prs) ^ 2
  sd_prs1 = sqrt(((summary(fit_caseonly)$coef[1, 2]) ^ 2 + (sd_prs) ^ 2 /
                    1000000) / (sd_prs) ^ 4)
  
  
  res = summary.caseonly(parms = c(beta_prs, beta_int),
                         sd = c(sd_prs1, sd_int))
  
  
  beta_int2 = fit_caseonly$coefficients[-1] / (sd(fit_caseonly$residuals)) ^
    2
  sd_int2 = summary(fit_caseonly)$coef[-1, 2] / (sd(fit_caseonly$residuals)) ^
    2
  
  beta_prs2 = (fit_caseonly$coefficients[1] - mean_prs) / (sd(fit_caseonly$residuals)) ^
    2
  sd_prs2 = sqrt(((summary(fit_caseonly)$coef[1, 2]) ^ 2 + (sd(
    fit_caseonly$residuals
  )) ^ 2 / 1000000) / (sd(fit_caseonly$residuals)) ^ 4)
  res2 = summary.caseonly(parms = c(beta_prs2, beta_int2),
                          sd = c(sd_prs2, sd_int2))
  rownames(res2) = paste0(rownames(res2), "_resid_sd")
  res = rbind(res, res2)
  
  return(res)
}


#Start simulation

sim_N = 1000

res = list()
res_spmle = list()
res_caseonly = list()


for (i in 1:sim_N) {
  set.seed(i)
  dat = simFit(
    ncontrol = 500,
    ncase = 500,
    beta0 = -6,
    betaG_normPRS = 0.45,
    betaE_bin = 0,
    betaE_norm = 0,
    betaGE_normPRS_bin = 0,
    betaGE_normPRS_norm = 0,
    E_bin_freq = 0.745,
    strata = T,
    interact_strata_prs = T,
    beta_strata = c(0.1, 0.3, 0.5),
    beta_strata_interact = c(0, 0.1, 0.2),
    mu_strata = c(1, 0.2, 0.3, 0.2)
  )
  
  dat0 = cbind(dat$D, dat$G, dat$E, dat$S)
  
  colnames(dat0) = c("D", "prs", "envir1", "envir2", "s1", "s2")
  dat0 = data.frame(dat0)
  dat0$s1 = factor(dat0$s1)
  dat0$envir1 = factor(dat0$envir1)
  
  res[[i]] = prs_e_function_gr(
    data = dat0,
    formula = D ~ prs + factor(s1) + s2 + factor(s1):prs + s2:prs,
    formula_prs = prs ~ factor(s1) + s2,
    facVar = c("s1"),
    initial_empirical = F,
    numDeriv = F,
    initial_eta_sigma = c(1, 0.4, 0.8, 0.9, 0.25, 0.5, 1),
    side0 = 2
  )
  
  #case only
  res_caseonly[[i]] = function_caseonly1(data_sim = dat, strata = F)
  
  #nonparametric
  s = dat$S
  race = dat$S[, 1]
  #Make sure only include two dummy variables for a factor of 3
  s <- cbind(s, ifelse(race == '2', 1, 0))
  s <- cbind(s, ifelse(race == '3', 1, 0))
  s = as.matrix(s)
  s <- apply(s, 2, as.numeric)
  
  spmleFull = spmle(
    D = dat$D,
    G = dat$G,
    E = s[, c(3, 4, 2)],
    pi1 = 0
  )
  res_spmle[[i]] = summary(spmleFull)$coefficients
  
}

save(res, res_caseonly, res_spmle, file = paste0(path, "sim_s6.rda"))
