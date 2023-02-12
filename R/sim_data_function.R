#-------------------------------------------------
#Simulate case-control data from a prospective cohort of 10^6 individuals
#The underlying disease model can incorporate multivariate independent and PRS-dependent environmental variables and their interaction terms with PRS
#The PRS normal model allows correlated environmental variables and stratification factors in the mean function and heteroskedasticity in the variance
#The number of cases and controls, log odds ratio parameters, regression coefficients can all be user-specified
#Feb 12, 2023
#-------------------------------------------------

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
                   mu_strata = c(1, 0.4, 0.8, 0.9)) {
  n_pop = 1000000
  
  if (strata == TRUE) {
    race = sample(c(1, 2, 3), size = n_pop, replace = T)
    sigma_sim = race
    sigma_sim[race == 1] = sigma_strata[1]
    sigma_sim[race == 2] = sigma_strata[2]
    sigma_sim[race == 3] = sigma_strata[3]
    
    strata_norm = rnorm(n = n_pop)
    race_dummy = array(0, c(n_pop, 2))
    race_dummy[, 1] <- ifelse(race == '2', 1, 0)
    race_dummy[, 2] <- ifelse(race == '3', 1, 0)
    S_sim.tmp = cbind(rep(1, n_pop), race_dummy, strata_norm)
    tmp_mu_sim = as.numeric(as.vector(S_sim.tmp %*% mu_strata))
    
    prs <- rnorm(n = n_pop, mean = tmp_mu_sim, sd = sigma_sim)
  } else {
    prs <- rnorm(n = n_pop)
  }
  
  E_sim_bin <- rbinom(n = n_pop, size = 1, prob = E_bin_freq)
  E_sim_norm <- rnorm(n = n_pop)
  
  
  
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
