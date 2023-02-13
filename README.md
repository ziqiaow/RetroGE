# RetroGE
A Retrospective Likelihood Method for Joint Modeling of Gene-Environment Interactions and Correlations using Polygenic Risk Score (PRS) in Case-Control Studies

**Reference**

* Ziqiao Wang, Wen Shi, Raymond J. Carrol, and Nilanjan Chatterjee (2023). Joint Modeling of Gene-Environment Correlations and Interactions
using Polygenic Risk Scores in Case-Control Studies.

## Example Analysis
We provide a simple example for running the retrospective likelihood method using simulated data. Complete R codes for the simulation study are in [Simulations](simulations/simulation.R). The R function of the proposed method is in [RetroGE](R/RetroGE.R). The R function to simulate data is available here [sim_data_function](R/sim_data_function.R).
```
rm(list = ls())
source("./R/RetroGE.R")
source("./R/sim_data_function.R")
```
First simulate some data based on a population disease model with two independent environmental variables $E_1$ (binary) and $E_2$ (continuous) and two environmental variables $S_1$ (categorical variable with 3 levels) and $S_2$ (continuous) that are correlated to PRS (Z), i.e., $$logit(Pr(D=1|Z,E,S)) = \alpha+\beta_ZZ+\beta_{E_1}E_1+\beta_{E_2}E_2+\beta_{S_{1(1)}}S_{1(1)}+\beta_{S_{2(1)}}S_{2(1)}+\beta_{S_2}S_{2}+\beta_{E_1Z}E_1Z+\beta_{E_2Z}E_2Z$$
PRS value in the population follows
$$Z|S \sim N(\eta_0+\eta_{S_{1(1)}}S_{1(1)}+\eta_{S_{2(1)}}S_{2(1)}+\eta_{S_{2}}S_{2},\sigma^2_{S_1})$$
This is the underlying model we assumed in [Simulations](simulations/simulation.R) Scenario 2 as described in the original article.
```
set.seed(02122023)
dat = simFit(
  ncontrol = 500,
  ncase = 500,
  beta0 = -5,
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
#Disease prevalence:  0.011445
```
View the simulated case-control data matrix.
```
#Convert the list to a data frame
dat0 = cbind(dat$D, dat$G, dat$E, dat$S)
colnames(dat0) = c("D", "prs", "envir1", "envir2", "s1", "s2")
dat0 = data.frame(dat0)
dat0$s1 = factor(dat0$s1)
dat0$envir1 = factor(dat0$envir1)
dim(dat0)
#[1] 1000    6
table(dat0$D) #500 cases and 500 controls
#  0   1 
#500 500 
head(dat0)
```
Fit the retrospective likelihood method
```
res = prs_e_function_gr(
    data = dat0,
    formula = D ~ prs + envir1 + envir2 + factor(s1) + s2 + envir1:prs + envir2:prs,
    formula_prs = prs ~ factor(s1) + s2,
    facVar = c("s1"),
    numDeriv = FALSE, #Use analytical score function (FALSE) for MLE or numerical gradient values (TRUE), note that the two results are similar but analytical is faster
    initial_empirical = TRUE, #Use the fitted values from logistic regression (disease model) and linear regression on the control samples (PRS model) for initial value inputs
    initial_eta_sigma = NULL,
    side0 = 2
  )
#After removing missing values, the number of observations is 1000 
#initial  value 1716.201381 
#final  value 1422.216099 
#converged
  ```

## Codes and Results of UK Biobank Data Analysis
The complete R codes and results for the data analysis of UK Biobank is available in R markdown.
* [Incident breast cancer for postmenopausal and premenopausal women](https://raw.githack.com/ziqiaow/RetroGE/main/results/UKB_breastcancer.html)
* [Incident colorectal cancer](https://raw.githack.com/ziqiaow/RetroGE/main/results/report_colorectal.html)

## Questions
Please feel free to email any questions/suggestions at zwang389@jhu.edu
