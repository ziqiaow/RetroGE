# RetroGE
A Retrospective Likelihood Method for Joint Modeling of Gene-Environment Interactions and Correlations using Polygenic Risk Score (PRS) in Case-Control Studies

**Reference**

* Ziqiao Wang, Wen Shi, Raymond J Carroll, Nilanjan Chatterjee, Joint Modeling of Gene-Environment Correlations and Interactions Using Polygenic Risk Scores in Case-Control Studies, American Journal of Epidemiology, 2024; https://doi.org/10.1093/aje/kwae081

## Example Analysis
We provide a simple example for running the retrospective likelihood method using simulated data. Complete R codes for the simulation study are in [simulations](simulations). The R function of the proposed method is in [RetroGE](R/RetroGE.R). The R function to simulate data is available here [sim_data_function](R/sim_data_function.R).
```
rm(list = ls())
library(devtools)
source_url("https://github.com/ziqiaow/RetroGE/blob/main/R/RetroGE.R?raw=TRUE")
source_url("https://github.com/ziqiaow/RetroGE/blob/main/R/sim_data_function.R?raw=TRUE")

#If directly downloaded the R files from Github to your local directory, you can load the R functions using source().
source("./R/RetroGE.R")
source("./R/sim_data_function.R")
```
First simulate a full cohort of 1000000 individuals based on a population disease model with two independent environmental variables $E_1$ (binary) and $E_2$ (continuous) and two environmental variables $S_1$ (categorical variable with 3 levels) and $S_2$ (continuous) that are correlated to PRS (Z), i.e., $$logit(Pr(D=1|Z,E,S)) = \alpha+\beta_ZZ+\beta_{E_1}E_1+\beta_{E_2}E_2+\beta_{S_{1(1)}}S_{1(1)}+\beta_{S_{2(1)}}S_{2(1)}+\beta_{S_2}S_{2}+\beta_{E_1Z}E_1Z+\beta_{E_2Z}E_2Z$$
PRS value in the population follows
$$Z|S \sim N(\eta_0+\eta_{S_{1(1)}}S_{1(1)}+\eta_{S_{2(1)}}S_{2(1)}+\eta_{S_{2}}S_{2},\sigma^2_{S_1})$$
This is the underlying model we assumed in [Scenario 2](simulation/simulation.R) as described in the original article.
```
set.seed(02122023)
dat = simFit(
  ncontrol = 500, #Simulate a cohort of 1000000 individuals, randomly select 500 controls and 500 cases
  ncase = 500,
  beta0 = -5,
  betaG_normPRS = 0.450,
  betaE_bin = 0.15,
  betaE_norm = -0.07,
  betaGE_normPRS_bin = 0,
  betaGE_normPRS_norm = -0.4,
  E_bin_freq = 0.745,
  strata = TRUE,
  interact_strata_prs = FALSE,
  beta_strata = c(0.1, 0.3, 0.5),
  mu_strata = c(1, 0.2, 0.3, 0.2),
  sigma_strata = c(0.25, 0.5, 1)
)
#Disease prevalence:  0.011445
```
View the simulated case-control data matrix, we randomly selected 500 cases and 500 controls.
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
dat0[1:5,]
#  D          prs envir1     envir2 s1         s2
#1 0  0.269289140      0  0.2074758  2  1.3675722
#2 0 -0.012508509      1 -0.5201500  3 -0.9360032
#3 0  1.025145443      1  1.1716256  3  0.2409941
#4 0 -0.264052650      1 -1.0113151  1  0.1802841
#5 0  0.008607869      1 -0.3353061  1  1.4849816
```
Fit the retrospective likelihood method to the case-control dataset.
```
startTime <- Sys.time()
res = prs_e_function_gr(
    data = dat0,
    formula = D ~ prs + envir1 + envir2 + factor(s1) + s2 + envir1:prs + envir2:prs, #Disease model
    formula_prs = prs ~ factor(s1) + s2, #PRS mean model
    facVar = c("s1"), #PRS variance
    numDeriv = FALSE, #Use analytical score function for MLE (more accurate and faster)
    initial_empirical = TRUE, #Use the fitted values from logistic regression (disease model) and linear regression on the control samples (PRS model) for initial value inputs
    initial_eta_sigma = NULL,
    side0 = 2
  )
#After removing missing values, the number of observations is 1000 
#initial  value 1716.201381 
#final  value 1422.216099 
#converged
endTime <- Sys.time()
```
Print running time.
```
print(endTime - startTime)
#Time difference of 0.09115911 secs
```
Output the final results. The output saves the fitted model formulas for both the disease model and PRS model, and the original data in 'model.info'. The retrospective likelihood method output is 'res_normal'. The standard logistic regression is also saved in the output as 'res_glm'. 
```
attributes(res)
#$names
#[1] "res_glm"    "res_normal" "model.info"
attributes(res$model.info)
#$names
#[1] "data"        "formula"     "formula_prs" "facVar"     

#Print the standard logistic regression results
res$res_glm
#                 Estimate Std. Error      z value     Pr(>|z|)
#(Intercept) -0.4191341519 0.18526102 -2.262397976 2.367283e-02
#prs          0.7714796054 0.23948822  3.221367633 1.275804e-03
#envir11      0.2536380216 0.17414419  1.456482807 1.452592e-01
#envir2       0.0003224133 0.06909941  0.004665934 9.962771e-01
#factor(s1)2 -0.1692431820 0.17499577 -0.967127259 3.334804e-01
#factor(s1)3 -0.0088460375 0.17725614 -0.049905393 9.601978e-01
#s2           0.5010261902 0.07239894  6.920352266 4.505219e-12
#prs:envir11 -0.3657490663 0.25574916 -1.430108562 1.526859e-01
#prs:envir2  -0.4391083120 0.09446156 -4.648539584 3.342935e-06

#Print the retrospective likelihood method results
res$res_normal
#                                Estimate  Std.Error      Z.value        Pvalue
#(Intercept)                 -0.400981115 0.18232555  -2.19925904  2.785951e-02
#prs                          0.626936447 0.14693917   4.26663937  1.984395e-05
#envir11                      0.222213575 0.16759097   1.32592810  1.848635e-01
#envir2                      -0.027940670 0.06714806  -0.41610541  6.773329e-01
#factor(s1)2                 -0.133973149 0.17401381  -0.76989954  4.413595e-01
#factor(s1)3                  0.004996946 0.17389114   0.02873606  9.770751e-01
#s2                           0.512355750 0.07118163   7.19786516  6.116258e-13
#prs:envir11                 -0.216526169 0.13755702  -1.57408307  1.154682e-01
#prs:envir2                  -0.318925205 0.05699702  -5.59547137  2.200232e-08
#eta_X.Intercept.            -0.150997513 0.01490880 -10.12807663  4.147267e-24
#eta_factor.s1.2              0.218147832 0.03219422   6.77599418  1.235537e-11
#eta_factor.s1.3              0.292126567 0.06963842   4.19490501  2.729860e-05
#eta_s2                       0.193974076 0.01164614  16.65565978  2.752914e-62
#sigma_stratadata[, facVar]1  0.242841841 0.01018782  23.83647848 1.398671e-125
#sigma_stratadata[, facVar]2  0.491772020 0.01888508  26.04023814 1.735463e-149
#sigma_stratadata[, facVar]3  0.978341123 0.03386884  28.88617162 1.781164e-183
```
In the printed retrospective likelihood method results, the 'Estimate' is the estimated log odds ratio for the disease model; for the PRS model (eta and sigma), the 'Estimate' is the estimated regression coefficient $\eta_S$ and estimated standard deviation $\sigma_{S_1}$.

## Codes and Results of UK Biobank Data Analysis
The complete R codes and results for the data analysis of UK Biobank is available in R markdown.
* [Incident breast cancer for postmenopausal and premenopausal women](https://raw.githack.com/ziqiaow/RetroGE/main/results/UKB_breastcancer.html)
* [Incident colorectal cancer](https://raw.githack.com/ziqiaow/RetroGE/main/results/report_colorectal.html)

Note that for real data applications, it is suggested to **standardize the continuous variables** that have largely varied scales for a stable numerical derivation. Simply use this R code to standardize the continuous variable x: 
```
scale(x, center = TRUE, scale = TRUE)
```

## Questions
Please send your questions/suggestions to zwang389@jhu.edu
