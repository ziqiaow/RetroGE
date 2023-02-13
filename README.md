# RetroGE
A Retrospective Likelihood Method for Joint Modeling of Gene-Environment Interactions and Correlations using Polygenic Risk Score (PRS) in Case-Control Studies

**Reference**

* Ziqiao Wang, Wen Shi, Raymond J. Carrol, and Nilanjan Chatterjee (2023). Joint Modeling of Gene-Environment Correlations and Interactions
using Polygenic Risk Scores in Case-Control Studies.

### Example Analysis
We provide a simple example for running the retrospective likelihood method using simulated data. Complete R codes for the simulation study are in [Simulations](simulations/simulation.R). The R function of the proposed method is in [RetroGE](R/RetroGE.R). The R function to generate data is available here [sim_data_function](R/sim_data_function.R).
```
rm(list = ls())
source("./R/RetroGE.R")
source("./R/sim_data_function.R")
```

First simulate some data. Suppose we simulate the data based on a full disease model with two independent environmental variables $E_1$ and $E_2$ and two environmental variables $S_1$ and $S_2$ that are correlated to PRS ($Z$), i.e., $logit(Pr(D=1|Z,E,S)) = \beta_ZZ+\beta_{E_1}E_1+\beta_{E_2}E_2+\beta_{S_{1(1)}}S_{1(1)}+\beta_{S_{2(1)}}S_{2(1)}+\beta_{S_2}S_{2}+\beta_{E_1Z}E_1Z+\beta_{E_2Z}E_2Z+\beta_{S_{1(1)}Z}S_{1(1)}Z+\beta_{S_{2(1)}Z}S_{2(1)}Z+\beta_{S_2Z}S_{2}Z$
```
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
```

### Codes and Results of UK Biobank Data Analysis
The complete R codes and results for the data analysis of UK Biobank is available in R markdown.
* [Incident breast cancer for postmenopausal and premenopausal women](https://raw.githack.com/ziqiaow/RetroGE/main/results/UKB_breastcancer.html)
* [Incident colorectal cancer](https://raw.githack.com/ziqiaow/RetroGE/main/results/report_colorectal.html)

### Questions
Please feel free to email any questions/suggestions at zwang389@jhu.edu
