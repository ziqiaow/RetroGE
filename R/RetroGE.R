#PRS x E functions
#1. The proposed retrospective likelihood method
#2. Case-only method designed for the simulation study


#1. The proposed retrospective likelihood method
#If there is no stratification factor s: input: formula_prs=prs~1, s.var=NULL
#facVar is the factor variable for sigma, can be NULL
RetroGE <- function(data=dat0,
                              formula=D ~prs+envir1 +envir2+factor(s1)+s2+envir1:prs+envir2:prs+factor(s1):prs+s2:prs, #Make sure that Disease is coded as 1 and control is coded as 0.
                              formula_prs=prs~factor(s1)+s2,
                              facVar=c("s1"), #The factor variable for sigma in constructing PRS, can be any factor variable with different levels as user defined. Only accept one variable input or NULL.
                              initial_empirical=T, #Use logistic regression to assign initial value for optim() for MLE and linear regression on the control samples for the stratification variables S.
                              initial_eta_sigma=c(1,0.4,0.8,0.9,0.25,0.5,1), #The initial value for the stratification variables S. If initial_empirical = T, then this condition will be ignored.
                              numDeriv=F, #Whether to use analytical score function for MLE or numerical gradient values, note that the two results are similar but analytical is faster
                              side0=2){ #this is the wald test side, default is 2-sided test
  data=data.frame(data)
  options(na.action="na.pass")
  data_input=model.matrix(formula,data=data)
  s.var=model.matrix(formula_prs,data=data)
  tmp=cbind(data_input,s.var)
  id=which(complete.cases(tmp))
  data_input=data_input[id,]
  data=data[id,]
  s.var=data.frame(s.var[id,])
  cat(paste("After removing missing values, the number of observations is",dim(data_input)[1],"\n"))
  data_input2=cbind(data_input,data)
  data_input2 <- data_input2[, !duplicated(colnames(data_input2))]
  
 
  
  
  prs.var=as.character(formula_prs[[2]])
  prs=as.vector(data[,which(colnames(data)==prs.var)])
  disease.name=as.character(formula[[2]]) #all.vars(formula)[1]

  fit <- glm(formula, data=data,family = binomial(), model = FALSE, x = FALSE, y = FALSE)
  res_glm=summary(fit)$coef
  
  
  if(is.null(facVar)==TRUE){
    n_strat=dim(s.var)[2]
    n_sigma=1
    X.strata=s.var
    sigma.strata=rep(1,dim(X.strata)[1])
    dim(sigma.strata) <- c(length(sigma.strata), 1)
    
  } else {
    n_strat=dim(s.var)[2]
    X.strata=s.var
    sigma.strata=model.matrix( ~ data[,facVar]-1)
    n_sigma=dim(sigma.strata)[2]
  } 
  
  
  
  if(initial_empirical==TRUE){
    prs_control=prs[which(data[,which(colnames(data)==disease.name)]==0)]
    strata_control=data.frame(X.strata[which(data[,which(colnames(data)==disease.name)]==0),])
    data_control=cbind(prs_control,strata_control)
    ff=as.formula(paste0("prs_control","~ ", paste(c(1,colnames(strata_control)[-1]), collapse=" + ")))
    
    fit2=lm(ff,data=data.frame(data_control))
    param0=c(fit$coefficients,fit2$coefficients,rep(sd(fit2$residuals),n_sigma))
    
  } else {
    param0=c(fit$coefficients,initial_eta_sigma)
  }
  
  names(param0)=c(names(fit$coefficients),paste0("eta_",colnames(X.strata)),paste0("sigma_strata",colnames(sigma.strata)))
  nbeta <- length(names(fit$coefficients))-1
  
  
  Z=data_input[,-1]
  D=data[,disease.name]
  X.strata=as.matrix(X.strata)
  
  
  
  
  loglilke_strat <- function(param){
    
    k=param[1]
    beta=param[2:(nbeta+1)]
    names(beta)=names(fit$coefficients)[-1]
    eta=param[c((length(names(fit$coefficients))+1):(length(names(fit$coefficients))+n_strat))]
    sigma=param[c((length(names(fit$coefficients))+n_strat+1):(length(names(fit$coefficients))+n_strat+n_sigma))]
    sigma[sigma<0]=0.01 #add a condition to make sigma > 0 
    dim(eta) <- c(length(eta), 1)
    dim(sigma) <- c(length(sigma), 1)
    tmp.eta <- as.numeric(as.vector(X.strata %*% eta))
    tmp.sd <- as.numeric(as.vector(sigma.strata %*% sigma))
    
    name_envir_int=colnames(data_input)[c(grep(":",colnames(data_input)))]
    name_rm1=paste0(prs.var,":")
    name_rm2=paste0(":",prs.var)
    name_envir_int=gsub(name_rm1,'',name_envir_int)
    name_envir_int=gsub(name_rm2,'',name_envir_int)
    
    envir_s_int=data_input2[,name_envir_int]
    envir_s_int=as.matrix(envir_s_int)
    beta_int=beta[grep(":",names(beta))]
    dim(beta_int) <- c(length(beta_int), 1)
    beta_p=beta[match(prs.var,names(beta))]
    envir_s=data_input[,-c(1,match(prs.var,colnames(data_input)),grep(":",colnames(data_input)))]
    beta_e_s = beta[-c(match(prs.var,names(beta)),grep(":",names(beta)))]
    dim(beta_e_s) <- c(length(beta_e_s), 1)
    
    
    f=dnorm(prs, mean =  tmp.eta, sd =  abs(tmp.sd), log = F)
    dim(beta) <- c(length(beta), 1)
    vec <- as.numeric(as.vector(exp(D * (k + Z %*% beta )) * f))
    denom_sum=1+exp(tmp.sd^2*(beta_p+ envir_s_int %*% beta_int)^2/2 +tmp.eta*(beta_p+ envir_s_int %*% beta_int) + k + envir_s %*% beta_e_s)
    denom_sum=as.vector(denom_sum)
    vec <- vec/denom_sum
    ret <- sum(log(vec))
    
    return(ret)
  }
  
  
  fn_gr <- function(param){
    
    k=param[1]
    beta=param[2:(nbeta+1)]
    names(beta)=names(fit$coefficients)[-1]
    eta=param[c((length(names(fit$coefficients))+1):(length(names(fit$coefficients))+n_strat))]
    sigma=param[c((length(names(fit$coefficients))+n_strat+1):(length(names(fit$coefficients))+n_strat+n_sigma))]
    dim(eta) <- c(length(eta), 1)
    dim(sigma) <- c(length(sigma), 1)
    tmp.eta <- as.numeric(as.vector(X.strata %*% eta))
    tmp.sd <- as.numeric(as.vector(sigma.strata %*% sigma))
    
    name_envir_int=colnames(data_input)[c(grep(":",colnames(data_input)))]
    name_rm1=paste0(prs.var,":")
    name_rm2=paste0(":",prs.var)
    name_envir_int=gsub(name_rm1,'',name_envir_int)
    name_envir_int=gsub(name_rm2,'',name_envir_int)
    
    envir_s_int=data_input2[,name_envir_int]
    envir_s_int=as.matrix(envir_s_int)
    beta_int=beta[grep(":",names(beta))]
    dim(beta_int) <- c(length(beta_int), 1)
    beta_p=beta[match(prs.var,names(beta))]
    envir_s=data_input[,-c(1,match(prs.var,colnames(data_input)),grep(":",colnames(data_input)))]
    beta_e_s = beta[-c(match(prs.var,names(beta)),grep(":",names(beta)))]
    dim(beta_e_s) <- c(length(beta_e_s), 1)
    
    
    f=dnorm(prs, mean =  tmp.eta, sd =  abs(tmp.sd), log = F)
    dim(beta) <- c(length(beta), 1)
    vec <- as.numeric(as.vector(exp(D * (k + Z %*% beta )) * f))
    denom_sum=1+exp(tmp.sd^2*(beta_p+ envir_s_int %*% beta_int)^2/2 +tmp.eta*(beta_p+ envir_s_int %*% beta_int) + k + envir_s %*% beta_e_s)
    denom_sum=as.vector(denom_sum)
    
    
    k.tmp=sum(D-(1-1/denom_sum))
    beta_p.tmp=sum(D*prs - (1-1/denom_sum)*(tmp.sd^2*as.vector(beta_p+ envir_s_int %*% beta_int)+tmp.eta))
    beta_e_s.tmp=apply(as.matrix(D*envir_s-(1-1/denom_sum)*envir_s),2,sum)
    beta_int.tmp=apply(D*envir_s_int*prs-(1-1/denom_sum)*(tmp.sd^2*as.vector(beta_p+ envir_s_int %*% beta_int)*envir_s_int+tmp.eta*envir_s_int),2,sum)
    
    
    eta.tmp=as.vector(1/tmp.sd^2*(prs-tmp.eta)) %*% X.strata- as.vector(as.vector(1-1/denom_sum)*as.vector(beta_p+ envir_s_int %*% beta_int)) %*% X.strata
    sigma.tmp=as.vector(-1/tmp.sd+(prs-tmp.eta)^2/tmp.sd^3-(1-1/denom_sum)*tmp.sd*as.vector(beta_p+ envir_s_int %*% beta_int)^2) %*% sigma.strata
    
    
    beta.tmp=c(beta_p.tmp,beta_e_s.tmp,beta_int.tmp)
    names(beta.tmp)=c(prs.var,names(fit$coefficients)[-c(1,match(prs.var,names(fit$coefficients)))])
    beta.tmp=beta.tmp[match(names(fit$coefficients)[-1],names(beta.tmp))]
    grr <- c(k.tmp,beta.tmp,eta.tmp,sigma.tmp)
    names(grr)=names(param)
    return(grr)
  }
  
  
  control <- list(fnscale = -1,trace = TRUE,
                  REPORT = 50,maxit=20000)
  
  
  if(numDeriv==T){
    ret <- optim(param0, loglilke_strat, method="BFGS",
                 control = control, hessian = TRUE)
  } else {
    ret <- optim(param0, loglilke_strat,gr=fn_gr, method="BFGS",
                 control = control, hessian = TRUE) }

  
  cov <- chol(-ret$hessian)
  cov <- chol2inv(cov)
  cnames <- names(param0)
  colnames(cov) <- cnames
  rownames(cov) <- cnames
  cov1=cov
  
  
  
  
  
  
  #Summarize the results
  res.sum=function (parms=ret$par, cov=cov1, sided) 
  {
    if (sided != 1) 
      sided <- 2
    cols <- c("Estimate", "Std.Error", "Z.value", "Pvalue")
    n <- length(parms)
    ret <- matrix(data = NA, nrow = n, ncol = 4)
    pnames <- names(parms)
    rownames(ret) <- pnames
    colnames(ret) <- cols
    ret[, 1] <- parms
    cols <- colnames(cov)
    cov <- sqrt(diag(cov))
    names(cov) <- cols
    if (is.null(pnames)) 
      pnames <- 1:n
    cov <- cov[pnames]
    ret[, 2] <- cov
    ret[, 3] <- parms/cov
    ret[, 4] <- sided * pnorm(abs(ret[, 3]), lower.tail = FALSE)
    ret
  }
  
  res_normal=res.sum(sided = side0)
  res=list(res_glm=res_glm,res_normal=res_normal)
  model <- list(data = data,formula=formula,formula_prs=formula_prs,facVar=facVar)
  res$model.info <- model
  return(res)
  
}




#2. Case-only method

summary.caseonly=function (parms, sd, sided = 2) 
{
  if (sided != 1) 
    sided <- 2
  cols <- c("Estimate", "Std.Error", "Z.value", "Pvalue")
  n <- length(parms)
  ret <- matrix(data = NA, nrow = n, ncol = 4)
  pnames <- c("prs",paste0("prs:",names(parms)[-1]))
  rownames(ret) <- pnames
  colnames(ret) <- cols
  ret[, 1] <- parms
  if (is.null(pnames)) 
    pnames <- 1:n
  cov <- sd
  ret[, 2] <- cov
  ret[, 3] <- parms/cov
  ret[, 4] <- sided * pnorm(abs(ret[, 3]), lower.tail = FALSE)
  ret
}

function_caseonly <- function(data_sim,strata=FALSE){
  D=data_sim$D
  prs=data_sim$G
  envir=data_sim$E
  mean_prs=data_sim$mean_prs
  sd_prs=data_sim$sd_prs
  if(strata==TRUE){
    fit_caseonly <- lm(prs[D==1] ~  envir[D==1,1] +envir[D==1,2] + factor(data_sim$S[D==1,1])+data_sim$S[D==1,2])
    
  } else {
    fit_caseonly <- lm(prs[D==1] ~  envir[D==1,1] +envir[D==1,2])
  }
  beta=fit_caseonly$coefficients
  beta_int=fit_caseonly$coefficients[-1]/(sd_prs)^2
  sd_int=summary(fit_caseonly)$coef[-1,2]/(sd_prs)^2
  
  beta_prs=(fit_caseonly$coefficients[1]-mean_prs)/(sd_prs)^2
  sd_prs1= sqrt(( (summary(fit_caseonly)$coef[1,2])^2 +(sd_prs)^2/1000000) / (sd_prs)^4)
  
  
  res=summary.caseonly(parms = c(beta_prs,beta_int),sd=c(sd_prs1,sd_int))
  
  
  beta_int2=fit_caseonly$coefficients[-1]/(sd(fit_caseonly$residuals))^2
  sd_int2=summary(fit_caseonly)$coef[-1,2]/(sd(fit_caseonly$residuals))^2
  
  beta_prs2=(fit_caseonly$coefficients[1]-mean_prs)/(sd(fit_caseonly$residuals))^2
  sd_prs2= sqrt(( (summary(fit_caseonly)$coef[1,2])^2 +(sd(fit_caseonly$residuals))^2/1000000) / (sd(fit_caseonly$residuals))^4)
  res2=summary.caseonly(parms = c(beta_prs2,beta_int2),sd=c(sd_prs2,sd_int2))
  rownames(res2)=paste0(rownames(res2),"_resid_sd")
  res=rbind(res,res2)
  
  return(res)
}
#Example:
#source("/users/zwang4/GXE/sim_data_function.R")
#dat=simFit(strata = T)
#dat1=simFit(strata = F)
#test=function_caseonly(data_sim = dat,strata = T)
#test2=function_caseonly(data_sim = dat1)


