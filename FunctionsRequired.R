####################################################################
# Library
####################################################################
library(copula)
library(xlsx)
library(ggplot2)
library(gridExtra)
####################################################################
# Density function of FGM copula.
Density_FGM=function(u,v,varphi){
  return(1+varphi*(1-2*u)*(1-2*v))
}
# Likelihood function
likelihood <- function(u, v, varphi) {
  return(prod(Density_FGM(u, v, varphi)))
}

# Derivative of the log-likelihood function.
# der: order of derivative.
Der_log_lik_FGM=function(varphi,u,v,der=1){
  dln=NA;  n=length(u)
  for (i in 1:n) {
    num=(1-2*u[i])*(1-2*v[i])
    den=1+varphi*num
    if(der==1){
      dln[i]=(num/den)
    }else if(der==2){
      dln[i]=-(num/den)^2
    }
  }
  return(sum(dln))
}

####################################################################
## Data simulation
####################################################################

# Generation of one random observation of FGM copula
rfgm=function(varphi){
  u1=runif(n = 1,0,1) 
  v2=runif(n = 1,0,1) # conditional probability
  A=varphi*(2*u1-1)-1
  B=sqrt(1-2*varphi*(2*u1-1)+varphi^2*(2*u1-1)^2+4*varphi*v2*(2*u1-1))
  u2=2*v2/(B-A)
  return(c(u1,u2))
}


## Generation of n random observation of FGM copula

rnfgm=function(varphi,n){
  datafgm=matrix(data = NA,nrow = 0,ncol = 2,dimnames = list(c(),c("u","v")))
  for (i in 1:n) {
    datafgm=rbind(datafgm,rfgm(varphi))
  }
  return(datafgm)
}
# Example
#dim(rnfgm(0.25,10))[1]

####################################################################
####################################################################
####################################################################

####################################################################
## Bayesian Inference
####################################################################

# Density of the Beta prior distribution
prior_density_beta <- function(phi, alpha, beta) {
  return((phi + 1)^(alpha - 1) * (1 - phi)^(beta - 1))
}

# Density of the Uniform prior distribution
prior_density_uniform <- function(phi) {
  return(1 / 2)  
}

# Density of the triangular prior distribution
prior_density_triangular <- function(phi, c) {
  if (phi >= -1 && phi < c) {
    return((phi + 1) / (c + 1))
  } else if (phi >= c && phi <= 1) {
    return((1 - phi) / (1 - c))
  } else {
    return(0)
  }
}

# Metropolis-Hastings Algorithm
# u and v are the sample observations.
# alpha_prior, beta_prior, and c_prior are the values of the hyperparameters.
# instrumental: determines the instrumental distribution that can be used, it can be "beta" or "uniform".
# precision: works when the instrumental distribution is the beta distribution.
# n_iter: length of the generated chain.
# init_varphi: seed with which the chain starts.

metropolis_hastings <- function(u, v, prior, alpha_prior=NA, beta_prior=NA, c_prior=NA, 
                                instrumental="beta",precision=20,
                                n_iter, init_varphi=runif(1,-1,1)) {
  varphi <- numeric(n_iter)
  varphi[1] <- init_varphi
  acceptance_rate <- 0
  
  for (i in 2:n_iter) {
    if(instrumental=="beta"){
      mu_current <- (varphi[i-1] + 1) / 2  # Transform current varphi to range [0, 1]
      # Calculate parameters of the instrumental Beta distribution
      alpha_inst <- mu_current * precision
      beta_inst <- (1 - mu_current) * precision
      
      # Proposed value of the instrumental Beta distribution
      varphi_proposal <- rbeta(1, alpha_inst, beta_inst)
      varphi_proposal <- 2 * varphi_proposal - 1  # Transform proposed to range [-1, 1]
      
      # Calculate parameters of the instrumental Beta distribution for the proposed.
      mu_proposal <- (varphi_proposal + 1) / 2
      alpha_inst_proposal <- mu_proposal * precision
      beta_inst_proposal <- (1 - mu_proposal) * precision
      
      # Density of the beta instrumental distribution with current and proposed value of varphi.
      proposal_density_current_given_proposal <- dbeta(mu_current, alpha_inst_proposal, beta_inst_proposal)
      proposal_density_proposal_given_current <- dbeta(mu_proposal, alpha_inst, beta_inst)
    } else if (instrumental=="uniform"){
      varphi_proposal=runif(n=1,-1,1)
      proposal_density_current_given_proposal <- 1/2
      proposal_density_proposal_given_current <- 1/2
    }
    
    # Selection of the a prior distribution
    if (prior == "beta") {
      prior_density_current <- prior_density_beta(varphi[i-1], alpha_prior, beta_prior)
      prior_density_proposal <- prior_density_beta(varphi_proposal, alpha_prior, beta_prior)
    } else if (prior == "uniform") {
      prior_density_current <- prior_density_uniform(varphi[i-1])
      prior_density_proposal <- prior_density_uniform(varphi_proposal)
    } else if (prior == "triangular") {
      prior_density_current <- prior_density_triangular(varphi[i-1], c_prior)
      prior_density_proposal <- prior_density_triangular(varphi_proposal, c_prior)
    }
    
    # Calculate the acceptance ratio
    acceptance_ratio <- (likelihood(u, v, varphi_proposal) * prior_density_proposal) / 
      (likelihood(u, v, varphi[i-1]) * prior_density_current)
    
    # Adjust the acceptance ratio with the instrumental density
    acceptance_ratio <- acceptance_ratio * proposal_density_current_given_proposal / proposal_density_proposal_given_current
    
    # Identification of errors
    if (is.nan(acceptance_ratio) || is.infinite(acceptance_ratio)) {
      stop(paste("Non-numeric acceptance_ratio detected at iteration",i, 
                 "with proposal", varphi_proposal, 
                 "and previous chain value", varphi[i - 1],"the acceptance ration is", acceptance_ratio, 
                 "the density current value is",proposal_density_current_given_proposal,
                 "the density proposal value is", proposal_density_proposal_given_current,
                 "the value of the parameter proposal and current are", alpha_inst_proposal, beta_inst_proposal,
                 alpha_inst, beta_inst))
    }
    if (runif(1) < acceptance_ratio) {
      varphi[i] <- varphi_proposal
      acceptance_rate <- acceptance_rate + 1
    } else {
      varphi[i] <- varphi[i-1]
    }
  }
  
  acceptance_rate <- acceptance_rate / (n_iter - 1)
  return(list(varphi = varphi, acceptance_rate = acceptance_rate))
}

####################################################################
####################################################################
####################################################################

####################################################################
# Moments estimate
####################################################################

me=function(dataset_target){
  pearson=cor(dataset_target,method = "pearson")[1,2]*3
  kendall=cor(dataset_target,method = "kendall")[1,2]*9/2
  spearman=cor(dataset_target,method = "spearman")[1,2]*3
  return(list(varphipearson=pearson,varphikendall=kendall,varphispearman=spearman))
}

####################################################################
# Estimates (point estimates, interval estimates, and credibility regions) 
# produced from a random sample.
####################################################################
# varphi: Dependence parameter.
# n: Sample size.
# nboot: Number of resamples performed in the bootstrap procedure.
# confidence: Desired confidence level for interval estimates; the probability 
# that a new sample will produce estimates within this interval.
# alpha_prior and beta_prior: Hyperparameters of the prior Beta distribution.
# c_prior: Hyperparameter of the prior Triangular distribution.
# n_thin: Thinning interval for the MCMC chain (every nth iteration is kept).

# Parameter estimation using the dataset_target sample.
Est_A_Samp=function(dataset_target,nboot=100,confidence = 0.95,alpha_prior, beta_prior, c_prior,
                    n_burne=1000,n_iter=5000,n_thin=1,type_sample="all"){
  ##############################################
  # Point estimates for the FGM dependence parameter
  ##############################################
  n=dim(dataset_target)[1]
  # Estimation of maximum log-likelihood
  MLE=safeUroot(f=function(ph){mapply(function(ph){
    Der_log_lik_FGM(ph,dataset_target[,1],dataset_target[,2])},ph)},c(-1,1))$root
  # Moment estimate
  EM=me(dataset_target)
  
  # Estimates using Bayes approach.
  ResTriangular=metropolis_hastings(u=dataset_target[,1], v=dataset_target[,2], 
                                    prior="triangular",c_prior = c_prior,instrumental="uniform",
                                    precision=20,n_iter=n_iter)$varphi
  ResBeta=metropolis_hastings(u=dataset_target[,1], v=dataset_target[,2], 
                              prior="beta", alpha_prior, beta_prior, 
                              instrumental="uniform",precision=20,n_iter=n_iter)$varphi
  ResUnif=metropolis_hastings(u=dataset_target[,1], v=dataset_target[,2], 
                              prior="uniform",instrumental="uniform",precision=20,n_iter=n_iter)$varphi
  
  ResTriangular=ResTriangular[seq((n_burne+1),n_iter,by=n_thin)]
  ResBeta=ResBeta[seq((n_burne+1),n_iter,by=n_thin)]
  ResUnif=ResUnif[seq((n_burne+1),n_iter,by=n_thin)]
  ##############################################
  # Interval estimation
  ##############################################
  
  alpha=1-confidence
  # Bootstrap interval
  vphitau=c();vphispe=c();vphilm=c()  #Bootstrap estimates for the dependence parameter using classical estimators.
  i=1 # bootstrap counter.
  
  if(type_sample=="all"){
    while(i<=nboot) {
      cho=sample(1:n,n,replace = T)     # Resampling the position in the sample.
      rdatan=dataset_target[cho,]       # Selecting elements from the sample.
      #Estimates
      philm=safeUroot(f=function(ph){mapply(function(ph){
        Der_log_lik_FGM(ph,rdatan[,1],rdatan[,2])},ph)},c(-1,1))$root
      
      vphilm[i]=philm
      vphim=me(rdatan)
      vphitau[i]=vphim$varphikendall
      vphispe[i]=vphim$varphispearman
      i=i+1
    }} else if (type_sample=="exclude"){
      ie=0 # counter of the number of times bootstrap resampling produces estimates outside [-1,1]. 
      while(i<=nboot) {
        cho=sample(1:n,n,replace = T)     # Resampling the position in the sample.
        rdatan=dataset_target[cho,]       # Selecting elements from the sample.
        #Estimates
        philm=safeUroot(f=function(ph){mapply(function(ph){
          Der_log_lik_FGM(ph,rdatan[,1],rdatan[,2])},ph)},c(-1,1))$root
        if(all(abs(philm)<=1)){
          vphilm[i]=philm
          vphim=me(rdatan)
          if(all(abs(vphim$varphikendall)<=1,abs(vphim$varphispearman)<=1)){
            vphitau=cbind(vphitau,vphim$varphikendall)
            vphispe=cbind(vphispe,vphim$varphispearman)
            i=i+1
          }}else{ie=ie+1}
      }}
  
  IntervalML=quantile(vphilm,probs=c(alpha/2,1-alpha/2))   # Quantile for estimates using Maximum likelihood.
  IntervalTau=quantile(vphitau,probs=c(alpha/2,1-alpha/2)) # Quantile for estimates using moment method.
  IntervalSpe=quantile(vphispe,probs=c(alpha/2,1-alpha/2)) # Quantile for estimates using moment method.
  
  # Asymptotic (1-alpha)100% confidence interval for the MLE
  variance_MLE=1/Der_log_lik_FGM(varphi = MLE,u = dataset_target[,1],v = dataset_target[,2],der = 2)
  Interval_A_MLE= MLE+c(qnorm(alpha/2),qnorm(1-alpha/2))*sqrt(abs(variance_MLE))
  
  # Asymptotic (1-alpha)100% confidence interval for the Tau Kendall
  Interval_A_Tau= EM$varphikendall+9/2*c(qnorm(alpha/2),qnorm(1-alpha/2))*sqrt(2*(2*n+5)/(9*n*(n-1)))
  
  # Asymptotic (1-alpha)100% confidence interval for the Spearman 
  Interval_A_Spe= EM$varphispearman+3*c(qnorm(alpha/2),qnorm(1-alpha/2))*sqrt(1/(n-1))
  
  # Credibility regions
  Cred_Reg_T=quantile(ResTriangular,probs=c(alpha/2,1-alpha/2))
  Cred_Reg_B=quantile(ResBeta,probs=c(alpha/2,1-alpha/2))
  Cred_Reg_U=quantile(ResUnif,probs=c(alpha/2,1-alpha/2))
  
  Summary=matrix(data = c(MLE,IntervalML,Interval_A_MLE,
                          EM$varphikendall,IntervalTau,Interval_A_Tau,
                          EM$varphispearman,IntervalSpe,Interval_A_Spe,
                          mean(ResTriangular),Cred_Reg_T,NA,NA,
                          mean(ResBeta),Cred_Reg_B,NA,NA,
                          mean(ResUnif),Cred_Reg_U,NA,NA),nrow =6 ,ncol =5,byrow = TRUE,
                 dimnames = list(list("LM","Tau","Spe","Triangular","Beta","Unif"),
                                 list("Mean","Low","Upper","AL","AU")))
  
  return(list(MLE=MLE,EM=EM,EBT=mean(ResTriangular),EBB=mean(ResBeta),EBU=mean(ResUnif),
              IML=IntervalML,ITau=IntervalTau,ISpe=IntervalSpe,CRT=Cred_Reg_T,IAMLE=Interval_A_MLE,
              IATAU=Interval_A_Tau,IASPE=Interval_A_Spe,
              CRB=Cred_Reg_B, CRU=Cred_Reg_U, vphilm=vphilm,
              vphitau=vphitau,vphispe=vphispe,Sims_T=ResTriangular,
              Sims_B=ResBeta,Sims_U=ResUnif,n_Chain_Length=length(ResTriangular),
              Summary=Summary))
}

#Example
#trial0=rnfgm(-0.75,600)
#trial1=Est_A_Samp(trial0,alpha_prior = 12,beta_prior = 14,c_prior = -0.75,type_sample = "exclude")
#trial1$Summary

# Parameter estimation using a randomly generated sample of size n and varphi dependence.
Sim_Est_A_Samp=function(varphi,n,nboot=100,confidence=0.95,alpha_prior, beta_prior, c_prior,
                        n_burne=1000,n_iter=5000,n_thin=1){
  dataset_target=rnfgm(varphi,n) # random sample
  
  # Estimation of maximum log-likelihood
  MLE=safeUroot(f=function(ph){mapply(function(ph){
    Der_log_lik_FGM(ph,dataset_target[,1],dataset_target[,2])},ph)},c(-1,1))$root
  # Moment estimate
  EM=me(dataset_target)
  
  # Selecting samples that have Maximum Likelihood and Moments estimates 
  # within the domain of the dependency parameter of the FGM copula.
  while(any(abs(MLE) > 1, abs(EM$varphikendall) > 1, abs(EM$varphispearman) > 1)){
    dataset_target=rnfgm(varphi,n)
    MLE=safeUroot(f=function(ph){mapply(function(ph){
      Der_log_lik_FGM(ph,dataset_target[,1],dataset_target[,2])},ph)},c(-1,1))$root
    EM=me(dataset_target)
  }
  results_Sim=Est_A_Samp(dataset_target,nboot,confidence,alpha_prior, beta_prior, c_prior,
                         n_burne,n_iter,n_thin,type_sample="exclude")
  return(results_Sim)
}

#Example
#tiral1=Sim_Est_A_Samp(-0.75,100,nboot=100,confidence = 0.95,alpha_prior=12, beta_prior=14, c_prior=-0.75,
#                       n_burne=1000,n_iter=5000,n_thin=2)

####################################################################
# Simulation study for evaluating the characteristics (mean squared error, 
# minimum, maximum, average length, coverage probability) of
# estimators from MLE, MM, and Bayesian approaches (Triangular, 
# Beta, Uniform).
####################################################################

# varphi: Dependence parameter.
# n: Sample size.
# N: Number of FGM samples generated.
# nboot: Number of resamples performed in the bootstrap procedure.
# confidence: Desired confidence level for interval estimates; the probability 
# that a new sample will produce estimates within this interval.
# alpha_prior and beta_prior: Hyperparameters of the prior Beta distribution.
# c_prior: Hyperparameter of the prior Triangular distribution.
# n_burne: Number of iterations discarded at the beginning of the MCMC chain (burn-in period).
# n_iter: Total number of iterations in the MCMC chain.
# n_thin: Thinning interval for the MCMC chain (every nth iteration is kept).

Sim_Est_N_Samp=function(varphi,n,nboot=100, N, confidence,alpha_prior,beta_prior,c_prior,n_burne=1000,
                        n_iter=5000,n_thin=1){
  PhiML=NA; PhiTau=NA; PhiSpe=NA; PhiBT=NA; PhiBB=NA; PhiBU=NA; Error_Sample_Boot=NA
  IntervalML=matrix(data=NA,nrow=0,ncol=2,dimnames = list(c(),c("L","R")))
  IntervalTau=IntervalML; IntervalSpe=IntervalML; IntervalIAMLE=IntervalML
  IntervalIATAU=IntervalML;IntervalIASPE=IntervalML
  Cred_Reg_T=IntervalML;Cred_Reg_B=IntervalML;Cred_Reg_U=IntervalML
  i=1
  while (i<=N) {
    results=Sim_Est_A_Samp(varphi,n,nboot,confidence,alpha_prior,beta_prior,c_prior,
                           n_burne,n_iter,n_thin)
    # Estimates
    PhiML[i]=results$MLE; PhiTau[i]=results$EM$varphikendall; PhiSpe[i]=results$EM$varphispearman
    PhiBT[i]=results$EBT; PhiBB[i]=results$EBB; PhiBU[i]=results$EBU
    # Intervals estimation
    IntervalML=rbind(IntervalML,results$IML)
    IntervalIAMLE=rbind(IntervalIAMLE,results$IAMLE)
    IntervalTau=rbind(IntervalTau,results$ITau)
    IntervalIATAU=rbind(IntervalIATAU,results$IATAU)
    IntervalSpe=rbind(IntervalSpe,results$ISpe)
    IntervalIASPE=rbind(IntervalIASPE,results$IASPE)
    Cred_Reg_T=rbind(Cred_Reg_T,results$CRT)
    Cred_Reg_B=rbind(Cred_Reg_B,results$CRB)
    Cred_Reg_U=rbind(Cred_Reg_U,results$CRU)
    #Error_Sample_Boot[i]=results$Er_S_B
    i=i+1}
  # Descriptive measures
  DesML=c(min(PhiML),mean(PhiML),sd(PhiML),max(PhiML))
  DesMTau=c(min(PhiTau),mean(PhiTau),sd(PhiTau),max(PhiTau))
  DesMSpe=c(min(PhiSpe),mean(PhiSpe),sd(PhiSpe),max(PhiSpe))
  DesBT=c(min(PhiBT),mean(PhiBT),sd(PhiBT),max(PhiBT))
  DesBB=c(min(PhiBB),mean(PhiBB),sd(PhiBB),max(PhiBB))
  DesBU=c(min(PhiBU),mean(PhiBU),sd(PhiBU),max(PhiBU))
  # Estimator bias
  Bias=c(varphi-DesML[2],varphi-DesMTau[2],varphi-DesMSpe[2],varphi-DesBT[2],varphi-DesBB[2],varphi-DesBU[2])
  
  # Average length 
  meanlenght=function(Interval){
    lenint=sum(Interval[,2]-Interval[,1])/nrow(Interval)
    return(lenint=lenint)
  }
  
  lenghtILM=meanlenght(IntervalML)
  lenghtIAMLE=meanlenght(IntervalIAMLE)
  lenghtITau=meanlenght(IntervalTau)
  lenghtIATAU=meanlenght(IntervalIATAU)
  lenghtISpe=meanlenght(IntervalSpe)
  lenghtIASPE=meanlenght(IntervalIASPE)
  lenghtRCT=meanlenght(Cred_Reg_T)
  lenghtRCB=meanlenght(Cred_Reg_B)
  lenghtRCU=meanlenght(Cred_Reg_U)
  
  # Coverage probability
  covprob=function(Interval){
    indicator=function(j){
      if(Interval[j,1]<=varphi & Interval[j,2]>=varphi){
        return(1)
      }else{return(0)}
    }
    resind=NA
    for(l in 1:nrow(Interval)){
      resind[l]=indicator(l)
    }
    covint=sum(resind)/nrow(Interval)
    
    return(covint=covint)
  }
  covprobILM=covprob(IntervalML)
  covprobIAMLE=covprob(IntervalIAMLE)
  covprobITau=covprob(IntervalTau)
  covprobIATAU=covprob(IntervalIATAU)
  covprobISpe=covprob(IntervalSpe)
  covprobIASPE=covprob(IntervalIASPE)
  covprobRCT=covprob(Cred_Reg_T)
  covprobRCB=covprob(Cred_Reg_B)
  covprobRCU=covprob(Cred_Reg_U)
  
  Descritive=matrix(data=c(c(DesML,Bias[1],var(PhiML)+(Bias[1])^2,lenghtILM,covprobILM,lenghtIAMLE,covprobIAMLE),
                           c(DesMTau,Bias[2],var(PhiTau)+(Bias[2])^2,lenghtITau,covprobITau,lenghtIATAU,covprobIATAU),
                           c(DesMSpe,Bias[3],var(PhiSpe)+(Bias[3])^2,lenghtISpe,covprobISpe,lenghtIASPE,covprobIASPE),
                           c(DesBT,Bias[4],var(PhiBT)+(Bias[4])^2,lenghtRCT,covprobRCT,NA,NA),
                           c(DesBB,Bias[5],var(PhiBB)+(Bias[5])^2,lenghtRCB,covprobRCB,NA,NA),
                           c(DesBU,Bias[6],var(PhiBU)+(Bias[6])^2,lenghtRCU,covprobRCU,NA,NA)),nrow=10,ncol=6,
                    dimnames = list(list("Min","Mean","SD","Max","Bias","MSE","BootLength","BootCoverage","ALength","ACoverage"),
                                    list("ML","Tau","Spe","Triangular","Beta","Unif")))
  return(Descritive)
}

#Example
#trial1=Sim_Est_N_Samp(varphi=-0.75,n=10,nboot=100, N=10, confidence=0.95,alpha_prior=12,beta_prior=14,
#               c_prior=-0.75,n_burne=1000,n_iter=5000,n_thin=1)
#trial1

###################################################################
# Obtaining hyperparameter values, modified Tovar's method
####################################################################
# x1, x2: Quantiles obtained from an elicited process.
# 1-alpha: Confidence level that the interval (x1, x2) contains the true value of the parameter.
# low: Lower limit of the Beta distribution.
# upp: Upper limit of the Beta distribution.
Mtovar=function(x1,x2,alp,low=-1,upp=1){
  tht0=(x1+x2)/2
  w=(tht0-low)/(upp-tht0)
  sig=sqrt(alp)*(x1-tht0)
  b= ((upp-low)^2*w-((w+1)^2*sig^2))/((w+1)^3*sig^2)
  a=w*b
  return(list(a=a,b=b,c=tht0))
}

################################################################################
## Graph construction
################################################################################
# m: Descriptive measure to be graphed.
# omega: Value of the dependence parameter varphi.
# omega_to_file: List containing the values of varphi and corresponding file names
# in .xlsx format. For example, omega_to_file = list("-0.75" = "Results-075.xlsx",
# "-0.25" = "Results-025.xlsx").
# useDependence: If TRUE, constructs graphs for the method cho_method = 7, 
# specifically for the Beta prior. Otherwise, constructs graphs for all methods.
# useAsymptotic: If TRUE, graphs asymptotic intervals of ML estimates, if FALSE, these intervals are not graphed.
# cho_method: The method to use for constructing the graphs. Default is 6.
GeneralGraph = function(m, omega = NULL, omega_to_file = NULL, useDependence = FALSE,useAsymtotic=FALSE,cho_method=6) {
  nameGra = c("Min", "Mean", "SD", "Max", "Bias","MSE", "Length", "Coverage")
  colName = ifelse(useDependence, "Dependence", "Method")
  if(useDependence) {labelSet =c("SND","WND", "WPD", "SPD")} else if(useAsymtotic) {
    labelSet = c("ML", "Tau-Kendall", "Spearman", "Triangular", "Beta", "Uniform","ML A","Tau A","Spe A")} else
    {labelSet = c("ML", "Tau-Kendall", "Spearman", "Triangular", "Beta", "Uniform")}
  namemethod = 1:length(labelSet)
  
  MA = matrix(data = NA, nrow = 0, ncol = 3, dimnames = list(NULL, c("Size", nameGra[m], colName)))
  size = c(10, seq(50, 1000, 50))
  
  loop_omega = if(useDependence) names(omega_to_file) else omega
  
  for(w in loop_omega) {
    filename = if(useDependence) {
      omega_to_file[[w]]
    } else {
      omega_to_file[[as.character(w)]]
    }
    
    if(is.null(filename)) {
      stop("Invalid value for varphi")
    }
    
    for(i in seq_along(size)) {
      A = read.xlsx(file = filename, sheetIndex = i)
      columnIndex = if(useDependence) cho_method else if(useAsymtotic) 2:10 else 2:7
      
      for(j in columnIndex) {
        if(j<=7){
          MA = rbind(MA, c(size[i], A[m, j], if(useDependence) w else if(useAsymtotic) namemethod[(j - 1)] else namemethod[(j - 1)]))
        }else if(j>7&&useAsymtotic){
          MA = rbind(MA, c(size[i], A[m+2, j-7+1], if(useDependence) w else if(useAsymtotic) namemethod[(j - 1)]))
        }
      }
    }
  }
  MA=apply(MA, 2, as.numeric)
  MA = as.data.frame(MA)
  MA[[colName]] = factor(MA[[colName]], labels = labelSet)
  return(MA)
}


# Example 

# omega_to_file = list(
#   "-0.75" = "2023Results-075.xlsx",
#   "-0.25" = "2023Results-025.xlsx",
#   "0.25" = "2023Results025.xlsx",
#   "0.75" = "2023Results075.xlsx")

# graph of all methods
# Bias=GeneralGraph(5, omega = "-0.75", omega_to_file = omega_to_file, useDependence = FALSE)
# Bias
# ggplot(Bias,aes(x=Size , y=Bias,color = Method))+geom_line()+labs(title="Dependence of -0.75")+
#   scale_x_continuous(breaks=c(10,seq(100,1000,100)))+xlab("Sample size")
# 
# # only graph descriptive measures obtained with the a priori Beta
# Bias=GeneralGraph(5, omega = NULL, omega_to_file = omega_to_file, useDependence = TRUE)
# Bias
# ggplot(Bias,aes(x=Size , y=Bias,color = Dependence))+geom_line()+labs(title=" ")+
#   scale_x_continuous(breaks=c(10,seq(100,1000,100)))+xlab("Sample size")



######################################################################
# Empirical marginal distribution
######################################################################
# c: Column index of the variable for which the empirical marginal distribution is computed.
# X: Data matrix where each row represents an observation and each column represents a variable.

Marginal <- function(c, X) {
  # Helper function to check if w is less than or equal to z.
  step.my <- function(w, z) {
    if (w <= z) {
      1
    } else {
      0
    }
  }
  
  N <- nrow(X)  # Number of observations in the data matrix
  fn <- NULL  # Initialize the empirical distribution vector
  
  for (i in 1:nrow(X)) {
    d <- 0  # Initialize the count for the current observation
    for (j in 1:nrow(X)) {
      e <- step.my(X[j, c], X[i, c])  # Apply the step function
      d <- e + d  # Accumulate the count
    }
    fn[i] <- 1 / (N + 1) * d  # Compute the empirical distribution value for the current observation
  }
  
  fn  # Return the empirical distribution vector
}
