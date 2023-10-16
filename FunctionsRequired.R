####################################################################
# Library
####################################################################
library(R2OpenBUGS)
library(copula)
library(xlsx)
library(ggplot2)
library(gridExtra)
####################################################################
# Density function of FGM copula.
Density_FGM=function(u,v,varphi){
  c=1+varphi*(1-2*u)*(1-2*v)
  return(c)
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
  sum(dln)
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
## Definition of models to be use in OpenBugs
####################################################################

# Triangular distribution as prior distribution
sink("modelT.txt")
cat("model
      {  
      const <-10000
      for(i in 1 : n) {
      z[i] <- 1
      z[i] ~ dbern(p[i])
      p[i] <- ( 1 + varphi*(1-2*u[i])*(1-2*v[i]) ) / const
      }
      y <- 0
      y ~ dloglik(rho)
      indicator1 <- step(varphi-c); indicator2 <- step(1+varphi); indicator3 <- step(1-varphi)
      LogP <- indicator2*(1-indicator1)*(varphi+1)/(c+1) + indicator1*indicator3*(1-varphi)/(1-c)
      varphi ~ dflat()
      rho <- log(LogP)
      }",fill=TRUE)
sink()

# Beta distribution as prior distribution.
sink("modelB.txt")
cat("model
      {  
      for(i in 1 : n) {
      z[i] <- 0
      z[i] ~ dloglik(logL[i])
      logL[i] <- log( 1 + varphi*(1-2*u[i])*(1-2*v[i]) )
      }
      theta ~ dbeta(a,b)
      varphi <- 2*theta-1
      }",fill=TRUE)
sink()

# Uniform distribution as prior distribution.
sink("modelU.txt")
cat("model
      {  
      for(i in 1 : n) {
      z[i] <- 0
      z[i] ~ dloglik(logL[i])
      logL[i] <- log( 1 + varphi*(1-2*u[i])*(1-2*v[i]) )
      }
      varphi ~ dunif(-1,1)
      }",fill=TRUE)
sink()

####################################################################
# Algorithmic to estimate the parameter varphi using Bayesian Method
####################################################################
# dataset_target: matrix of two columns.
# a and b are hiperparameters of prior beta distribution.
# c is the hiperparameters of prior triangular distribution.
# model: file name .txt containing the post_est distribution in OpenBugs format
post_est=function(dataset_target,model,a=0,b=0,c=0,n_burnd=1000,n_iter=5000,n_thin=1,n_chain=1){
  parameter=c("varphi")
  if(model=="modelU.txt"){ # Uniform distribution
    Data <- list (u=dataset_target[,1], v=dataset_target[,2], n=nrow(dataset_target))
    Inits=NULL
  }else if(model=="modelB.txt"){ # Beta distribution
    Inits=function(){list(varphi=2*rbeta(1,1,1)-1)}
    Data <- list (u=dataset_target[,1], v=dataset_target[,2], n=nrow(dataset_target),a=a, b=b)
  }else if(model=="modelT.txt"){ # Triangular  distribution
    Inits=function(){list(varphi=runif(1,-1,1))}
    Data <- list (u=dataset_target[,1], v=dataset_target[,2], n=nrow(dataset_target), c=c)
  }
  results.sim <- bugs(data=Data, inits=Inits, parameters = parameter,
                   model.file=model, n.iter=n_iter, n.burnin = n_burnd,
                   n.chains=n_chain,n.thin = n_thin)
  return(results.sim)
}

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
# Estimates (point, interval, credibility region) produced by a random sample.
####################################################################
# varphi: dependence parameter
# n: sample size
# nboot: Number of resamples performed in the bootstrap.
# confidence: Probability that a new sample produces estimates within an interval. 
# a and b are hiperparameters of prior beta distribution.
# c is the hiperparameters of prior triangular distribution.

Est_A_Samp=function(dataset_target,nboot=200,confidence,a,b,c,n_burnd=1000,n_iter=5000,n_thin=1,n_chain=1){
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
  ResTriangular=post_est(dataset_target,"modelT.txt", c=c,n_burnd=n_burnd,n_iter=n_iter,n_thin=n_thin,n_chain=n_chain)
  ResBeta=post_est(dataset_target,"modelB.txt",a=a, b=b,n_burnd=n_burnd,n_iter=n_iter,n_thin=n_thin,n_chain=n_chain)
  ResUnif=post_est(dataset_target,"modelU.txt",n_burnd=n_burnd,n_iter=n_iter,n_thin=n_thin,n_chain=n_chain)
  ##############################################
  # Interval estimation
  ##############################################
  
  alpha=1-confidence
  
  # Bootstrap interval
  vphitau=c();vphispe=c();vphilm=c()  #Bootstrap estimates for the dependence parameter using classical estimators.
  i=1 # bootstrap counter.
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
  }
  IntervalML=quantile(vphilm,probs=c(alpha/2,1-alpha/2))   # Quantile for estimates using Maximum likelihood.
  IntervalTau=quantile(vphitau,probs=c(alpha/2,1-alpha/2)) # Quantile for estimates using moment method.
  IntervalSpe=quantile(vphispe,probs=c(alpha/2,1-alpha/2)) # Quantile for estimates using moment method.
  
  # Asymptotic (1-alpha)100% confidence interval for the MLE
  variance_MLE=1/Der_log_lik_FGM(varphi = MLE,u = dataset_target[,1],v = dataset_target[,2],der = 2)
  Interval_A_MLE= MLE+c(qnorm(alpha/2),qnorm(1-alpha/2))*sqrt(abs(variance_MLE))
  
  # Credibility regions
  Cred_Reg_T=c(ResTriangular$summary[1,3],ResTriangular$summary[1,7])
  Cred_Reg_B=c(ResBeta$summary[1,3],ResBeta$summary[1,7])
  Cred_Reg_U=c(ResUnif$summary[1,3],ResUnif$summary[1,7])
  
  return(list(MLE=MLE,EM=EM,EBT=ResTriangular$mean$varphi,EBB=ResBeta$mean$varphi,EBU=ResUnif$mean$varphi,
              IML=IntervalML,ITau=IntervalTau,ISpe=IntervalSpe,CRT=Cred_Reg_T,IAMLE=Interval_A_MLE,
              CRB=Cred_Reg_B, CRU=Cred_Reg_U, vphilm=vphilm,
              vphitau=vphitau,vphispe=vphispe,Er_S_B=ie,summ_Trian=ResTriangular$summary,
              summ_Beta=ResBeta$summary,summ_Unif=ResUnif$summary,DIC_T=ResTriangular$DIC,
              DIC_B=ResBeta$DIC,DIC_U=ResUnif$DIC,Sims_T=ResTriangular$sims.list$varphi,
              Sims_B=ResBeta$sims.list$varphi,Sims_U=ResUnif$sims.list$varphi))
}

#Example
#Est_A_Samp(dataset_target,10,nboot=200,confidence = 0.95,a=1,b=1,c=0,n_burnd=1000,n_iter=5000,n_thin=10,n_chain=1)

Sim_Est_A_Samp=function(varphi,n,nboot=200,confidence,a,b,c,n_burnd=1000,n_iter=5000,n_thin=1,n_chain=1){
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
  results_Sim=Est_A_Samp(dataset_target,nboot=200,confidence,a,b,c,n_burnd=1000,n_iter=5000,n_thin=1,n_chain=1)
  return(results_Sim)
}
#Example
#Sim_Est_A_Samp(0.25,10,nboot=200,confidence = 0.95,a=1,b=1,c=0,n_burnd=1000,n_iter=5000,n_thin=10,n_chain=1)
                    
####################################################################
# Simulation study for the characteristics (mean squared error, 
# minimum, maximum, average length, coverage probability) of
# estimators from MLE, MM, and Bayesian approaches (Triangular, 
# Beta, Uniform)
####################################################################

# varphi: dependence parameter
# n: sample size
# N: Number of FGM samples generated.
# nbood: Number of resamples performed in the bootstrap.
# confidence: Probability that a new sample produces estimates within an interval. 
# a and b are hiperparameters of prior beta distribution.
# c is the hiperparameters of prior triangular distribution.

Sim_Est_N_Samp=function(varphi,n,nboot=200, N, confidence,a,b,c,n_burnd=1000,n_iter=5000,n_thin=1,n_chain=1){
  PhiML=NA; PhiTau=NA; PhiSpe=NA; PhiBT=NA; PhiBB=NA; PhiBU=NA; Error_Sample_Boot=NA
  IntervalML=matrix(data=NA,nrow=0,ncol=2,dimnames = list(c(),c("L","R")))
  IntervalTau=IntervalML; IntervalSpe=IntervalML; IntervalIAMLE=IntervalML
  Cred_Reg_T=IntervalML;Cred_Reg_B=IntervalML;Cred_Reg_U=IntervalML
  i=0
  while (i<=N) {
    results=Sim_Est_A_Samp(varphi,n,nboot,confidence,a,b,c,n_burnd=n_burnd,n_iter=n_iter,n_thin=n_thin,n_chain=n_chain)
    # Estimates
    PhiML[i]=results$MLE; PhiTau[i]=results$EM$varphikendall; PhiSpe[i]=results$EM$varphispearman
    PhiBT[i]=results$EBT; PhiBB[i]=results$EBB; PhiBU[i]=results$EBU
    # Intervals estimation
    IntervalML=rbind(IntervalML,results$IML)
    IntervalIAMLE=rbind(IntervalIAMLE,results$IAMLE)
    IntervalTau=rbind(IntervalTau,results$ITau)
    IntervalSpe=rbind(IntervalSpe,results$ISpe)
    Cred_Reg_T=rbind(Cred_Reg_T,results$CRT)
    Cred_Reg_B=rbind(Cred_Reg_B,results$CRB)
    Cred_Reg_U=rbind(Cred_Reg_U,results$CRU)
    Error_Sample_Boot[i]=results$Er_S_B
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
  lenghtISpe=meanlenght(IntervalSpe)
  lenghtRCT=meanlenght(Cred_Reg_T)
  lenghtRCB=meanlenght(Cred_Reg_B)
  lenghtRCU=meanlenght(Cred_Reg_U)
  
  # Coverage probability
  covprob=function(Interval){
    indicator=function(i){
      if(Interval[i,1]<=varphi & Interval[i,2]>=varphi){
        return(1)
      }else{return(0)}
    }
    resind=NA
    for(i in 1:nrow(Interval)){
      resind[i]=indicator(i)
    }
    covint=sum(resind)/nrow(Interval)
    
    return(covint=covint)
  }
  covprobILM=covprob(IntervalML)
  covprobIAMLE=covprob(IntervalIAMLE)
  covprobITau=covprob(IntervalTau)
  covprobISpe=covprob(IntervalSpe)
  covprobRCT=covprob(Cred_Reg_T)
  covprobRCB=covprob(Cred_Reg_B)
  covprobRCU=covprob(Cred_Reg_U)
  
  Descritive=matrix(data=c(c(DesML,Bias[1],lenghtIAMLE,covprobIAMLE),
                           c(DesML,Bias[1],lenghtILM,covprobILM),
                           c(DesMTau,Bias[2],lenghtITau,covprobITau),
                           c(DesMSpe,Bias[3],lenghtISpe,covprobISpe),
                           c(DesBT,Bias[4],lenghtRCT,covprobRCT),
                           c(DesBB,Bias[5],lenghtRCB,covprobRCB),
                           c(DesBU,Bias[6],lenghtRCU,covprobRCU),
                           c(min(Error_Sample_Boot),mean(Error_Sample_Boot),
                             sd(Error_Sample_Boot),max(Error_Sample_Boot),0,0,0)),nrow=7,ncol=8,
                    dimnames = list(list("Min","Mean","SD","Max","Bias","Length","Coverage"),
                                    list("ML_Asym","ML_Boot","Tau_Boot","Spe_Boot","Triangular","Beta","Unif","Error_S_B")))
  return(Descritive)
}

#Example
#Sim_Est_N_Samp(0.25,n=10,nboot=200, N=5, confidence=0.95,a=1,b=1,c=0,n_burnd=1000,n_iter=5000,n_thin=10,n_chain=1)

####################################################################
# Obtaining hyperparameter values, Tovar's method.
####################################################################
# x1, x2 are quantiles obtained from an elicited process 
# 1-alph is the confidence that the interval (x1,x2) contains 
# the true value of the parameter
Mtovar=function(x1,x2,alp){
  tht0=(x1+x2)/2
  w=(tht0+1)/(1-tht0)
  sig=sqrt(alp)*(x1-tht0)
  b= (4*w-((w+1)*sig)^2)/((w+1)^3*sig^2)
  a=w*b
  return(list(a=a,b=b,c=tht0))
}


################################################################################
## Graph construction
################################################################################
# m: descriptive measure to graph
# omega: value of dependence parameter varphi.
# omega_to_file: list that contains the value of varphi and file name
# in format .xlsx. For example, omega_to_file=list("-0.75"="Results-075.xlsx",
# "-0.25"="Results-025.xlsx")
# useDependence: If this is True, the graphs of the method cho_method=7 
# are constructed, in this case of the a priori Beta. therwise, the graphs 
# of all methods are constructed.
# m: descriptive measure that will be graphed. 
GeneralGraph = function(m, omega = NULL, omega_to_file = NULL, useDependence = FALSE,useAsymtotic=FALSE,cho_method=7) {
  nameGra = c("Min", "Mean", "SD", "Max", "Bias", "Length", "Coverage")
  colName = ifelse(useDependence, "Dependence", "Method")
  if(useDependence) {labelSet =c("SND","WND", "WPD", "SPD")} else if(useAsymtotic) {
    labelSet = c("ML_Asymptotic","ML", "Tau-Kendall", "Spearman", "Triangular", "Beta", "Uniform")} else
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
      columnIndex = if(useDependence) cho_method else if(useAsymtotic) 2:8 else 3:8
      
      for(j in columnIndex) {
        MA = rbind(MA, c(size[i], A[m, j], if(useDependence) w else if(useAsymtotic) namemethod[(j - 1)] else namemethod[(j - 2)]))
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
# empirical marginal distribution
######################################################################
Marginal=function(c,X){
  step.my=function(w,z){
    if(w<=z){
      1
    }else {0}
  }
  N=nrow(X); fn=NULL
  for(i in 1:nrow(X)){
    d=0
    for(j in 1:nrow(X)){
      e=step.my(X[j,c],X[i,c])
      d=e+d
    } 
    fn[i]=1/(N+1) * d
  }
  fn
}
