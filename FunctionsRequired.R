###################################################################################################
#                                                                                                 #
# Description:                                                                                    #
# This R script contains a comprehensive set of functions for simulating and estimating           #
# the dependence parameter (φ) of the Farlie-Gumbel-Morgenstern (FGM) copula. It includes         #
# routines for generating random samples, computing the likelihood and its derivatives,           #
# and performing parameter estimation using classical (MLE, moments, Blomqvist's beta)            #
# and Bayesian methods.                                                                           #
#                                                                                                 #
# The script implements the Metropolis-Hastings algorithm with various prior distributions        #
# (uniform, beta, triangular) and instrumental proposals. Diagnostic tools for monitoring         #
# MCMC convergence are provided, including trace plots, histograms, ACF plots, and                #
# cumulative averages. Simulation studies are also included to assess bias, MSE, coverage,        #
# and average length of the credible/confidence intervals across different scenarios.             #
#                                                                                                 #
# Author: Llerzy Torres Ome                                                                       #
# Creation Date: july 1, 2025                                                                     #
#                                                                                                 #
# Main functionalities:                                                                           #
# 1. Density_FGM: FGM copula density function.                                                    #
# 2. Der_log_lik_FGM: Derivatives of the log-likelihood function.                                 #
# 3. rfgm / rnfgm: Generation of FGM random samples (1 or n observations).                        #
# 4. me: Moment-based estimators (Pearson, Kendall, Spearman).                                    #
# 5. metropolis_hastings: MCMC sampling with flexible priors/instrumentals.                       #
# 6. Graphs: MCMC convergence diagnostics.                                                        #
# 7. Est_A_Samp: Estimation from a single sample (point, interval, credibility).                  #
# 8. Sim_Est_A_Samp: Parameter estimation with filtering of invalid samples.                      #
# 9. Sim_Est_N_Samp: Repeated simulations for performance evaluation.                             #
# 10. Mtovar: Prior elicitation for the beta distribution (Tovar’s method).                       #
# 11. GeneralGraph: Aggregated plot generation for simulation results.                            #
# 12. Marginal: Empirical marginal distribution estimation.                                       #
#                                                                                                 #
###################################################################################################

###################################################################################################
# Library
###################################################################################################
library(copula)
library(xlsx)
library(readxl)
library(ggplot2)
library(gridExtra)
library(correlation)
library(betafunctions)
library(parallel)
####################################################################
# Density function of FGM copula.
Density_FGM = function(u,v,varphi){
  return(1+varphi*(1-2*u)*(1-2*v))
}
# Likelihood function
likelihood = function(u, v, varphi) {
  return(prod(Density_FGM(u, v, varphi)))
}

# Derivative of the log-likelihood function.
# der: order of derivative.
Der_log_lik_FGM = function(varphi,u,v,der = 1){
  dln = NA;  n = length(u)
  for (i in 1:n) {
    num = (1-2*u[i])*(1-2*v[i])
    den = 1+varphi*num
    if(der == 1){
      dln[i] = (num/den)
    }else if(der==2){
      dln[i] =-(num/den)^2
    }
  }
  return(sum(dln))
}

###################################################################################################
## Data simulation
###################################################################################################

# Generation of one random observation of FGM copula
rfgm = function(varphi){
  u1 = runif(n = 1,0,1) 
  v2 = runif(n = 1,0,1) # conditional probability
  A = varphi*(2*u1-1)-1
  B = sqrt(1-2*varphi*(2*u1-1)+varphi^2*(2*u1-1)^2+4*varphi*v2*(2*u1-1))
  u2 = 2*v2/(B-A)
  return(c(u1,u2))
}

## Generation of n random observation of FGM copula
rnfgm = function(varphi,n){
  datafgm = matrix(data = NA,nrow = 0,ncol = 2,dimnames = list(c(),c("u","v")))
  for (i in 1:n) {
    datafgm = rbind(datafgm,rfgm(varphi))
  }
  return(datafgm)
}
# Example
# dim(rnfgm(0.25,10))[1]

###################################################################################################
## Bayesian Inference
###################################################################################################

# Density of the Beta prior distribution
prior_density_beta = function(phi, alpha, beta) {
  return((phi + 1)^(alpha - 1) * (1 - phi)^(beta - 1))
}

# Density of the Uniform prior distribution
prior_density_uniform = function(phi) {
  return(1 / 2)  
}

# Density of the triangular prior distribution
prior_density_triangular = function(phi, c) {
  if (phi >= -1 && phi < c) {
    return((phi + 1) / (c + 1))
  } else if (phi >= c && phi <= 1) {
    return((1 - phi) / (1 - c))
  } else {
    return(0)
  }
}

###################################################################################################
# Metropolis-Hastings Algorithm
###################################################################################################
# u and v are the sample observations.
# alpha_prior, beta_prior, and c_prior are the values of the hyperparameters.
# instrumental: determines the instrumental distribution that can be used, it can be "beta" or "uniform".
# precision: works when the instrumental distribution is the beta distribution.
# n_iter: length of the generated chain.
# init_varphi: seed with which the chain starts.
metropolis_hastings = function(u, v, prior, alpha_prior = NA, beta_prior = NA, c_prior = NA, 
                                instrumental = "beta",precision=20,
                                n_iter, init_varphi = runif(1,-1,1), burn_in = 0 , thinning = 1) {
  
  total_iter = burn_in + n_iter * thinning
  varphi <- numeric(total_iter)
  varphi[1] <- init_varphi
  acceptance_samples <- 0
  
  for (i in 2:total_iter) {
    if(instrumental =="beta"){
      
      alpha_current = (varphi[i - 1] + 1) / (2) * precision
      beta_current = (1 - varphi[i - 1]) / (2) * precision
      
      varphi_proposal = rBeta.4P(n = 1, l = -1, u = 1, alpha = alpha_current, beta = beta_current)
      
      alpha_inst_proposal = (varphi_proposal + 1) / (2) * precision
      beta_inst_proposal = (1 - varphi_proposal) / (2) * precision
      
      # Density of the beta instrumental distribution with current and proposed value of varphi.
      proposal_density_current_given_proposal = dBeta.4P(x = varphi[i-1],l = -1,u = 1, alpha_inst_proposal, beta_inst_proposal)
      proposal_density_proposal_given_current = dBeta.4P(x = varphi_proposal,l = -1,u = 1, alpha_current, beta_current)
    } else if (instrumental == "uniform"){
      varphi_proposal = runif(n = 1,-1,1)
      proposal_density_current_given_proposal = 1/2
      proposal_density_proposal_given_current = 1/2
    }
    
    # Selection of the a prior distribution
    if (prior == "beta") {
      prior_density_current = prior_density_beta(varphi[i-1], alpha_prior, beta_prior)
      prior_density_proposal = prior_density_beta(varphi_proposal, alpha_prior, beta_prior)
    } else if (prior == "uniform") {
      prior_density_current = prior_density_uniform(varphi[i-1])
      prior_density_proposal = prior_density_uniform(varphi_proposal)
    } else if (prior == "triangular") {
      prior_density_current = prior_density_triangular(varphi[i-1], c_prior)
      prior_density_proposal = prior_density_triangular(varphi_proposal, c_prior)
    }
    
    # Calculate the acceptance ratio
    acceptance_ratio = (likelihood(u, v, varphi_proposal) * prior_density_proposal) / 
      (likelihood(u, v, varphi[i-1]) * prior_density_current)
    
    # Adjust the acceptance ratio with the instrumental density
    acceptance_ratio = acceptance_ratio * proposal_density_current_given_proposal / proposal_density_proposal_given_current
    
    # Identification of errors
    if (is.nan(acceptance_ratio) || is.infinite(acceptance_ratio)) {
      stop(paste("Non-numeric acceptance_ratio detected at iteration",i, 
                 "with proposal", varphi_proposal, 
                 "and previous chain value", varphi[i - 1],"the acceptance ration is", acceptance_ratio, 
                 "the density current value is",proposal_density_current_given_proposal,
                 "the density proposal value is", proposal_density_proposal_given_current,
                 "the value of the parameter proposal and current are", alpha_inst_proposal, beta_inst_proposal,
                 alpha_current, beta_current))
    }
    if (runif(1) < acceptance_ratio) {
      varphi[i] = varphi_proposal
      acceptance_samples = acceptance_samples + 1
    } else {
      varphi[i] = varphi[i-1]
    }
  }
  
  acceptance_rate <- acceptance_samples / (total_iter - 1)
  return(list(varphi = varphi, acceptance_rate = acceptance_rate))
}

###################################################################################################
# Function that plots the histogram, density, trace, and convergence control using the average.
###################################################################################################
# "nameaxisy" is the name of the vertical axis for trace and convergence monitoring.
# "width" is the width for the confidence intervals of convergence monitoring.
# "lscatt" is an increment to the minimum value generated. It allows plotting the line at an "lscatt" distance from the trace to enhance visualization.
# "uscatt" is an increment to the maximum value generated. It allows plotting the line at an "uscatt" distance from the trace to enhance visualization.
Graphs = function(dataset, nameaxisy, width = 10, lscatt = 0.05, uscatt = 0.05) {
  # Function to apply a consistent theme
  custom_theme = theme_minimal(base_size = 10) +
    theme(
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  # Histogram with density
  l = length(dataset[, 1])
  hist = ggplot(dataset, aes(x = dataset[, 1])) + 
    geom_histogram(aes(y = after_stat(density)), colour = 1, fill = "white", binwidth = 0.02) +
    geom_density(lwd = 1.2, linetype = 2, colour = 2, fill = 4, alpha = 0.25) +
    labs(title = "(a) Histogram and Density") + ylab("Density") +
    xlab(expression("Chain Values of" ~ phi^(t))) +
    custom_theme
  
  # Trace plot with maximum and minimum
  trace = ggplot(dataset, aes(x = 1:l, y = dataset[, 1])) +
    geom_line() + xlab("Chain Iterations (t)") +
    ylab(expression("Chain of" ~ phi^(t))) +
    ylim(c(min(dataset) - lscatt, max(dataset) + uscatt)) +
    geom_hline(aes(yintercept = min(dataset[, 1])), colour = "red", linetype = 2) +
    annotate("text", x = l - l/10, y = min(dataset[, 1]), label = round(min(dataset[, 1]), 3), vjust = 2, colour = "red") +
    geom_hline(aes(yintercept = max(dataset[, 1])), colour = "red", linetype = 2) +
    annotate("text", x = l - l/10, y = max(dataset[, 1]), label = round(max(dataset[, 1]), 3), vjust = -1, colour = "red") +
    labs(title = substitute(list("(b) Trace of the random sample of size", n), list(n = l))) +
    custom_theme
  
  # Acf plot
  alfa = 0.05
  lim = qnorm((1 - alfa / 2)) / sqrt(l)
  acf_values = acf(dataset, plot = FALSE)
  acf_data = data.frame(Lag = acf_values$lag[-1],  # Remove the first lag value (always 0)
                        ACF = acf_values$acf[-1])  # Remove the first ACF value (always 1)
  acfplot = ggplot(acf_data, aes(x = Lag, y = ACF)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = c(lim, -lim), linetype = "dashed") +
    labs(title = "(d) Autocorrelation Function",
         x = "Lag",
         y = "ACF") +
    custom_theme  
  
  # Convergence control using averaging
  dataset$estintden = cumsum(dataset[, 1]) / (1:l)
  dataset$esterrden = sqrt(cumsum((dataset[, 1] - dataset$estintden)^2)) / (1:l)
  
  mean_X1_X2 = ggplot(dataset, aes(x = 1:l, y = estintden)) + geom_line() +
    geom_line(aes(x = 1:l, y = estintden - 1.95 * esterrden, colour = "Upper")) +
    geom_line(aes(x = 1:l, y = estintden + 1.95 * esterrden, colour = "Lower")) +
    ylim(mean(dataset$estintden) + width * c(-dataset$esterrden[l], dataset$esterrden[l])) +
    ylab(expression("Accumulated Average of" ~ phi^(t))) +
    xlab("Chain Iterations (t)") + geom_hline(yintercept = mean(dataset[, 1]), colour = "red", linetype = 2) +
    annotate("text", x = l - l/10, y = mean(dataset[, 1]), label = round(mean(dataset[, 1]), 3), vjust = -2, colour = "red") +
    labs(title = "(c) Convergence Control using Averaging", color = "Bounds") + 
    scale_shape_discrete(name = " ") + custom_theme
  
  # Plots of histogram, trace, convergence control, and acf.
  grid.arrange(hist, trace, mean_X1_X2, acfplot, 
               ncol = 2, nrow = 2, widths = c(4, 4), heights = c(2, 2), layout_matrix = rbind(c(1, 2), c(3, 4)))
}


###################################################################################################
# Moments estimate
###################################################################################################

me = function(dataset_target){
  pearson = cor(dataset_target,method = "pearson")[1,2]*3
  kendall = cor(dataset_target,method = "kendall")[1,2]*9/2
  spearman = cor(dataset_target,method = "spearman")[1,2]*3
  return(list(varphipearson = pearson,varphikendall = kendall,varphispearman = spearman))
}

###################################################################################################
# Estimates (point estimates, interval estimates, and credibility regions) 
# produced from a random sample.
###################################################################################################
# varphi: Dependence parameter.
# n: Sample size.
# nboot: Number of resamples performed in the bootstrap procedure.
# confidence: Desired confidence level for interval estimates; the probability 
# that a new sample will produce estimates within this interval.
# alpha_prior and beta_prior: Hyperparameters of the prior Beta distribution.
# c_prior: Hyperparameter of the prior Triangular distribution.
# thinning: Thinning interval for the MCMC chain (every nth iteration is kept).

# Parameter estimation using the dataset_target sample.
Est_A_Samp = function(dataset_target,nboot = 100,confidence = 0.95,alpha_prior, beta_prior, c_prior,
                    burn_in = 1000, n_iter = 5000, thinning = 1,type_sample = "all", precision = 20){
  total_iter = burn_in + n_iter * thinning
  ##############################################
  # Point estimates for the FGM dependence parameter
  ##############################################
  n = dim(dataset_target)[1]
  # Estimation of maximum log-likelihood
  MLE = safeUroot(f=function(ph){mapply(function(ph){
    Der_log_lik_FGM(ph,dataset_target[,1],dataset_target[,2])},ph)},c(-1,1))$root
  # Moment estimate
  EM = me(dataset_target)
  
  # Inversion of Blomqvist’s beta
  
  Blomqvist_Beta = 4*cor_test(data = as.data.frame(dataset_target), x="u",y="v", method = "blomqvist")$r
  
  # Estimates using Bayes approach.
  ResTriangular_all = metropolis_hastings(u=dataset_target[,1], v=dataset_target[,2], 
                                    prior="triangular", c_prior = c_prior, instrumental="uniform",
                                    precision=20, n_iter=n_iter, burn_in = burn_in, thinning = thinning)
  ResBeta_all = metropolis_hastings(u=dataset_target[,1], v=dataset_target[,2], 
                              prior="beta", alpha_prior, beta_prior, 
                              instrumental="uniform",precision=20,n_iter=n_iter, 
                              burn_in = burn_in , thinning = thinning)
  ResBetaBeta_all = metropolis_hastings(u=dataset_target[,1], v=dataset_target[,2], 
                                    prior="beta", alpha_prior, beta_prior, 
                                    instrumental="beta",precision = precision,n_iter=n_iter, 
                                    burn_in = burn_in , thinning = thinning)
  ResUnif_all = metropolis_hastings(u=dataset_target[,1], v=dataset_target[,2], 
                              prior="uniform",instrumental="uniform",precision=20,n_iter=n_iter, 
                              burn_in = burn_in, thinning = thinning)
  
  ResTriangular = ResTriangular_all$varphi[seq((burn_in+1), total_iter, by=thinning)]
  ResBeta = ResBeta_all$varphi[seq((burn_in+1), total_iter, by=thinning)]
  ResBetaBeta = ResBetaBeta_all$varphi[seq((burn_in+1), total_iter, by=thinning)]
  ResUnif = ResUnif_all$varphi[seq((burn_in+1), total_iter, by=thinning)]
  ##############################################
  # Interval estimation
  ##############################################
  
  alpha = 1-confidence
  # Bootstrap interval
  vphitau = c();vphispe=c();vphilm=c();vBlomqvist_Beta=c()  #Bootstrap estimates for the dependence parameter using classical estimators.
  i = 1 # bootstrap counter.
  
  if(type_sample == "all"){
    while(i <= nboot) {
      cho = sample(1:n,n,replace = T)     # Resampling the position in the sample.
      rdatan = dataset_target[cho,]       # Selecting elements from the sample.
      #Estimates
      philm = safeUroot(f=function(ph){mapply(function(ph){
        Der_log_lik_FGM(ph,rdatan[,1],rdatan[,2])},ph)},c(-1,1))$root
      
      vphilm[i] = philm
      vphim = me(rdatan)
      vphitau[i] = vphim$varphikendall
      vphispe[i] = vphim$varphispearman
      vBlomqvist_Beta[i] = 4*cor_test(data = as.data.frame(rdatan), x="u",y="v", method = "blomqvist")$r
      i = i+1
    }} else if (type_sample == "exclude"){
      ie = 0 # counter of the number of times bootstrap resampling produces estimates outside [-1,1]. 
      while(i <= nboot) {
        cho = sample(1:n,n,replace = T)     # Resampling the position in the sample.
        rdatan = dataset_target[cho,]       # Selecting elements from the sample.
        #Estimates
        philm = safeUroot(f = function(ph){mapply(function(ph){
          Der_log_lik_FGM(ph,rdatan[,1],rdatan[,2])},ph)},c(-1,1))$root
        if(all(abs(philm) <= 1)){
          vphilm[i] = philm
          vphim = me(rdatan)
          SvBlomqvist_Beta = 4*cor_test(data = as.data.frame(rdatan), x="u",y="v", method = "blomqvist")$r
          if(all(abs(vphim$varphikendall) <= 1,abs(vphim$varphispearman)<=1,abs(SvBlomqvist_Beta)<=1)){
            vphitau = cbind(vphitau,vphim$varphikendall)
            vphispe = cbind(vphispe,vphim$varphispearman)
            vBlomqvist_Beta[i] = 4*cor_test(data = as.data.frame(rdatan), x="u",y="v", method = "blomqvist")$r
            i = i+1
          }}else{ie = ie+1}
      }}
  
  IntervalML = quantile(vphilm,probs=c(alpha/2,1-alpha/2))   # Quantile for estimates using Maximum likelihood.
  IntervalTau = quantile(vphitau,probs=c(alpha/2,1-alpha/2)) # Quantile for estimates using moment method.
  IntervalSpe = quantile(vphispe,probs=c(alpha/2,1-alpha/2)) # Quantile for estimates using moment method.
  IntervalBlomqvist_Beta = quantile(vBlomqvist_Beta,probs=c(alpha/2,1-alpha/2)) # Quantile for estimates using moment method.
  
  # Asymptotic (1-alpha)100% confidence interval for the MLE
  variance_MLE = 1/Der_log_lik_FGM(varphi = MLE,u = dataset_target[,1],v = dataset_target[,2],der = 2)
  Interval_A_MLE = MLE + c(qnorm(alpha/2),qnorm(1-alpha/2))*sqrt(abs(variance_MLE))
  
  # Asymptotic (1-alpha)100% confidence interval for the Tau Kendall
  # Interval_A_Tau = EM$varphikendall + 9/2*c(qnorm(alpha/2),qnorm(1-alpha/2))*sqrt(2*(2*n+5)/(9*n*(n-1)))
  Interval_A_Tau = EM$varphikendall + c(qnorm(alpha/2),qnorm(1-alpha/2))*
    sqrt( 4/9-(184/2025)* EM$varphikendall^2) / ((2/9)*sqrt(n))
  
  # Asymptotic (1-alpha)100% confidence interval for the Spearman 
  # Interval_A_Spe = EM$varphispearman + 3*c(qnorm(alpha/2),qnorm(1-alpha/2))*sqrt(1/(n-1))
  Interval_A_Spe = EM$varphispearman + c(qnorm(alpha/2),qnorm(1-alpha/2))*
    sqrt(1-(11/45)*(EM$varphispearman)^2) / ((1/3)*sqrt(n))
  
  # Asymptotic (1-alpha)100% confidence interval for the Blomqvist_Beta 
  Interval_A_Blomqvist_Beta = Blomqvist_Beta + c(qnorm(alpha/2),qnorm(1-alpha/2))*
    sqrt(1-(Blomqvist_Beta)^2/16) / ((1/4)*sqrt(n))
  
  # Credibility regions
  Cred_Reg_T = quantile(ResTriangular,probs=c(alpha/2,1-alpha/2))
  Cred_Reg_B = quantile(ResBeta,probs=c(alpha/2,1-alpha/2))
  Cred_Reg_BB = quantile(ResBetaBeta,probs=c(alpha/2,1-alpha/2))
  Cred_Reg_U = quantile(ResUnif,probs=c(alpha/2,1-alpha/2))
  
  Summary = matrix(data = c(MLE,IntervalML,Interval_A_MLE, NA,NA,
                          EM$varphikendall,IntervalTau,Interval_A_Tau, NA,NA,
                          EM$varphispearman,IntervalSpe,Interval_A_Spe, NA,NA,
                          Blomqvist_Beta,IntervalBlomqvist_Beta,Interval_A_Blomqvist_Beta,NA,NA,
                          mean(ResTriangular),Cred_Reg_T,NA,NA, ResTriangular_all$acceptance_rate,1,
                          mean(ResBeta),Cred_Reg_B,NA,NA, ResBeta_all$acceptance_rate,1,
                          mean(ResBetaBeta),Cred_Reg_BB,NA,NA, ResBetaBeta_all$acceptance_rate,2,
                          mean(ResUnif),Cred_Reg_U,NA,NA, ResUnif_all$acceptance_rate, 1),nrow =8 ,ncol =7,byrow = TRUE,
                 dimnames = list(list("LM","Tau","Spe","Blomqvist Beta","Triangular","Beta","BetaBeta","Unif"),
                                 list("Mean","Low","Upper","AL","AU","Acceptance","Instrumental")))
  
  return(list(MLE=MLE,EM=EM,Blomqvist_Beta = Blomqvist_Beta,
              EBT=mean(ResTriangular),EBB=mean(ResBeta),EBBB=mean(ResBetaBeta),EBU=mean(ResUnif),
              IML=IntervalML,ITau=IntervalTau,ISpe=IntervalSpe,IBB=IntervalBlomqvist_Beta,
              CRT=Cred_Reg_T,IAMLE=Interval_A_MLE,
              IATAU=Interval_A_Tau,IASPE=Interval_A_Spe,IABB=Interval_A_Blomqvist_Beta,
              CRB=Cred_Reg_B,CRBB=Cred_Reg_BB, CRU=Cred_Reg_U, vphilm=vphilm,
              vphitau=vphitau,vphispe=vphispe,vphibb=vBlomqvist_Beta,Sims_T=ResTriangular,
              Sims_B=ResBeta,Sims_BB=ResBetaBeta,Sims_U=ResUnif,n_Chain_Length=length(ResTriangular),
              Summary=Summary))
}

#Example
#trial0=rnfgm(-0.75,600)
#trial1=Est_A_Samp(trial0,alpha_prior = 12,beta_prior = 14,c_prior = -0.75,type_sample = "exclude")
#trial1$Summary[5:8,6]

# Parameter estimation using a randomly generated sample of size n and varphi dependence.
Sim_Est_A_Samp = function(varphi,n,nboot = 100,confidence = 0.95,alpha_prior, beta_prior, c_prior,
                        burn_in = 1000,n_iter = 5000,thinning = 1,precision = 20){
  dataset_target=rnfgm(varphi,n) # random sample
  
  # Estimation of maximum log-likelihood
  MLE = safeUroot(f=function(ph){mapply(function(ph){
    Der_log_lik_FGM(ph,dataset_target[,1],dataset_target[,2])},ph)},c(-1,1))$root
  # Moment estimate
  EM = me(dataset_target)
  # Inversion of Blomqvist’s beta
  Blomqvist_Beta = 4*cor_test(data = as.data.frame(dataset_target), x="u",y="v", method = "blomqvist")$r
  # Selecting samples that have Maximum Likelihood and Moments estimates 
  # within the domain of the dependency parameter of the FGM copula.
  while(any(abs(MLE) > 1, abs(EM$varphikendall) > 1, abs(EM$varphispearman) > 1, abs(Blomqvist_Beta)>1)){
    dataset_target=rnfgm(varphi,n)
    MLE=safeUroot(f = function(ph){mapply(function(ph){
      Der_log_lik_FGM(ph,dataset_target[,1],dataset_target[,2])},ph)},c(-1,1))$root
    EM=me(dataset_target)
    Blomqvist_Beta = 4*cor_test(data = as.data.frame(dataset_target), x="u",y="v", method = "blomqvist")$r
  }
  results_Sim = Est_A_Samp(dataset_target,nboot,confidence,alpha_prior, beta_prior, c_prior,
                        burn_in,n_iter,thinning,type_sample="exclude", precision=precision)
  return(results_Sim)
}

#Example
#tiral1=Sim_Est_A_Samp(-0.75,100,nboot=100,confidence = 0.95,alpha_prior=12, beta_prior=14, c_prior=-0.75,
#                       burn_in=1000,n_iter=5000,thinning=2)

###################################################################################################
# Simulation study for evaluating the characteristics (mean squared error, 
# minimum, maximum, average length, coverage probability) of
# estimators from MLE, MM, Blomqvist_Beta, and Bayesian approaches (Triangular, 
# Beta, Uniform).
###################################################################################################

# varphi: Dependence parameter.
# n: Sample size.
# N: Number of FGM samples generated.
# nboot: Number of resamples performed in the bootstrap procedure.
# confidence: Desired confidence level for interval estimates; the probability 
# that a new sample will produce estimates within this interval.
# alpha_prior and beta_prior: Hyperparameters of the prior Beta distribution.
# c_prior: Hyperparameter of the prior Triangular distribution.
# burn_in: Number of iterations discarded at the beginning of the MCMC chain (burn-in period).
# n_iter: Total number of iterations in the MCMC chain.
# thinning: Thinning interval for the MCMC chain (every nth iteration is kept).

Sim_Est_N_Samp = function(varphi,n,nboot = 100, N, confidence,alpha_prior,beta_prior,c_prior,burn_in = 1000,
                        n_iter = 5000,thinning = 1, precision = 20){
  PhiML= NA; PhiTau = NA; PhiSpe = NA;PhiB_B = NA; PhiBT = NA; PhiBB = NA; PhiBBB = NA; PhiBU = NA; Error_Sample_Boot = NA
  IntervalML = matrix(data=NA,nrow=0,ncol=2,dimnames = list(c(),c("L","R")))
  IntervalTau = IntervalML; IntervalSpe=IntervalML; IntervalBB=IntervalML
  IntervalIAMLE = IntervalML;IntervalIATAU=IntervalML;IntervalIASPE=IntervalML;
  IntervalABB = IntervalML;Cred_Reg_T=IntervalML;Cred_Reg_B=IntervalML;Cred_Reg_BB=IntervalML;Cred_Reg_U=IntervalML
  acceptance = matrix(data=NA,nrow=0,ncol=4,dimnames = list(c(),c("Triangular","Beta","BetaBeta","Unif")))
  i = 1
  while (i<=N) {
    results = Sim_Est_A_Samp(varphi,n,nboot,confidence,alpha_prior,beta_prior,c_prior,
                           burn_in,n_iter,thinning, precision = precision)
    # Estimates
    PhiML[i] = results$MLE; PhiTau[i]=results$EM$varphikendall; PhiSpe[i]=results$EM$varphispearman
    PhiB_B[i] = results$Blomqvist_Beta;
    PhiBT[i] = results$EBT; PhiBB[i]=results$EBB; PhiBBB[i]=results$EBBB; PhiBU[i]=results$EBU
    # Intervals estimation
    IntervalML = rbind(IntervalML,results$IML)
    IntervalIAMLE = rbind(IntervalIAMLE,results$IAMLE)
    IntervalTau = rbind(IntervalTau,results$ITau)
    IntervalIATAU = rbind(IntervalIATAU,results$IATAU)
    IntervalSpe = rbind(IntervalSpe,results$ISpe)
    IntervalIASPE = rbind(IntervalIASPE,results$IASPE)
    IntervalBB = rbind(IntervalBB,results$IBB)
    IntervalABB = rbind(IntervalABB,results$IABB)
    
    Cred_Reg_T = rbind(Cred_Reg_T,results$CRT)
    Cred_Reg_B = rbind(Cred_Reg_B,results$CRB)
    Cred_Reg_BB = rbind(Cred_Reg_B,results$CRBB)
    Cred_Reg_U = rbind(Cred_Reg_U,results$CRU)
    acceptance = rbind(acceptance, results$Summary[5:8,6])
    #Error_Sample_Boot[i]=results$Er_S_B
    i = i+1}
  # Descriptive measures
  DesML = c(min(PhiML),mean(PhiML),sd(PhiML),max(PhiML))
  DesMTau =c(min(PhiTau),mean(PhiTau),sd(PhiTau),max(PhiTau))
  DesMSpe =c(min(PhiSpe),mean(PhiSpe),sd(PhiSpe),max(PhiSpe))
  DesMBB = c(min(PhiB_B),mean(PhiB_B),sd(PhiB_B),max(PhiB_B))
  DesBT = c(min(PhiBT),mean(PhiBT),sd(PhiBT),max(PhiBT))
  DesBB = c(min(PhiBB),mean(PhiBB),sd(PhiBB),max(PhiBB))
  DesBBB = c(min(PhiBBB),mean(PhiBBB),sd(PhiBBB),max(PhiBBB))
  DesBU = c(min(PhiBU),mean(PhiBU),sd(PhiBU),max(PhiBU))
  # Estimator bias
  Bias = c(varphi-DesML[2],varphi-DesMTau[2],varphi-DesMSpe[2],varphi-DesMBB[2],
         varphi-DesBT[2],varphi-DesBB[2],varphi-DesBBB[2],varphi-DesBU[2])
  
  # Average length 
  meanlenght = function(Interval){
    lenint = sum(Interval[,2]-Interval[,1])/nrow(Interval)
    return(lenint = lenint)
  }
  
  lenghtILM = meanlenght(IntervalML)
  lenghtIAMLE = meanlenght(IntervalIAMLE)
  lenghtITau = meanlenght(IntervalTau)
  lenghtIATAU = meanlenght(IntervalIATAU)
  lenghtISpe = meanlenght(IntervalSpe)
  lenghtIASPE = meanlenght(IntervalIASPE)
  lenghtIBB = meanlenght(IntervalBB)
  lenghtIABB = meanlenght(IntervalABB)
  
  lenghtRCT = meanlenght(Cred_Reg_T)
  lenghtRCB = meanlenght(Cred_Reg_B)
  lenghtRCBB = meanlenght(Cred_Reg_BB)
  lenghtRCU = meanlenght(Cred_Reg_U)
  
  # Coverage probability
  covprob = function(Interval){
    indicator = function(j){
      if(Interval[j,1]<=varphi & Interval[j,2]>=varphi){
        return(1)
      }else{return(0)}
    }
    resind = NA
    for(l in 1:nrow(Interval)){
      resind[l] = indicator(l)
    }
    covint = sum(resind)/nrow(Interval)
    
    return(covint = covint)
  }
  covprobILM = covprob(IntervalML)
  covprobIAMLE = covprob(IntervalIAMLE)
  covprobITau = covprob(IntervalTau)
  covprobIATAU = covprob(IntervalIATAU)
  covprobISpe = covprob(IntervalSpe)
  covprobIASPE = covprob(IntervalIASPE)
  covprobIBB = covprob(IntervalBB)
  covprobIABB = covprob(IntervalABB)
  
  covprobRCT = covprob(Cred_Reg_T)
  covprobRCB = covprob(Cred_Reg_B)
  covprobRCBB = covprob(Cred_Reg_BB)
  covprobRCU = covprob(Cred_Reg_U)
  
  Descritive = matrix(data=c(c(DesML,Bias[1],var(PhiML)+(Bias[1])^2,lenghtILM,covprobILM,lenghtIAMLE,covprobIAMLE,NA),
                           c(DesMTau,Bias[2],var(PhiTau)+(Bias[2])^2,lenghtITau,covprobITau,lenghtIATAU,covprobIATAU,NA),
                           c(DesMSpe,Bias[3],var(PhiSpe)+(Bias[3])^2,lenghtISpe,covprobISpe,lenghtIASPE,covprobIASPE,NA),
                           c(DesMBB,Bias[4],var(PhiB_B)+(Bias[4])^2,lenghtIBB,covprobIBB,lenghtIABB,covprobIABB,NA),
                           c(DesBT,Bias[5],var(PhiBT)+(Bias[5])^2,lenghtRCT,covprobRCT,NA,NA,mean(acceptance[,1])),
                           c(DesBB,Bias[6],var(PhiBB)+(Bias[6])^2,lenghtRCB,covprobRCB,NA,NA,mean(acceptance[,2])),
                           c(DesBBB,Bias[7],var(PhiBBB)+(Bias[7])^2,lenghtRCBB,covprobRCBB,NA,NA,mean(acceptance[,3])),
                           c(DesBU,Bias[8],var(PhiBU)+(Bias[8])^2,lenghtRCU,covprobRCU,NA,NA,mean(acceptance[,4]))),
                    nrow=11,ncol=8,
                    dimnames = list(list("Min","Mean","SD","Max","Bias","MSE","BootLength","BootCoverage","ALength","ACoverage","Acceptance"),
                                    list("ML","Tau","Spe","Blomqvist_Beta","Triangular","Beta","BetaBeta","Unif")))
  return(Descritive)
}

#Example
#trial1 = Sim_Est_N_Samp(varphi=-0.75,n=10,nboot=100, N=10, confidence=0.95,alpha_prior=2.0625,beta_prior=14.4375,
#               c_prior=-0.75,burn_in=5000,n_iter=5000,thinning=5,precision = 20)

#round(trial1,4)

###################################################################################################
# Obtaining hyperparameter values, modified Tovar's method
###################################################################################################
# x1, x2: Quantiles obtained from an elicited process.
# 1-alpha: Confidence level that the interval (x1, x2) contains the true value of the parameter.
# low: Lower limit of the Beta distribution.
# upp: Upper limit of the Beta distribution.
Mtovar = function(x1,x2,alp,low=-1,upp=1){
  tht0  = (x1+x2)/2
  w     = (tht0-low)/(upp-tht0)
  sig   = sqrt(alp)*(x1-tht0)
  b     = ((upp-low)^2*w-((w+1)^2*sig^2))/((w+1)^3*sig^2)
  a     = w*b
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
# cho_method: The method to use for constructing the graphs. Default is 7.
GeneralGraph = function(m, omega = NULL, omega_to_file = NULL, useDependence = FALSE,useAsymtotic = FALSE,cho_method = 7 ) {
  nameGra = c("Min", "Mean", "SD", "Max", "Bias","MSE", "Length", "Coverage")
  colName = ifelse(useDependence, "Dependence", "Method")
  
  if(useDependence) 
    {labelSet = c("SND","WND", "WPD", "SPD")} 
  else if(useAsymtotic) {
    labelSet = c("ML", "Tau-Kendall", "Spearman", "Blomqvist Beta", "Triangular-Unif", "Beta-Unif", "Beta-Beta", "Unif-Unif","ML A","Tau A","Spe A","Blomqvist Beta A")} 
  else
    {labelSet = c("ML", "Tau-Kendall", "Spearman", "Blomqvist Beta", "Triangular-Unif","Beta-Unif", "Beta-Beta", "Unif-Unif")}
  
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
      columnIndex = if(useDependence) cho_method else if(useAsymtotic) 2:13 else 2:9
      
      for(j in columnIndex) {
        if(j<=9){
          MA = rbind(MA, c(size[i], A[m, j], if(useDependence) w else if(useAsymtotic) namemethod[(j - 1)] else namemethod[(j - 1)]))
        }else if(j>9&&useAsymtotic){
          MA = rbind(MA, c(size[i], A[m+2, j-9+1], if(useDependence) w else if(useAsymtotic) namemethod[(j - 1)]))
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
# Bias = GeneralGraph(5, omega = "-0.75", omega_to_file = omega_to_file, useDependence = FALSE)
# Bias
# ggplot(Bias,aes(x=Size , y=Bias,color = Method))+geom_line()+labs(title="Dependence of -0.75")+
#   scale_x_continuous(breaks=c(10,seq(100,1000,100)))+xlab("Sample size")
# 
# # only graph descriptive measures obtained with the a priori Beta
# Bias = GeneralGraph(5, omega = NULL, omega_to_file = omega_to_file, useDependence = TRUE)
# Bias
# ggplot(Bias,aes(x=Size , y=Bias,color = Dependence))+geom_line()+labs(title=" ")+
#   scale_x_continuous(breaks=c(10,seq(100,1000,100)))+xlab("Sample size")

###################################################################################################
# Empirical marginal distribution
###################################################################################################
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


###################################################################################################
# Empirical joint distribution
###################################################################################################

Joint_Dist = function(X, Y) {
  # Helper function to check if w is less than or equal to z.
  step.my = function(w, z) {
    if (all(w[1] <= z[1],w[2]<=z[2])) {
      1
    } else {
      0
    }
  }
  
  N = nrow(X)  # Number of observations in the data matrix
  fn = NULL  # Initialize the vector for storing empirical distribution values
  
  for (i in 1:nrow(Y)) {
    d = 0  # Initialize the count for the current observation in Y
    for (j in 1:nrow(X)) {
      e = step.my(X[j, ], Y[i, ])  # Apply the step function to each pair (X_j, Y_i)
      d = e + d  # Accumulate the count
    }
    fn[i] = 1 / (N + 1) * d  # Compute the empirical joint distribution value for Y[i, ]
  }
  
  fn  # Return the empirical joint distribution vector
}

Joint_Dist_Parallel = function(X, Y) {
  # Helper function to check if w is less than or equal to z (component-wise)
  step.my = function(w, z) {
    if (all(w[1] <= z[1], w[2] <= z[2])) {
      1
    } else {
      0
    }
  }
  
  N = nrow(X) # Number of observations in the data matrix
  numCores = detectCores() - 1  # Use all but one of the available cores
  cl = makeCluster(numCores)
  
  # Export variables and functions to the cluster environment
  clusterExport(cl, varlist = c("X", "Y", "step.my"), envir = environment())
  
  # Compute empirical values in parallel
  fn = parSapply(cl, 1:nrow(Y), function(i) {
    d = 0
    for (j in 1:nrow(X)) {
      d = d + step.my(X[j, ], Y[i, ])
    }
    1 / (N + 1) * d
  })
  
  stopCluster(cl)
  return(fn)
}