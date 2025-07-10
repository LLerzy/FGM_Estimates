###################################################################################################
# Load necessary functions.
###################################################################################################
source("FunctionsRequired.R")

###################################################################################################
# Obtaining hyperparameter values
###################################################################################################
TSND=Mtovar(-1,-0.5,0.4) #Strong Negative Dependence
TWND=Mtovar(-0.5,0,0.4) #Weak Negative Dependence
TWPD=Mtovar(0,0.5,0.4) #Weak Positive Dependence
TSPD=Mtovar(0.5,1,0.4) #Strong Positive Dependence

###################################################################################################
# MCMC chain convergence monitoring
###################################################################################################

set.seed(123)
trial=rnfgm(-0.75,50)
####
E_ResBeta_B=metropolis_hastings(u=trial[,1], v=trial[,2], 
                                      prior="beta", alpha_prior=TSND$a, beta_prior=TSND$b,  c_prior = TSND$c,
                                      instrumental="beta",precision=20,n_iter=5000,burn_in = 5000, thinning = 5)
E_ResBeta_B$acceptance_rate
Graphs(as.data.frame(E_ResBeta_B$varphi[seq(5000+1,length(E_ResBeta_B$varphi), by = 5)]), 
       "X1", width = 40, lscatt = 0.3, uscatt = 0.3)
####
E_ResBeta_U=metropolis_hastings(u=trial[,1], v=trial[,2], 
                                      prior="beta", alpha_prior=TSND$a, beta_prior=TSND$b,  c_prior = TSND$c,
                                      instrumental="uniform",precision=20,n_iter=5000,burn_in = 5000, thinning = 5)
E_ResBeta_U$acceptance_rate
Graphs(as.data.frame(E_ResBeta_U$varphi[seq(5000+1,length(E_ResBeta_U$varphi), by = 5)]), 
       "X1", width = 40, lscatt = 0.3, uscatt = 0.3)
####
E_ResUnif_U=metropolis_hastings(u=trial[,1], v=trial[,2], 
                                      prior="uniform", alpha_prior=TSND$a, beta_prior=TSND$b, c_prior = TSND$c,
                                      instrumental="uniform",precision=20,n_iter=5000,burn_in = 5000, thinning = 5)
E_ResUnif_U$acceptance_rate
Graphs(as.data.frame(E_ResUnif_U$varphi[seq(5000+1,length(E_ResUnif_U$varphi), by = 5)]), 
       "X1", width = 40, lscatt = 0.3, uscatt = 0.3)

####
E_ResTriang_U=metropolis_hastings(u=trial[,1], v=trial[,2], 
                                        prior="triangular", alpha_prior=TSND$a, beta_prior=TSND$b, c_prior = TSND$c,
                                        instrumental="uniform",precision=20,n_iter=5000,burn_in = 5000, thinning = 5)
E_ResTriang_U$acceptance_rate
Graphs(as.data.frame(E_ResTriang_U$varphi[seq(5000+1,length(E_ResTriang_U$varphi), by = 5)]), 
       "X1", width = 40, lscatt = 0.3, uscatt = 0.3)

###################################################################################################
# Running in parallel and saving results
###################################################################################################
library(parallel)
size=c(10,seq(50,1000,50))
Sizes=c(1:21)
name=as.character(size)

size=c(100)
Sizes=c(1)
name=as.character(size)
###################################################################################################
# Strong Negative Dependence, -0.75
###################################################################################################
n.cores <- round(detectCores() -1) #Number of cores
cl <- makeCluster(n.cores)
clusterExport(cl, c('Sizes','size','name','Sim_Est_N_Samp','Sim_Est_A_Samp','Est_A_Samp','rfgm','rnfgm','safeUroot',
                    'Der_log_lik_FGM','write.xlsx','me','metropolis_hastings','prior_density_triangular',
                    'prior_density_beta','prior_density_uniform','TSND','likelihood','Density_FGM','cor_test','dBeta.4P','rBeta.4P'))
start_time <- Sys.time()

wss <- parSapply(cl = cl,     # Cluster
                 Sizes,
                 function(i){
                   example=Sim_Est_N_Samp(varphi=-0.75,n=size[i],nboot=100, N=1000, confidence=0.95, 
                                          alpha_prior = TSND$a, beta_prior = TSND$b, c_prior = TSND$c,
                                          burn_in = 5000,n_iter = 5000,thinning = 5)
                   write.xlsx(x=example, file="Results-075.xlsx",sheetName =name[i],col.names = TRUE,
                              row.names = TRUE,append = TRUE)
                   return(i)
                 })
finish_time <- Sys.time()
stopCluster(cl)
#Processing time
finish_time - start_time

###################################################################################################
# Negative Weak Dependence, -0.25
###################################################################################################
n.cores <- round(detectCores() -1) #Number of cores
cl <- makeCluster(n.cores)
clusterExport(cl, c('Sizes','size','name','Sim_Est_N_Samp','Sim_Est_A_Samp','Est_A_Samp','rfgm','rnfgm','safeUroot',
                    'Der_log_lik_FGM','write.xlsx','me','metropolis_hastings','prior_density_triangular',
                    'prior_density_beta','prior_density_uniform','TWND','likelihood','Density_FGM','cor_test','dBeta.4P','rBeta.4P'))

start_time <- Sys.time()
wss <- parSapply(cl = cl,     # Cluster
                 Sizes, 
                 function(i){
                   example=Sim_Est_N_Samp(varphi=-0.25,n=size[i],nboot=100, N=1000, confidence=0.95, 
                                          alpha_prior = TWND$a, beta_prior = TWND$b, c_prior = TWND$c,
                                          burn_in = 5000,n_iter = 5000,thinning = 5)
                   write.xlsx(x=example, file="Results-025.xlsx",sheetName =name[i],col.names = TRUE,
                              row.names = TRUE,append = TRUE)
                   return(i)
                 })
finish_time <- Sys.time()
stopCluster(cl)

#Processing time
finish_time - start_time

###################################################################################################
# Positive Weak Dependence, 0.25
###################################################################################################
n.cores <- round(detectCores() -1) #Number of cores
cl <- makeCluster(n.cores)
clusterExport(cl, c('Sizes','size','name','Sim_Est_N_Samp','Sim_Est_A_Samp','Est_A_Samp','rfgm','rnfgm','safeUroot',
                    'Der_log_lik_FGM','write.xlsx','me','metropolis_hastings','prior_density_triangular',
                    'prior_density_beta','prior_density_uniform','TWPD','likelihood','Density_FGM','cor_test','dBeta.4P','rBeta.4P'))

start_time <- Sys.time()
wss <- parSapply(cl = cl,     # Cluster
                 Sizes, 
                 function(i){
                   example=Sim_Est_N_Samp(varphi=0.25,n=size[i],nboot=100, N=1000, confidence=0.95, 
                                          alpha_prior = TWPD$a, beta_prior = TWPD$b, c_prior = TWPD$c,
                                          burn_in = 5000,n_iter = 5000,thinning = 5)
                   write.xlsx(x=example, file="Results025.xlsx",sheetName =name[i],col.names = TRUE,
                              row.names = TRUE,append = TRUE)
                   return(i)
                 })
finish_time <- Sys.time()
stopCluster(cl)

#Processing time
finish_time - start_time

###################################################################################################
# Positive Strong Dependence, 0.75
###################################################################################################

n.cores <- round(detectCores() -1) #Number of cores
cl <- makeCluster(n.cores)
clusterExport(cl, c('Sizes','size','name','Sim_Est_N_Samp','Sim_Est_A_Samp','Est_A_Samp','rfgm','rnfgm','safeUroot',
                    'Der_log_lik_FGM','write.xlsx','me','metropolis_hastings','prior_density_triangular',
                    'prior_density_beta','prior_density_uniform','TSPD','likelihood','Density_FGM','cor_test','dBeta.4P','rBeta.4P'))

start_time <- Sys.time()
wss <- parSapply(cl = cl,     # Cluster
                 Sizes, 
                 function(i){
                   example=Sim_Est_N_Samp(varphi=0.75,n=size[i],nboot=100, N=1000, confidence=0.95, 
                                          alpha_prior = TSPD$a, beta_prior = TSPD$b, c_prior = TSPD$c,
                                          burn_in = 5000,n_iter = 5000,thinning = 5)
                   write.xlsx(x=example, file="Results075.xlsx",sheetName =name[i],col.names = TRUE,
                              row.names = TRUE,append = TRUE)
                   return(i)
                 })
finish_time <- Sys.time()
stopCluster(cl)

#Processing time
finish_time - start_time

###################################################################################################
# Graph of Results Obtained.
###################################################################################################
omega_to_file = list(
  "-0.75" = "Results-075.xlsx",
  "-0.25" = "Results-025.xlsx",
  "0.25" = "Results025.xlsx",
  "0.75" = "Results075.xlsx")

# Plot for Bias
Bias=GeneralGraph(5, omega = "-0.75", omega_to_file = omega_to_file, useDependence = FALSE)
Biasn075=ggplot(Bias,aes(x=Size , y=Bias,color = Method, shape = Method, linetype = Method))+geom_line()+geom_point(size = 2)+
  labs(title="a. Dependence of -0.75")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+xlab("") + ylab("Bias")+
  theme(legend.position = "none") +
  scale_shape_manual(values = c(16, 17, 18, 19, 15, 8, 5, 6))+
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "dotted", "twodash"))

Bias=GeneralGraph(5, omega = "-0.25", omega_to_file = omega_to_file, useDependence = FALSE)
Biasn025=ggplot(Bias,aes(x=Size , y=Bias,color = Method, shape = Method, linetype = Method))+geom_line()+geom_point(size = 2)+
  labs(title="b. Dependence of -0.25")+
  theme(legend.position = "none") +
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+xlab("") + ylab("")+
  scale_shape_manual(values = c(16, 17, 18, 19, 15, 8, 5, 6))+
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "dotted", "twodash"))

Bias=GeneralGraph(5, omega = "0.25", omega_to_file = omega_to_file, useDependence = FALSE)
Bias025=ggplot(Bias,aes(x=Size , y=Bias,color = Method, shape = Method, linetype = Method))+geom_line()+geom_point(size = 2)+
  labs(title="c. Dependence of 0.25")+
  theme(legend.position = "none") +
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+xlab("Sample size") + ylab("Bias")+
  scale_shape_manual(values = c(16, 17, 18, 19, 15, 8, 5, 6))+
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "dotted", "twodash"))

Bias=GeneralGraph(5, omega = "0.75", omega_to_file = omega_to_file, useDependence = FALSE)
Bias075=ggplot(Bias,aes(x=Size , y=Bias,color = Method, shape = Method, linetype = Method))+geom_line()+geom_point(size = 2)+
  labs(title="d. Dependence of 0.75")+
  theme(legend.position = "none") +
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+xlab("Sample size") + ylab("")+
  scale_shape_manual(values = c(16, 17, 18, 19, 15, 8, 5, 6))+
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "dotted", "twodash"))

# Function to extract the legend from a plot
g_legend <- function(a.gplot) {
  tmp <- ggplotGrob(a.gplot)
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Extract the legend from the mean plot
legend <- g_legend(Biasn075 + theme(legend.position = "right", legend.title = element_blank()))

# Combine the plots without legend with the legend on the right
grid.arrange(
  arrangeGrob(Biasn075, Biasn025, Bias025, Bias075, nrow=2, ncol = 2, heights = c(9, 9)),
  legend, 
  ncol = 2, 
  widths = c(7, 1)
)

# Plot for MSE
MSE=GeneralGraph(6, omega = "-0.75", omega_to_file = omega_to_file, useDependence = FALSE)
Msen075=ggplot(MSE,aes(x=Size , y=MSE,color = Method, shape = Method, linetype = Method))+geom_line()+geom_point(size = 2)+
  labs(title="a. Dependence of -0.75") +
  theme(legend.position = "none") +
  scale_x_continuous(breaks=c(10,seq(100,1000,100))) + xlab("") + ylab("MSE") +
  scale_shape_manual(values = c(16, 17, 18, 19, 15, 8, 5, 6))+
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "dotted", "twodash"))

MSE=GeneralGraph(6, omega = "-0.25", omega_to_file = omega_to_file, useDependence = FALSE) 
Msen025=ggplot(MSE,aes(x=Size , y=MSE,color = Method, shape = Method, linetype = Method))+geom_line()+geom_point(size = 2)+
  labs(title="b. Dependence of -0.25") +
  theme(legend.position = "none") +
  scale_x_continuous(breaks=c(10,seq(100,1000,100))) + xlab("") + ylab("") +
  scale_shape_manual(values = c(16, 17, 18, 19, 15, 8, 5, 6))+
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "dotted", "twodash"))

MSE=GeneralGraph(6, omega = "0.25", omega_to_file = omega_to_file, useDependence = FALSE) 
Mse025=ggplot(MSE,aes(x=Size , y=MSE,color = Method, shape = Method, linetype = Method))+geom_line()+geom_point(size = 2)+
  labs(title="c. Dependence of 0.25") +
  theme(legend.position = "none") +
  scale_x_continuous(breaks=c(10,seq(100,1000,100))) + xlab("Sample size") +
  scale_shape_manual(values = c(16, 17, 18, 19, 15, 8, 5, 6))+
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "dotted", "twodash"))

MSE=GeneralGraph(6, omega = "0.75", omega_to_file = omega_to_file, useDependence = FALSE) 
Mse075=ggplot(MSE,aes(x=Size , y=MSE,color = Method, shape = Method, linetype = Method))+geom_line()+geom_point(size = 2)+
  labs(title="d. Dependence of 0.75") +
  theme(legend.position = "none") +
  scale_x_continuous(breaks=c(10,seq(100,1000,100))) + xlab("Sample size") + ylab("") +
  scale_shape_manual(values = c(16, 17, 18, 19, 15, 8, 5, 6))+
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "dotted", "twodash"))

# Extract the legend from the mean plot
legend <- g_legend(Msen075 + theme(legend.position = "right", legend.title = element_blank()))

# Combine the plots without legend with the legend on the right
grid.arrange(
  arrangeGrob(Msen075, Msen025, Mse025, Mse075, nrow=2, ncol = 2, heights = c(9, 9)),
  legend, 
  ncol = 2, 
  widths = c(7, 1)
)

# Plot for Average Length
Length=GeneralGraph(7, omega = "-0.75", omega_to_file = omega_to_file, useDependence = FALSE,useAsymtotic = TRUE)
lengn075=ggplot(Length,aes(x=Size , y=Length,color = Method, shape = Method, linetype = Method))+geom_line()+geom_point(size = 2) +
  labs(title="a. Dependence of -0.75") +
  theme(legend.position = "none") +
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Average length") + xlab("")+
  scale_shape_manual(values = c(16, 17, 18, 19, 15, 8, 5, 6, 4, 3 , 2, 1))+
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "dotted", "twodash", "dashed", "dotted", "dotdash","dashed"))

Length=GeneralGraph(7, omega = "-0.25", omega_to_file = omega_to_file, useDependence = FALSE,useAsymtotic = TRUE)
lengn025=ggplot(Length,aes(x=Size , y=Length,color = Method, shape = Method, linetype = Method))+geom_line()+geom_point(size = 2)+
  labs(title="b. Dependence of -0.25") +
  theme(legend.position = "none") +
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Average length") + xlab("") + ylab("") +
  scale_shape_manual(values = c(16, 17, 18, 19, 15, 8, 5, 6, 4, 3 , 2, 1))+
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "dotted", "twodash", "dashed", "dotted", "dotdash","dashed"))

Length=GeneralGraph(7, omega = "0.25", omega_to_file = omega_to_file, useDependence = FALSE,useAsymtotic = TRUE)
leng025=ggplot(Length,aes(x=Size , y=Length,color = Method, shape = Method, linetype = Method))+geom_line()+geom_point(size = 2)+
  labs(title="c. Dependence of -0.25") +
  theme(legend.position = "none") +
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Average length") + xlab("Sample size") +
  scale_shape_manual(values = c(16, 17, 18, 19, 15, 8, 5, 6, 4, 3 , 2, 1))+
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "dotted", "twodash", "dashed", "dotted", "dotdash","dashed"))

Length=GeneralGraph(7, omega = "0.75", omega_to_file = omega_to_file, useDependence = FALSE,useAsymtotic = TRUE)
leng075=ggplot(Length,aes(x=Size , y=Length,color = Method, shape = Method, linetype = Method))+geom_line()+geom_point(size = 2)+
  labs(title="d. Dependence of 0.75") +
  theme(legend.position = "none") +
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Average length") + xlab("Sample size") + ylab("") +
  scale_shape_manual(values = c(16, 17, 18, 19, 15, 8, 5, 6, 4, 3 , 2, 1))+
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "dotted", "twodash", "dashed", "dotted", "dotdash","dashed"))

# Extract the legend from the mean plot
legend <- g_legend(lengn075 + theme(legend.position = "right", legend.title = element_blank()))

# Combine the plots without legend with the legend on the right
grid.arrange(
  arrangeGrob(lengn075, lengn025, leng025, leng075, nrow=2, ncol = 2, heights = c(9, 9)),
  legend, 
  ncol = 2, 
  widths = c(7, 1)
)

# Plot for Probability Coverage
Coverage=GeneralGraph(8, omega = "-0.75", omega_to_file = omega_to_file, useDependence = FALSE,useAsymtotic = TRUE)
Coveragen075=ggplot(Coverage,aes(x=Size , y=Coverage,color = Method, shape = Method, linetype = Method))+geom_line()+geom_point(size = 2)+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Probability") + xlab("") +
  labs(title="a. Dependence of -0.75") +
  theme(legend.position = "none") +
  scale_shape_manual(values = c(16, 17, 18, 19, 15, 8, 5, 6, 4, 3 , 2, 1))+
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "dotted", "twodash", "dashed", "dotted", "dotdash", "dashed"))

Coverage=GeneralGraph(8, omega = "-0.25", omega_to_file = omega_to_file, useDependence = FALSE,useAsymtotic = TRUE)
Coveragen025=ggplot(Coverage,aes(x=Size , y=Coverage,color = Method, shape = Method, linetype = Method))+geom_line()+geom_point(size = 2)+
  labs(title="b. Dependence of -0.25") +
  theme(legend.position = "none") +
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Probability") + xlab("") + ylab("") +
  scale_shape_manual(values = c(16, 17, 18, 19, 15, 8, 5, 6, 4, 3 , 2,1))+
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "dotted", "twodash", "dashed", "dotted", "dotdash","dashed"))

Coverage=GeneralGraph(8, omega = "0.25", omega_to_file = omega_to_file, useDependence = FALSE,useAsymtotic = TRUE)
Coverage025=ggplot(Coverage,aes(x=Size , y=Coverage,color = Method, shape = Method, linetype = Method))+geom_line()+geom_point(size = 2)+
  labs(title="c. Dependence of 0.25") +
  theme(legend.position = "none") +
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Probability") + xlab("Sample size") +
  scale_shape_manual(values = c(16, 17, 18, 19, 15, 8, 5, 6, 4, 3 , 2,1))+
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "dotted", "twodash", "dashed", "dotted", "dotdash","dashed"))

Coverage=GeneralGraph(8, omega = "0.75", omega_to_file = omega_to_file, useDependence = FALSE,useAsymtotic = TRUE)
Coverage075=ggplot(Coverage,aes(x=Size , y=Coverage,color = Method, shape = Method, linetype = Method))+geom_line()+geom_point(size = 2)+
  labs(title="d. Dependence of 0.75") +
  theme(legend.position = "none") +
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Probability") + xlab("Sample size") + ylab("") +
  scale_shape_manual(values = c(16, 17, 18, 19, 15, 8, 5, 6, 4, 3 , 2,1)) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "dotted", "twodash", "dashed", "dotted", "dotdash","dashed"))

# Extract the legend from the mean plot
legend <- g_legend(Coveragen075 + theme(legend.position = "right", legend.title = element_blank()))

# Combine the plots without legend with the legend on the right
grid.arrange(
  arrangeGrob(Coveragen075, Coveragen025, Coverage025, Coverage075, nrow=2, ncol = 2, heights = c(9, 9)),
  legend, 
  ncol = 2, 
  widths = c(7, 1)
)

###################################################################################################
#Graphs for the Beta prior distribution
###################################################################################################

# Plot for Bias
Bias=GeneralGraph(5, omega = NULL, omega_to_file = omega_to_file, useDependence = TRUE,cho_method = 8)
BiasBeta=ggplot(Bias,aes(x=Size , y=Bias,color = Dependence, shape = Dependence, linetype = Dependence))+geom_line() + geom_point(size = 2) +
  labs(title="") + theme(legend.position = "none") +
  scale_x_continuous(breaks=c(10,seq(100,1000,100))) + xlab("") +
  scale_shape_manual(values = c(16, 17, 18, 19)) +
  scale_linetype_manual(values = c("longdash", "twodash", "dotted", "twodash"))


# Plot for MSE
MSE=GeneralGraph(6, omega = NULL, omega_to_file = omega_to_file, useDependence = TRUE,cho_method = 8) #just to initialize.
MsemBeta=ggplot(MSE,aes(x=Size , y=MSE , color = Dependence, shape = Dependence, linetype = Dependence))+geom_line() + geom_point(size = 2) +
  labs(title="") + theme(legend.position = "none") +
  scale_x_continuous(breaks=c(10,seq(100,1000,100))) + xlab("") +
  scale_shape_manual(values = c(16, 17, 18, 19)) +
  scale_linetype_manual(values = c("longdash", "twodash", "dotted", "twodash"))

# Plot for Average Length
Length=GeneralGraph(7, omega = NULL, omega_to_file = omega_to_file, useDependence = TRUE,cho_method = 8)
lengmBeta=ggplot(Length,aes(x=Size , y=Length,color = Dependence, shape = Dependence, linetype = Dependence))+geom_line() + geom_point(size = 2) +
  labs(title="") + theme(legend.position = "none") +
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Average length") + xlab("Sample size") +
  scale_shape_manual(values = c(16, 17, 18, 19)) +
  scale_linetype_manual(values = c("longdash", "twodash", "dotted", "twodash"))

# Plot for Probability Coverage
Coverage=GeneralGraph(8, omega = NULL, omega_to_file = omega_to_file, useDependence = TRUE,cho_method = 8)
CoveragemBeta=ggplot(Coverage,aes(x=Size , y=Coverage,color = Dependence, shape = Dependence, linetype = Dependence))+geom_line() + geom_point(size = 2) +
  labs(title="") + theme(legend.position = "none") +
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Probability") + xlab("Sample size") +
  scale_shape_manual(values = c(16, 17, 18, 19)) +
  scale_linetype_manual(values = c("longdash", "twodash", "dotted", "twodash"))

# Extract the legend from the mean plot
legend <- g_legend(BiasBeta + theme(legend.position = "right", legend.title = element_blank()))

# Combine the plots without legend with the legend on the right
grid.arrange(
  arrangeGrob(BiasBeta, MsemBeta, lengmBeta, CoveragemBeta, nrow=2, ncol = 2, heights = c(9, 9)),
  legend, 
  ncol = 2, 
  widths = c(7, 1)
)
