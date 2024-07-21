###############################################
# Load necessary functions.
###############################################
source("FunctionsRequired.R")

####################################################################
# Obtaining hyperparameter values
####################################################################
TSND=Mtovar(-1,-0.5,0.4) #Strong Negative Dependence
TWND=Mtovar(-0.5,0,0.4) #Weak Negative Dependence
TWPD=Mtovar(0,0.5,0.4) #Weak Positive Dependence
TSPD=Mtovar(0.5,1,0.4) #Strong Positive Dependence

####################################################################
# Running in parallel and saving results
####################################################################
library(parallel)
size=c(10,seq(50,1000,50))
Sizes=c(1:21)
name=as.character(size)

####################################
# Strong Negative Dependence, -0.75
####################################
n.cores <- round(detectCores() -1) #Number of cores
cl <- makeCluster(n.cores)
clusterExport(cl, c('Sizes','size','name','Sim_Est_N_Samp','Sim_Est_A_Samp','Est_A_Samp','rfgm','rnfgm','safeUroot',
                    'Der_log_lik_FGM','write.xlsx','me','metropolis_hastings','prior_density_triangular',
                    'prior_density_beta','prior_density_uniform','TSND','likelihood','Density_FGM'))
start_time <- Sys.time()

wss <- parSapply(cl = cl,     # Cluster
                 Sizes,
                 function(i){
                   example=Sim_Est_N_Samp(varphi=-0.75,n=size[i],nboot=100, N=1000, confidence=0.95, 
                                          alpha_prior = TSND$a, beta_prior = TSND$b, c_prior = TSND$c,
                                          n_burne = 10000,n_iter = 30000,n_thin = 5)
                   write.xlsx(x=example, file="Results-075.xlsx",sheetName =name[i],col.names = TRUE,
                              row.names = TRUE,append = TRUE)
                   return(i)
                 })
finish_time <- Sys.time()
stopCluster(cl)
#Processing time
finish_time - start_time

####################################
# Negative Weak Dependence, -0.25
####################################
n.cores <- round(detectCores() -1) #Number of cores
cl <- makeCluster(n.cores)
clusterExport(cl, c('Sizes','size','name','Sim_Est_N_Samp','Sim_Est_A_Samp','Est_A_Samp','rfgm','rnfgm','safeUroot',
                    'Der_log_lik_FGM','write.xlsx','me','metropolis_hastings','prior_density_triangular',
                    'prior_density_beta','prior_density_uniform','TWND','likelihood','Density_FGM'))

start_time <- Sys.time()
wss <- parSapply(cl = cl,     # Cluster
                 Sizes, 
                 function(i){
                   example=Sim_Est_N_Samp(varphi=-0.25,n=size[i],nboot=100, N=1000, confidence=0.95, 
                                          alpha_prior = TWND$a, beta_prior = TWND$b, c_prior = TWND$c,
                                          n_burne = 10000,n_iter = 30000,n_thin = 5)
                   write.xlsx(x=example, file="Results-025.xlsx",sheetName =name[i],col.names = TRUE,
                              row.names = TRUE,append = TRUE)
                   return(i)
                 })
finish_time <- Sys.time()
stopCluster(cl)

#Processing time
finish_time - start_time

####################################
# Positive Weak Dependence, 0.25
####################################
n.cores <- round(detectCores() -1) #Number of cores
cl <- makeCluster(n.cores)
clusterExport(cl, c('Sizes','size','name','Sim_Est_N_Samp','Sim_Est_A_Samp','Est_A_Samp','rfgm','rnfgm','safeUroot',
                    'Der_log_lik_FGM','write.xlsx','me','metropolis_hastings','prior_density_triangular',
                    'prior_density_beta','prior_density_uniform','TWPD','likelihood','Density_FGM'))

start_time <- Sys.time()
wss <- parSapply(cl = cl,     # Cluster
                 Sizes, 
                 function(i){
                   example=Sim_Est_N_Samp(varphi=0.25,n=size[i],nboot=100, N=1000, confidence=0.95, 
                                          alpha_prior = TWPD$a, beta_prior = TWPD$b, c_prior = TWPD$c,
                                          n_burne = 10000,n_iter = 30000,n_thin = 5)
                   write.xlsx(x=example, file="Results025.xlsx",sheetName =name[i],col.names = TRUE,
                              row.names = TRUE,append = TRUE)
                   return(i)
                 })
finish_time <- Sys.time()
stopCluster(cl)

#Processing time
finish_time - start_time

####################################
# Positive Strong Dependence, 0.75
####################################

n.cores <- round(detectCores() -1) #Number of cores
cl <- makeCluster(n.cores)
clusterExport(cl, c('Sizes','size','name','Sim_Est_N_Samp','Sim_Est_A_Samp','Est_A_Samp','rfgm','rnfgm','safeUroot',
                    'Der_log_lik_FGM','write.xlsx','me','metropolis_hastings','prior_density_triangular',
                    'prior_density_beta','prior_density_uniform','TSPD','likelihood','Density_FGM'))

start_time <- Sys.time()
wss <- parSapply(cl = cl,     # Cluster
                 Sizes, 
                 function(i){
                   example=Sim_Est_N_Samp(varphi=0.75,n=size[i],nboot=100, N=1000, confidence=0.95, 
                                          alpha_prior = TSPD$a, beta_prior = TSPD$b, c_prior = TSPD$c,
                                          n_burne = 10000,n_iter = 30000,n_thin = 5)
                   write.xlsx(x=example, file="Results075.xlsx",sheetName =name[i],col.names = TRUE,
                              row.names = TRUE,append = TRUE)
                   return(i)
                 })
finish_time <- Sys.time()
stopCluster(cl)

#Processing time
finish_time - start_time

####################################
# Graph of Results Obtained.
####################################
omega_to_file = list(
  "-0.75" = "Results-075.xlsx",
  "-0.25" = "Results-025.xlsx",
  "0.25" = "Results025.xlsx",
  "0.75" = "Results075.xlsx")

# Plot for Bias
Bias=GeneralGraph(5, omega = "-0.75", omega_to_file = omega_to_file, useDependence = FALSE)
Biasn075=ggplot(Bias,aes(x=Size , y=Bias,color = Method))+geom_line()+labs(title="Dependence of -0.75")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+xlab("Sample size")

Bias=GeneralGraph(5, omega = "-0.25", omega_to_file = omega_to_file, useDependence = FALSE)
Biasn025=ggplot(Bias,aes(x=Size , y=Bias,color = Method))+geom_line()+labs(title="Dependence of -0.25")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+xlab("Sample size")

Bias=GeneralGraph(5, omega = "0.25", omega_to_file = omega_to_file, useDependence = FALSE)
Bias025=ggplot(Bias,aes(x=Size , y=Bias,color = Method))+geom_line()+labs(title="Dependence of 0.25")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+xlab("Sample size")

Bias=GeneralGraph(5, omega = "0.75", omega_to_file = omega_to_file, useDependence = FALSE)
Bias075=ggplot(Bias,aes(x=Size , y=Bias,color = Method))+geom_line()+labs(title="Dependence of 0.75")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+xlab("Sample size")

grid.arrange(Biasn075, Biasn025, Bias025, Bias075,nrow = 2,ncol=2,layout_matrix = rbind(c(1, 2),c(3, 4)))

# Plot for MSE
MSE=GeneralGraph(6, omega = "-0.75", omega_to_file = omega_to_file, useDependence = FALSE)
Msen075=ggplot(MSE,aes(x=Size , y=MSE,color = Method))+geom_line()+labs(title="Dependence of -0.75")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+xlab("Sample size")

MSE=GeneralGraph(6, omega = "-0.25", omega_to_file = omega_to_file, useDependence = FALSE) 
Msen025=ggplot(MSE,aes(x=Size , y=MSE,color = Method))+geom_line()+labs(title="Dependence of -0.25")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+xlab("Sample size")

MSE=GeneralGraph(6, omega = "0.25", omega_to_file = omega_to_file, useDependence = FALSE) 
Mse025=ggplot(MSE,aes(x=Size , y=MSE,color = Method))+geom_line()+labs(title="Dependence of 0.25")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+xlab("Sample size")

MSE=GeneralGraph(6, omega = "0.75", omega_to_file = omega_to_file, useDependence = FALSE) 
Mse075=ggplot(MSE,aes(x=Size , y=MSE,color = Method))+geom_line()+labs(title="Dependence of 0.75")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+xlab("Sample size")

grid.arrange(Msen075, Msen025, Mse025, Mse075, nrow = 2,ncol=2,layout_matrix = rbind(c(1, 2),
                                                                                     c(3, 4)))

# Plot for Average Length
Length=GeneralGraph(7, omega = "-0.75", omega_to_file = omega_to_file, useDependence = FALSE,useAsymtotic = TRUE)
lengn075=ggplot(Length,aes(x=Size , y=Length,color = Method))+geom_line()+labs(title="Dependence of -0.75")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Average length")+xlab("Sample size")
Length=GeneralGraph(7, omega = "-0.25", omega_to_file = omega_to_file, useDependence = FALSE,useAsymtotic = TRUE)
lengn025=ggplot(Length,aes(x=Size , y=Length,color = Method))+geom_line()+labs(title="Dependence of -0.25")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Average length")+xlab("Sample size")
Length=GeneralGraph(7, omega = "0.25", omega_to_file = omega_to_file, useDependence = FALSE,useAsymtotic = TRUE)
leng025=ggplot(Length,aes(x=Size , y=Length,color = Method))+geom_line()+labs(title="Dependence of 0.25")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Average length")+xlab("Sample size")
Length=GeneralGraph(7, omega = "0.75", omega_to_file = omega_to_file, useDependence = FALSE,useAsymtotic = TRUE)
leng075=ggplot(Length,aes(x=Size , y=Length,color = Method))+geom_line()+labs(title="Dependence of 0.75")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Average length")+xlab("Sample size")

grid.arrange(lengn075, lengn025, leng025, leng075, nrow = 2,ncol=2,layout_matrix = rbind(c(1, 2),
                                                                                         c(3, 4)))

# Plot for Probability Coverage
Coverage=GeneralGraph(8, omega = "-0.75", omega_to_file = omega_to_file, useDependence = FALSE,useAsymtotic = TRUE)
Coveragen075=ggplot(Coverage,aes(x=Size , y=Coverage,color = Method))+geom_line()+labs(title="Dependence of -0.75")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Probability")+xlab("Sample size")
Coverage=GeneralGraph(8, omega = "-0.25", omega_to_file = omega_to_file, useDependence = FALSE,useAsymtotic = TRUE)
Coveragen025=ggplot(Coverage,aes(x=Size , y=Coverage,color = Method))+geom_line()+labs(title="Dependence of -0.25")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Probability")+xlab("Sample size")
Coverage=GeneralGraph(8, omega = "0.25", omega_to_file = omega_to_file, useDependence = FALSE,useAsymtotic = TRUE)
Coverage025=ggplot(Coverage,aes(x=Size , y=Coverage,color = Method))+geom_line()+labs(title="Dependence of 0.25")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Probability")+xlab("Sample size")
Coverage=GeneralGraph(8, omega = "0.75", omega_to_file = omega_to_file, useDependence = FALSE,useAsymtotic = TRUE)
Coverage075=ggplot(Coverage,aes(x=Size , y=Coverage,color = Method))+geom_line()+labs(title="Dependence of 0.75")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Probability")+xlab("Sample size")

grid.arrange(Coveragen075, Coveragen025, Coverage025, Coverage075, nrow = 2,ncol=2,layout_matrix = rbind(c(1, 2),
                                                                                                         c(3, 4)))


###############################################################################
#Graphs for the Beta prior distribution
###############################################################################

# Plot for Bias
Bias=GeneralGraph(5, omega = NULL, omega_to_file = omega_to_file, useDependence = TRUE)
BiasBeta=ggplot(Bias,aes(x=Size , y=Bias,color = Dependence))+geom_line()+labs(title=" ")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+xlab("Sample size")

# Plot for MSE
MSE=GeneralGraph(6, omega = NULL, omega_to_file = omega_to_file, useDependence = TRUE) #just to initialize.
MsemBeta=ggplot(MSE,aes(x=Size , y=MSE , color = Dependence))+geom_line()+labs(title=" ")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+xlab("Sample size")

# Plot for Average Length
Length=GeneralGraph(7, omega = NULL, omega_to_file = omega_to_file, useDependence = TRUE)
lengmBeta=ggplot(Length,aes(x=Size , y=Length,color = Dependence))+geom_line()+labs(title=" ")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Average length")+xlab("Sample size")

# Plot for Probability Coverage
Coverage=GeneralGraph(8, omega = NULL, omega_to_file = omega_to_file, useDependence = TRUE)
CoveragemBeta=ggplot(Coverage,aes(x=Size , y=Coverage,color = Dependence))+geom_line()+labs(title=" ")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Probability")+xlab("Sample size")

grid.arrange(BiasBeta, MsemBeta, lengmBeta, CoveragemBeta, nrow = 2,ncol=2,layout_matrix = rbind(c(1, 2),
                                                                                                 c(3, 4)))
