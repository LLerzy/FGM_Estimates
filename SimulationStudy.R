###############################################
# Load necessary functions.
###############################################
source("FunctionsRequired.R")

####################################################################
# Obtaining hyperparameter values
####################################################################
Mtovar(-1,-0.5,0.05) #Strong Negative Dependence
Mtovar(-0.5,0,0.05) #Weak Negative Dependence
Mtovar(0,0.5,0.05) #Weak Positive Dependence
Mtovar(0.5,1,0.05) #Strong Positive Dependence

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
n.cores <- round(detectCores() / 2) #Number of cores
cl <- makeCluster(n.cores)
clusterExport(cl, c('Sizes','size','name','Sim_Est_N_Samp','Sim_Est_A_Samp','Est_A_Samp','rfgm','rnfgm','safeUroot',
                    'Der_log_lik_FGM','me','post_est','bugs','write.xlsx'))
start_time <- Sys.time()

wss <- parSapply(cl = cl,     # Cluster
                 Sizes,
                 function(i){
                   example=Sim_Est_N_Samp(varphi=-0.75,n=size[i],nboot=200, N=1000, confidence=0.95, 
                                          a=17.375, b=121.625, c=-0.75,n_burnd = 100,n_iter = 2000,n_thin = 10)
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
n.cores <- round(detectCores() / 2) #Number of cores
cl <- makeCluster(n.cores)
clusterExport(cl, c('Sizes','size','name','Sim_Est_N_Samp','Sim_Est_A_Samp','Est_A_Samp','rfgm','rnfgm','safeUroot',
                    'Der_log_lik_FGM','me','post_est','bugs','write.xlsx'))

start_time <- Sys.time()
wss <- parSapply(cl = cl,     # Cluster
                 Sizes, 
                 function(i){
                   example=Sim_Est_N_Samp(varphi=-0.25,n=size[i],nboot=200, N=1000, confidence=0.95, 
                                          a=112.125, b=186.875, c=-0.25,n_burnd = 100,n_iter = 2000,n_thin = 10)
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
n.cores <- round(detectCores() / 2) #Number of cores
cl <- makeCluster(n.cores)
clusterExport(cl, c('Sizes','size','name','Sim_Est_N_Samp','Sim_Est_A_Samp','Est_A_Samp','rfgm','rnfgm','safeUroot',
                    'Der_log_lik_FGM','me','post_est','bugs','write.xlsx'))

start_time <- Sys.time()
wss <- parSapply(cl = cl,     # Cluster
                 Sizes, 
                 function(i){
                   example=Sim_Est_N_Samp(varphi=0.25,n=size[i],nboot=200, N=1000, confidence=0.95, 
                                          a=186.875, b=112.125, c=0.25,n_burnd = 100,n_iter = 2000,n_thin = 10)
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

n.cores <- round(detectCores() / 2) #Number of cores
cl <- makeCluster(n.cores)
clusterExport(cl, c('Sizes','size','name','Sim_Est_N_Samp','Sim_Est_A_Samp','Est_A_Samp','rfgm','rnfgm','safeUroot',
                    'Der_log_lik_FGM','me','post_est','bugs','write.xlsx'))

start_time <- Sys.time()
wss <- parSapply(cl = cl,     # Cluster
                 Sizes, 
                 function(i){
                   example=Sim_Est_N_Samp(varphi=0.75,n=size[i],nboot=200, N=1000, confidence=0.95, 
                                          a=121.625, b=17.375, c=0.75,n_burnd = 100,n_iter = 2000,n_thin = 10)
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
  "-0.75" = "Full-Results-075.xlsx",
  "-0.25" = "Full-Results-025.xlsx",
  "0.25" = "Full-Results025.xlsx",
  "0.75" = "Full-Results075.xlsx")

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
MSE=GeneralGraph(3, omega = "-0.75", omega_to_file = omega_to_file, useDependence = FALSE) #just to initialize.
MSE[,2]=MSE$SD^2+GeneralGraph(5, omega = "-0.75", omega_to_file = omega_to_file, useDependence = FALSE)$Bias^2
colnames(MSE)=list("Size","MSE","Method")
Msen075=ggplot(MSE,aes(x=Size , y=MSE,color = Method))+geom_line()+labs(title="Dependence of -0.75")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+xlab("Sample size")

MSE=GeneralGraph(3, omega = "-0.25", omega_to_file = omega_to_file, useDependence = FALSE) #just to initialize.
MSE[,2]=MSE$SD^2+GeneralGraph(5, omega = "-0.25", omega_to_file = omega_to_file, useDependence = FALSE)$Bias^2
colnames(MSE)=list("Size","MSE","Method")
Msen025=ggplot(MSE,aes(x=Size , y=MSE,color = Method))+geom_line()+labs(title="Dependence of -0.25")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+xlab("Sample size")

MSE=GeneralGraph(3, omega = "0.25", omega_to_file = omega_to_file, useDependence = FALSE) #just to initialize.
MSE[,2]=MSE$SD^2+GeneralGraph(5, omega = "0.25", omega_to_file = omega_to_file, useDependence = FALSE)$Bias^2
colnames(MSE)=list("Size","MSE","Method")
Mse025=ggplot(MSE,aes(x=Size , y=MSE,color = Method))+geom_line()+labs(title="Dependence of 0.25")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+xlab("Sample size")

MSE=GeneralGraph(3, omega = "0.75", omega_to_file = omega_to_file, useDependence = FALSE) #just to initialize.
MSE[,2]=MSE$SD^2+GeneralGraph(5, omega = "0.75", omega_to_file = omega_to_file, useDependence = FALSE)$Bias^2
colnames(MSE)=list("Size","MSE","Method")
Mse075=ggplot(MSE,aes(x=Size , y=MSE,color = Method))+geom_line()+labs(title="Dependence of 0.75")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+xlab("Sample size")

grid.arrange(Msen075, Msen025, Mse025, Mse075, nrow = 2,ncol=2,layout_matrix = rbind(c(1, 2),
                                                                                     c(3, 4)))

# Plot for Average Length
Length=GeneralGraph(6, omega = "-0.75", omega_to_file = omega_to_file, useDependence = FALSE,useAsymtotic = TRUE)
lengn075=ggplot(Length,aes(x=Size , y=Length,color = Method))+geom_line()+labs(title="Dependence of -0.75")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Average length")+xlab("Sample size")
Length=GeneralGraph(6, omega = "-0.25", omega_to_file = omega_to_file, useDependence = FALSE,useAsymtotic = TRUE)
lengn025=ggplot(Length,aes(x=Size , y=Length,color = Method))+geom_line()+labs(title="Dependence of -0.25")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Average length")+xlab("Sample size")
Length=GeneralGraph(6, omega = "0.25", omega_to_file = omega_to_file, useDependence = FALSE,useAsymtotic = TRUE)
leng025=ggplot(Length,aes(x=Size , y=Length,color = Method))+geom_line()+labs(title="Dependence of 0.25")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Average length")+xlab("Sample size")
Length=GeneralGraph(6, omega = "0.75", omega_to_file = omega_to_file, useDependence = FALSE,useAsymtotic = TRUE)
leng075=ggplot(Length,aes(x=Size , y=Length,color = Method))+geom_line()+labs(title="Dependence of 0.75")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Average length")+xlab("Sample size")

grid.arrange(lengn075, lengn025, leng025, leng075, nrow = 2,ncol=2,layout_matrix = rbind(c(1, 2),
                                                                                         c(3, 4)))

# Plot for Probability Coverage
Coverage=GeneralGraph(7, omega = "-0.75", omega_to_file = omega_to_file, useDependence = FALSE,useAsymtotic = TRUE)
Coveragen075=ggplot(Coverage,aes(x=Size , y=Coverage,color = Method))+geom_line()+labs(title="Dependence of -0.75")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Probability")+xlab("Sample size")
Coverage=GeneralGraph(7, omega = "-0.25", omega_to_file = omega_to_file, useDependence = FALSE,useAsymtotic = TRUE)
Coveragen025=ggplot(Coverage,aes(x=Size , y=Coverage,color = Method))+geom_line()+labs(title="Dependence of -0.25")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Probability")+xlab("Sample size")
Coverage=GeneralGraph(7, omega = "0.25", omega_to_file = omega_to_file, useDependence = FALSE,useAsymtotic = TRUE)
Coverage025=ggplot(Coverage,aes(x=Size , y=Coverage,color = Method))+geom_line()+labs(title="Dependence of 0.25")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Probability")+xlab("Sample size")
Coverage=GeneralGraph(7, omega = "0.75", omega_to_file = omega_to_file, useDependence = FALSE,useAsymtotic = TRUE)
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
MSE=GeneralGraph(3, omega = NULL, omega_to_file = omega_to_file, useDependence = TRUE) #just to initialize.
MSE[,2]=MSE$SD^2+GeneralGraph(5, omega = NULL, omega_to_file = omega_to_file, useDependence = TRUE)$Bias^2
colnames(MSE)=list("Size","MSE","Dependence")
MsemBeta=ggplot(MSE,aes(x=Size , y=MSE , color = Dependence))+geom_line()+labs(title=" ")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+xlab("Sample size")

# Plot for Average Length
Length=GeneralGraph(6, omega = NULL, omega_to_file = omega_to_file, useDependence = TRUE)
lengmBeta=ggplot(Length,aes(x=Size , y=Length,color = Dependence))+geom_line()+labs(title=" ")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Average length")+xlab("Sample size")

# Plot for Probability Coverage
Coverage=GeneralGraph(7, omega = NULL, omega_to_file = omega_to_file, useDependence = TRUE)
CoveragemBeta=ggplot(Coverage,aes(x=Size , y=Coverage,color = Dependence))+geom_line()+labs(title=" ")+
  scale_x_continuous(breaks=c(10,seq(100,1000,100)))+ylab("Probability")+xlab("Sample size")

grid.arrange(BiasBeta, MsemBeta, lengmBeta, CoveragemBeta, nrow = 2,ncol=2,layout_matrix = rbind(c(1, 2),
                                                                                                 c(3, 4)))













