###################################################################################################
# Load necessary functions.
###################################################################################################
source("FunctionsRequired.R")

###################################################################################################
# Dengue Data
###################################################################################################

# Import Data from Excel
Dengue = read.xlsx("DengueDataDx.xlsx", sheetIndex = "dengue", header = T)
Dengue$Dx = as.factor(Dengue$dx)
levels(Dengue$Dx) = c("Negative","Positive")
Y = Dengue[,c(1,2,4)]

###################################################################################################
# Diagrams of Leukocytes and Platelets variables.
###################################################################################################
# Histogram of leukocytes with 1379 values
ggplot(Y, aes(x=leuco_totales, fill=Dx)) + 
  geom_histogram(aes(y=..density..), position="identity", alpha=0.5, color="black")+
  geom_density(alpha=0.2)+
  labs(title=" ",x="Leukocytes", y = "Density")

# Boxplot of leukocytes with 1379 values
ggplot(Y, aes(x=1,y=leuco_totales, fill=Dx))+
  geom_boxplot()+ labs(title=" ", x="Dx", y = "Leukocytes")

#Histogram of platelets with 1379 values
ggplot(Y, aes(x=plaquetas, fill=Dx)) + 
  geom_histogram(aes(y=..density..), position="identity", alpha=0.5, color="black")+
  geom_density(alpha=0.2)+
  labs(title=" ",x="Platelets", y = "Density")

#Boxplot of platelets with 1379 values
ggplot(Y, aes(x=1,y=plaquetas, fill=Dx)) + 
  geom_boxplot()+
  labs(title=" ",x="Dx", y = "Platelets")

###################################################################################################
# Correlation between random variable of Y (leukocytes, platelets)
###################################################################################################
cor(Y[,c(1,2)],use = "everything", method = "kendall") # All patients.
cor(Y[which(Y$Dx == "Negative"),c(1,2)],use = "everything", method = "kendall") # With negative diagnosis.
cor(Y[which(Y$Dx == "Positive"),c(1,2)],use = "everything", method = "kendall") # With positive diagnosis.

cor(Y[,c(1,2)],use = "everything", method = "spearman") # All patients.
cor(Y[which(Y$Dx == "Negative"),c(1,2)],use = "everything", method = "spearman") # With negative diagnosis.
cor(Y[which(Y$Dx == "Positive"),c(1,2)],use = "everything", method = "spearman") # With positive diagnosis.

# Hypothesis test for correlation coefficient.
cor.test(Y[which(Y$Dx == "Positive"),1],Y[which(Y$Dx == "Positive"),2],method = c("kendall"))
cor.test(Y[which(Y$Dx == "Positive"),1],Y[which(Y$Dx == "Positive"),2],method = c("spearman"))

# Pearson Correlation Summary:
# t-statistic: 1.4075, df: 742, p-value: 0.1597
# The p-value is greater than the 0.05 significance level, indicating insufficient evidence to reject the null hypothesis that the true correlation is zero.
# 95% confidence interval includes zero (-0.0204 to 0.1230), further confirming non-significance.
# Sample correlation estimate: 0.0516, suggesting a weak, non-significant linear relationship between the two variables for the "Positive" diagnosis group.
cor.test(Y[which(Y$Dx == "Positive"),1],Y[which(Y$Dx == "Positive"),2],method = c("pearson"))

###################################################################################################
# Dot plot for cumulative probabilities.
###################################################################################################
# Marginal distribution by leukocyte count and platelet count
U1 = Marginal(1,Y)
U2 = Marginal(2,Y)
U = data.frame(matrix(data = c(U1,U2),nrow = 1380, ncol = 2, byrow = F))
U$X3 = as.factor(Y$Dx)
names(U) = c("X1","X2","Dx")
levels(U$Dx) = c("Non Dengue Case","Dengue Case")

# Scatter plot of positive cases
ggplot(data = U, aes(x = X1,y = X2, color = Dx)) + geom_point()+
  xlab(expression(u[1]))+ylab(expression(u[2]))+labs(title = " ")
#names(U)=c("X1","X2","X3")

###################################################################################################
# Obtaining quantile interval for the dependency parameter varphi
###################################################################################################
n = dim(U[which(U$Dx == "Dengue Case"), c(1, 2)])[1]
Boot_Varphi_Samples = matrix(data = NA, nrow = n, ncol = 0)
Boot_Tau = NA

# Perform 200 resamples of size 744 and calculate Kendall's Tau for each resample.
for (i in 1:100) {
  Boot_Varphi_Samples = cbind(Boot_Varphi_Samples, U[which(U$Dx == "Dengue Case"), c(1, 2)][sample(1:n, size = n, replace = TRUE), c(1, 2)])
  Boot_Tau[i] = cor(Boot_Varphi_Samples[, c((2 * i - 1), (2 * i))], use = "everything", method = "kendall")[1, 2]
}

# Density plot for the bootstrap estimator of Kendall's Tau
boothistration = ggplot(as.data.frame(Boot_Tau), aes(x = Boot_Tau)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  geom_density(lwd = 1.2,
               linetype = 2,
               colour = 2, fill = 4, alpha = 0.25) +
  labs(title = " ") +
  ylab("Density") +
  xlab(substitute(va, list(va = "Kendall's Tau")))
boothistration

# Select quantiles associated with 2.5% and 97.5% for Kendall's Tau
quantile_Tau = round(quantile(Boot_Tau, probs = c(0.025, 0.975)),4)

# Quantile interval for Varphi
round(9 / 2 * quantile_Tau,4)

###################################################################################################
# Obtaining hyperparameter values
###################################################################################################
Hyp_Tov = Mtovar(round(9 / 2 * quantile_Tau[1],4),round(9 / 2 * quantile_Tau[2],4),0.4)
Hyp_Tov

###################################################################################################
# Results obtained from point estimates and intervals.
###################################################################################################
data_dengue = U[which(U$Dx == "Dengue Case"),c(1,2)]
colnames(data_dengue) = c("u","v")
Est_Dengue = Est_A_Samp(dataset_target = data_dengue, nboot = 100,confidence = 0.95,
                      alpha_prior = Hyp_Tov$a, beta_prior =Hyp_Tov$b, c_prior = Hyp_Tov$c,burn_in = 5000,
                      n_iter = 30000,thinning = 5,type_sample = "exclude",precision=20)
Est_Dengue$Summary

round(
  matrix(data = c(Est_Dengue$MLE, Est_Dengue$EM$varphikendall,Est_Dengue$EM$varphispearman,Est_Dengue$Blomqvist_Beta,
                  Est_Dengue$EBT,Est_Dengue$EBB,Est_Dengue$EBU,var(as.numeric(Est_Dengue$vphilm)),
                  var(as.numeric(Est_Dengue$vphitau)),var(as.numeric(Est_Dengue$vphispe)),var(as.numeric(Est_Dengue$vphibb)),var(Est_Dengue$Sims_T),
                  var(Est_Dengue$Sims_B),var(Est_Dengue$Sims_U),
                  Est_Dengue$IML[[1]],Est_Dengue$ITau[[1]],Est_Dengue$ISpe[[1]],Est_Dengue$IBB[[1]],
                  Est_Dengue$CRT[[1]],Est_Dengue$CRB[[1]],Est_Dengue$CRU[[1]],
                  Est_Dengue$IML[[2]],Est_Dengue$ITau[[2]],Est_Dengue$ISpe[[2]],Est_Dengue$IBB[[2]],
                  Est_Dengue$CRT[[2]],Est_Dengue$CRB[[2]],Est_Dengue$CRU[[2]],
                  Est_Dengue$IML[[2]]-Est_Dengue$IML[[1]],Est_Dengue$ITau[[2]]-Est_Dengue$ITau[[1]],
                  Est_Dengue$ISpe[[2]]-Est_Dengue$ISpe[[1]],Est_Dengue$IBB[[2]]-Est_Dengue$IBB[[1]],Est_Dengue$CRT[[2]]-Est_Dengue$CRT[[1]],
                  Est_Dengue$CRB[[2]]-Est_Dengue$CRB[[1]],Est_Dengue$CRU[[2]]-Est_Dengue$CRU[[1]]),ncol = 5,nrow = 7,
         dimnames = list(list("MLE","Kendall","Spearman","Blomqvist's Beta","Triangular","Beta","Uniform"),
                         list("Estimate","Variance","Lower","Upper","Range"))),
  digits = 3)


###################################################################################################
# Posterior density plot
###################################################################################################
length(Est_Dengue$Sims_T)

PlotDens=data.frame(
  Posterior=factor(rep(c("Triangular", "Beta","Uniform"), each=length(Est_Dengue$Sims_T))),
  varphi=c(Est_Dengue$Sims_T,Est_Dengue$Sims_B,Est_Dengue$Sims_U))

ggplot(PlotDens,aes(x=varphi,fill=Posterior))+geom_density(alpha=0.4)+ 
  labs(title="",x=expression(phi), y = "Density")

###################################################################################################
# Goodness-of-fit test
###################################################################################################

# Joint marginal distribution (only for individuals with "Dengue Case")
U_DengueCase = U[which(U$Dx=="Dengue Case"),c(1,2)]*(744+1)/744

# Empirical joint distribution without parallelization
start_time = Sys.time()
U_DengueCase_Joint = Joint_Dist(U_DengueCase,U_DengueCase) 
finish_time = Sys.time()
finish_time - start_time

# Empirical joint distribution with parallelization
start_time = Sys.time()
U_DengueCase_Joint_P = Joint_Dist_Parallel(U_DengueCase,U_DengueCase) 
finish_time = Sys.time()
finish_time - start_time

# Grid of points in (0,1)x(0,1) used to evaluate the test statistic
grid_sample = expand.grid(x = seq(from = 0.01, to = 0.99, length.out = 30), y = seq(from = 0.01, to = 0.99, length.out = 30))

#####
# Step 1: Compute the goodness-of-fit test statistic using the observed data sample
#####
Est_DengueCase = me(U_DengueCase)$varphikendall # Estimate φ using Kendall's tau
Delta = abs(Joint_Dist_Parallel(U_DengueCase,grid_sample)-Density_FGM(grid_sample$x ,grid_sample$y, varphi = Est_DengueCase))
Delta0 = max(Delta) # Supremum of the difference

#####
# Step 2: Generate FGM samples of size 744 using the estimated φ from Kendall's tau
#####
N = 1000  # Number of bootstrap replications
samples_bootstrap = rnfgm(Est_DengueCase,nrow(U_DengueCase))
for (i in 1 : (N-1)) {
  samples_bootstrap = cbind(samples_bootstrap,rnfgm(Est_DengueCase,nrow(U_DengueCase)))
}

# For each generated sample, compute the φ estimate using Kendall's tau
Est_Samples_DengueCase = NA
for (i in 1 : N) {
  Est_Samples_DengueCase[i] = cor(samples_bootstrap[,c(2*i-1,2*i)], method = "kendall")[1,2]*9/2
}
Est_Samples_DengueCase

#####
# Step 3: Compute the test statistic for each generated sample
#####

start_time = Sys.time()
emp_dist = Joint_Dist_Parallel(samples_bootstrap[, c(1,2)], grid_sample)
for (i in 2:length(Est_Samples_DengueCase)) {
  idx_cols = c(2 * i - 1, 2 * i)
  emp_dist = cbind(emp_dist,Joint_Dist_Parallel(samples_bootstrap[, idx_cols], grid_sample))
  print(i)
}
finish_time = Sys.time()
finish_time - start_time

# Store the empirical copula results for each simulated dataset using the 30x30 grid generated above
write.xlsx(x = emp_dist, file = "gof_test_results.xlsx",sheetName = "Empirical_Data",col.names = TRUE,
           row.names = TRUE,append = TRUE)

# Read stored results
emp_dist = read_excel("gof_test_results.xlsx", sheet = "Empirical_Data")
emp_dist = emp_dist[,2:1001]
#ncol(emp_dist)

# Compute the goodness-of-fit statistic for all bootstrap samples
start_time = Sys.time()
numCores = detectCores() - 1 # Detect available cores
cl = makeCluster(numCores)

# Export required variables and functions to the cluster
clusterExport(cl, varlist = c("samples_bootstrap", "grid_sample", 
                              "Est_Samples_DengueCase", "Density_FGM","emp_dist"), envir = environment())

# Parallel computation of test statistics (supremum norm)
Delta_samples = parSapply(cl, 1:length(Est_Samples_DengueCase), function(i) {
  idx_cols = c(2 * i - 1, 2 * i)
  emp_distM = emp_dist[,i]
  theo_dist = Density_FGM(grid_sample$x, grid_sample$y, varphi = Est_Samples_DengueCase[i])
  max(abs(emp_distM - theo_dist))
})
stopCluster(cl)
finish_time = Sys.time()
finish_time - start_time

#####
# Quantiles for the bootstrap sample of the Delta statistic
#####

# Delta_samples  # (Stored previously)

# p-value computed using the bootstrap sample; does not consider an explicit alpha level.
sum(Delta_samples > Delta0 ) / N # p-valor

# The null hypothesis is not rejected because the test statistic Delta0 
# is greater than the 90th, 95th, and 99th percentiles of the bootstrap sample.
Delta0 > quantile(Delta_samples, probs = c(0.9,0.95,0.99))

# Histogram with density overlay
ggplot(as.data.frame(Delta_samples), aes(x = as.data.frame(Delta_samples)$Delta_samples)) + 
  geom_histogram(aes(y = after_stat(density)), colour = 1, fill = "white", binwidth = 0.02) +
  geom_density(lwd = 1.2, linetype = 2, colour = 2, fill = 4, alpha = 0.25) +
  labs(title = " ") + ylab("Density") +
  xlab(expression(Delta~"*"^(b))) +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 15),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(10, 10, 10, 10)
  ) + geom_vline(aes(xintercept = 1.3), colour = "red", linetype = 1,lwd=1.2) +
  annotate("text", x = 1.28, y = 4, label = "1.3", vjust = 2, colour = "red",size=5) +
  geom_vline(aes(xintercept = 1.483), colour = "blue", linetype = 1,lwd=1.2) +
  annotate("text", x = 1.51, y = 4, label = "1.483", vjust = 2, colour = "blue",size=5) 
