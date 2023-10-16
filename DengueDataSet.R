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

###############################################################################
#Dengue Data
###############################################################################

#Import Data from Excel
Dengue <- read.xlsx("DengueDataDx.xlsx", sheetIndex = "dengue", header = T)
Dengue$Dx=as.factor(Dengue$dx)
levels(Dengue$Dx)=c("Negative","Positive")
Y=Dengue[,c(1,2,4)]

######################################################################
# Diagrams of Leukocytes and Platelets variables.
######################################################################
#Histogram of leukocytes with 1379 values
ggplot(Y, aes(x=leuco_totales, fill=Dx)) + 
  geom_histogram(aes(y=..density..), position="identity", alpha=0.5, color="black")+
  geom_density(alpha=0.2)+
  labs(title=" ",x="Leukocytes", y = "Density")

#Boxplot of leukocytes with 1379 values
ggplot(Y, aes(x=1,y=leuco_totales, fill=Dx))+
  geom_boxplot()+ labs(title=" ",x="Dx", y = "Leukocytes")


#Histogram of platelets with 1379 values
ggplot(Y, aes(x=plaquetas, fill=Dx)) + 
  geom_histogram(aes(y=..density..), position="identity", alpha=0.5, color="black")+
  geom_density(alpha=0.2)+
  labs(title=" ",x="Platelets", y = "Density")

#Boxplot of platelets with 1379 values
ggplot(Y, aes(x=1,y=plaquetas, fill=Dx)) + 
  geom_boxplot()+
  labs(title=" ",x="Dx", y = "Platelets")

######################################################################
# Correlation between random variable of Y.
######################################################################
cor(Y[,c(1,2)],use = "everything", method = "kendall")
cor(Y[which(Y$Dx=="Negative"),c(1,2)],use = "everything", method = "kendall")
cor(Y[which(Y$Dx=="Positive"),c(1,2)],use = "everything", method = "kendall")

cor(Y[,c(1,2)],use = "everything", method = "spearman")
cor(Y[which(Y$Dx=="Negative"),c(1,2)],use = "everything", method = "spearman")
cor(Y[which(Y$Dx=="Positive"),c(1,2)],use = "everything", method = "spearman")


cor.test(Y[which(Y$Dx=="Positive"),1],Y[which(Y$Dx=="Positive"),2],method=c("kendall"))
cor.test(Y[which(Y$Dx=="Positive"),1],Y[which(Y$Dx=="Positive"),2],method=c("spearman"))

# Pearson Correlation Summary:
# t-statistic: 1.4075, df: 742, p-value: 0.1597
# The p-value is greater than the 0.05 significance level, indicating insufficient evidence to reject the null hypothesis that the true correlation is zero.
# 95% confidence interval includes zero (-0.0204 to 0.1230), further confirming non-significance.
# Sample correlation estimate: 0.0516, suggesting a weak, non-significant linear relationship between the two variables for the "Positive" diagnosis group.
cor.test(Y[which(Y$Dx=="Positive"),1],Y[which(Y$Dx=="Positive"),2],method=c("pearson"))

######################################################################
# Dot plot for cumulative probabilities.
######################################################################
# Marginal distribution by leukocyte count and platelet count
U1=Marginal(1,Y)
U2=Marginal(2,Y)
U=data.frame(matrix(data = c(U1,U2),nrow=1380, ncol=2, byrow = F))
U$X3=as.factor(Y$Dx)
names(U)=c("X1","X2","Dx")
levels(U$Dx)=c("Non Dengue Case","Dengue Case")

#Scatter plot of positive cases
ggplot(data = U, aes(x=X1,y = X2, color=Dx)) + geom_point()+
  xlab(expression(u[1]))+ylab(expression(u[2]))+labs(title = " ")
#names(U)=c("X1","X2","X3")



######################################################################
# Results obtained from point estimates and intervals.
######################################################################
Est_Dengue=Est_A_Samp(dataset_target = U[which(U$Dx=="Dengue Case"),c(1,2)],nboot = 200,confidence = 0.95,
                      a=186.875, b=112.125, c = 0.25,n_burnd = 100,
                      n_iter = 10000,n_chain = 1,n_thin = 100)

# Certainly. When one evaluates model performance using the Deviance Information Criterion (DIC), 
# a lower DIC value is generally preferred. This suggests that the model with the lowest DIC value
# provides a better fit to the observed data while penalizing for model complexity. 
# In this specific case, the model with a DIC of -36.27 would likely be the best choice 
# for representing the data, as it has the lowest DIC value.

round(
  matrix(data = c(Est_Dengue$MLE, Est_Dengue$EM$varphikendall,Est_Dengue$EM$varphispearman,
                Est_Dengue$EBT,Est_Dengue$EBB,Est_Dengue$EBU, var(as.numeric(Est_Dengue$vphilm)),
                var(as.numeric(Est_Dengue$vphitau)),var(as.numeric(Est_Dengue$vphispe)),Est_Dengue$summ_Trian[1,2]^2,
                Est_Dengue$summ_Beta[1,2]^2,Est_Dengue$summ_Unif[1,2]^2,
                Est_Dengue$IML[[1]],Est_Dengue$ITau[[1]],Est_Dengue$ISpe[[1]],
                Est_Dengue$CRT[[1]],Est_Dengue$CRB[[1]],Est_Dengue$CRU[[1]],
                Est_Dengue$IML[[2]],Est_Dengue$ITau[[2]],Est_Dengue$ISpe[[2]],
                Est_Dengue$CRT[[2]],Est_Dengue$CRB[[2]],Est_Dengue$CRU[[2]],
                Est_Dengue$IML[[2]]-Est_Dengue$IML[[1]],Est_Dengue$ITau[[2]]-Est_Dengue$ITau[[1]],
                Est_Dengue$ISpe[[2]]-Est_Dengue$ISpe[[1]],Est_Dengue$CRT[[2]]-Est_Dengue$CRT[[1]],
                Est_Dengue$CRB[[2]]-Est_Dengue$CRB[[1]],Est_Dengue$CRU[[2]]-Est_Dengue$CRU[[1]],
                0,0,0,Est_Dengue$DIC_T,Est_Dengue$DIC_B,Est_Dengue$DIC_U),ncol = 6,nrow = 6,
       dimnames = list(list("MLE","Kendall","Spearman","Triangular","Beta","Uniform"),
                       list("Estimate","Variance","Lower","Upper","Range","DIC"))),
  digits = 3)


######################################################################
# Posterior density plot
######################################################################

PlotDens=data.frame(
  Posterior=factor(rep(c("Triangular", "Beta","Uniform"), each=9900)),
  varphi=c(Est_Dengue$Sims_T,Est_Dengue$Sims_B,Est_Dengue$Sims_U))

ggplot(PlotDens,aes(x=varphi,fill=Posterior))+geom_density(alpha=0.4)+ 
  labs(title="Posterior Density",x=expression(phi), y = "Density")

