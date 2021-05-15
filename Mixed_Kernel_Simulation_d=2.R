#The mksmld2 function is used to evaluate the performance of linear interaction model, nonlinear interaction model and mixed interaction model that are applied to test significant main and interaction effects in (generalized) linear model when there are only two main effects.
#R version: R/3.6.3.

#Description
#This function can be used to evaluate the performance of linear interaction model, nonlinear interaction model and mixed interaction model when there are only two main effects. It is a general funciton which can be conveniently used to simulate the scenarios that users set. This function is applied in our simualtion studies. We list the scenarios in our paper after the function. MASS and mvtnorm package are needed to be installed before running this function.

library(mvtnorm)
library(MASS)


#Usage
mksmld2<-function(B,n,C,tau,beta0,alpha,beta,w,family=c("gaussian","binomial"),...){
  
  #Arguments
  #B  number of iterations.
  #n  sample size.
  #C  vector of weights which control the proportion of linear and nonlinear effect when generating the data with mixing interaction effects. 
      #It is between 0 and 1. When C=1, linear interaction dominates and when C=0, nonlinear interaction dominates
      #A mixing effect is generated when C takes value between 0 and 1. 
  #tau  number of standard deviation of random error terms.
  #beta0  number of intercept of generated regression model.
  #alpha  vector of coefficients of main effects when generating the data. 
  #beta  number of coefficient of interaction effect when generating the data. 
  #w  vector of weights used to analyze the data with a mixed kernel. p-values are combined with the Cauchy comination method under different w values.
  #family   'gaussian' or 'binomial' corresponding to continuous or binary response variables.  
  # ... : other arguments passed to glm.
  
  #Value
  #powerkernel  power or size vector of interaction effect
  #poweroverall  power or size vector of three interaction models
  #poweroverkernel  list of poweroverall and powerkernel
  
  #There are two main effects, so only one interaction effect in each interaction model. Taking the main and interaction effects into account, the column number of power and p-value matrix is 3.
  d=2
  
  #test_P is the matrix of p-value obtained from B times of GLM. 
  test_Pl=test_Pg=test_Pm=matrix(0,B,3)
  #overalltest_P is the matrix of p-value obtained from B times of overall test.
  overalltest_Pl=overalltest_Pg=overalltest_Pm=matrix(0,B,1)
  #powerl, powerg and powerm are the power and size matrix.
  power1=power2=power3=powerall=matrix(0,length(C),3)
  
  for (t in 1:length(C))
  {
    c=C[t]
    for (i in 1:B)
    {
      # assuming independence between gene variables   
      Rho=0
      sigma=matrix(0,d,d)
      for (k in 1:d)
      {
        for (h in 1:d)
        {
          sigma[k,h]=Rho^(abs(k-h))
        }
      }
      diag(sigma)=1
      
      #The main effect variables are generated from a multivariate normal distribution with mean 0 and covariance sigma.
      
      x=mvrnorm(n,rep(0,d),sigma)
      x1=x[,1]
      x2=x[,2]
      
      #Using linear kernel,gaussian kernel and a kernel that is similar to guassian kernel to define the interaction effect, we combine three types of interaction term with main effects matrix to constuct three types of independent variable matrix.
      z1=scale(x1*x2) #linear kernel
      kernell<-cbind(x1,x2,z1) 
      z2=scale(exp(-0.5*((x1-x2)^2))) #Gaussian kernel
      kernelg<-cbind(x1,x2,z2)
      z.2=scale(exp(-0.5*(abs(x1-x2)))) ##The reason we construc this kernel is that we will use this kernel as a nonlinear kernel to generate response variable to avoid using the same kernel to generate and analyze the data.
      
      #generating the response variable.
      if (family=="binomial") {
        #Response variable is binary, then error does not exist.
        tau<-NULL
        Y0<-beta0+x%*%alpha+c*beta*z1+(1-c)*beta*z.2
        Y=rbinom(n,1,1/(1+exp(-Y0)))
      } else {
        #Response variable is continious.
        mu<-0
        tau<-tau
        error<-rnorm(n,mu,tau)
        Y<-beta0+x%*%alpha+c*beta*z1+(1-c)*beta*z.2+error 
      }
      
      ########Generalized Linear Models (GLM)########
        #GLM in linear interaction model
        fitlinear <- glm (Y~kernell,family=family)
        test_Pl[i,]<-summary(fitlinear)$coefficients[2:4,4]
        anoval <- anova(fitlinear,test="Chisq")
        overalltest_Pl[i,]=anoval[2,5]
          
        #GLM in nonlinear interaction model
        fitgaussian <- glm (Y~kernelg,family=family)
        test_Pg[i,]<-summary(fitgaussian)$coefficients[2:4,4]
        anovag <- anova(fitgaussian,test="Chisq")
        overalltest_Pg[i,]=anovag[2,5]
          
          
        #GLM in mixed interaction model
        pval=matrix(0,length(w),ncol(kernell))
        pvaloverall=matrix(0,length(w),1)
        for (o in 1:length(w)){
        #establishing the mixed interaction vector according to each w, we combine it  with x to constitute the predictor variables to fit glm.
        z3=w[o]*z1+(1-w[o])*z2
        kernelmix<-cbind(x1,x2,z3)
        fitmix <- glm (Y~kernelmix,family = family)
        pval[o,]=summary(fitmix)$coefficients[2:4,4]
        anovam <- anova(fitmix,test="Chisq")
        pvaloverall[o,]=anovam[2,5]
        }
      
      #using Cauchy combined method to combine p-values obtained with HDI under different w values
      CauchyP = function(p)
      {
        cct = sum(tan(pi*(0.5-p))/length(p))
        return(1-pcauchy(cct))
      }
      
      test_Pm[i,]<-as.matrix(apply(pval,2,CauchyP))
      overalltest_Pm[i,]<-as.matrix(apply(pvaloverall,2,CauchyP))
      }
    
    ########Power and size calculation########
    #calculate the power or size of x1.
    power1[t,]<-c(length(which(test_Pl[,1]<0.05)),length(which(test_Pg[,1]<0.05)),length(which(test_Pm[,1]<0.05)))/B
    #calculate the power or size of x2.
    power2[t,]<-c(length(which(test_Pl[,2]<0.05)),length(which(test_Pg[,2]<0.05)),length(which(test_Pm[,2]<0.05)))/B
    #calculate the power or size of kernel.
    power3[t,]<-c(length(which(test_Pl[,3]<0.05)),length(which(test_Pg[,3]<0.05)),length(which(test_Pm[,3]<0.05)))/B
    #calculate the power or size of model.
    powerall[t,]=c(length(which(overalltest_Pl[,1]<0.05)),length(which(overalltest_Pg[,1]<0.05)),length(which(overalltest_Pm[,1]<0.05)))/B
    colnames(power1)<-colnames(power2)<-colnames(power3)<-colnames(powerall)<-c("L","NL","M")
    }
  
  powerx1<-cbind(C,power1)
  powerx2<-cbind(C,power2)
  powerkernel<-cbind(C,power3)
  poweroverall<-cbind(C,powerall)
  poweroverkernel<-list(powerkernel=powerkernel,poweroverall=poweroverall)
  
  if (sum(alpha)==0 & beta==0) {
    return(poweroverall)
  } else {
  if (sum(alpha)==0 & beta!=0){
    return(powerkernel)
  } else { 
    return(poweroverkernel)
    }
  }
}

#Different simulation scenarios
#As in our paper,we set the iterations of all scenarios as 1000

######################################################
## scenario when the response variable is continuous##
######################################################
B=1000
C=seq(0,1,by=0.2)
tau=1
beta0=0.5
alpha=rep(0,2) # for power, change alpha=c(0.2,0.2)
beta=0 # for power, change beta=0.2
w=seq(0,1,by=0.2)
sim1_n100=mksmld2(B,n=100,C,tau,beta0,alpha,beta,w,family=c("gaussian"))
sim1_n200=mksmld2(B,n=200,C,tau,beta0,alpha,beta,w,family=c("gaussian"))

######################################################
## scenario when the response variable is binary    ##
######################################################
B=1000
C=seq(0,1,by=0.2)
beta0=0.2
alpha=rep(0,2) # for power, change alpha=c(0.43,0.43)
beta=0 # for power, change beta=0.43
w=seq(0,1,by=0.2)
sim2_n100=mksmld2(B,n=100,C,tau,beta0,alpha,beta,w,family=c("binomial"))
sim2_n200=mksmld2(B,n=200,C,tau,beta0,alpha,beta,w,family=c("binomial"))

