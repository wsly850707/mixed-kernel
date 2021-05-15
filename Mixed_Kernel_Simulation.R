#The mksml function is used to evaluate the performance of linear interaction model, nonlinear interaction model and mixed interaction model that are applied to test significant main and interaction effects in high-dimensional (generalized) linear model. 
#R version: R/3.6.3.

#Description
#This function can be used to evaluate the performance of linear interaction model, nonlinear interaction model and mixed interaction model. This simulation code is for the d>2 case implemented with the HDI procedure. It is a general funciton which can be conveniently used to simulate the scenarios that users set. This function is applied in our simualtion studies. We list the scenarios in our paper after the function. MASS, mvtnorm and HDI package is needed to be installed before running this function.

library(mvtnorm)
library(MASS)
library(hdi)

#Usage
mksml<-function(B,n,d,C,tau,beta0,alpha,beta,w,family=c("gaussian","binomial"),method=c("lasso","ridge"),...){
  
  #Arguments
  #B  number of iterations.
  #n  sample size. 
  #d  number of gene variables which give d main effects and d(d-1)/2 two-way interaction effects.
  #C  vector of weights which control the proportion of linear and nonlinear effect when generating the data with mixing interaction effects. 
      #It is between 0 and 1. When C=1, linear interaction dominates and when C=0, nonlinear interaction dominates
      #A mixing effect is generated when C takes value between 0 and 1. 
  #tau  the standard deviation of the error term.
  #beta0  intercept.
  #alpha  vector of coefficients of main effects when generating the data. 
  #beta  vector of coefficients of interaction effects when generating the data. Its length equals to d*(d-1)/2. 
  #w  vector of weights used to analyze the data with a mixed kernel. p-values are combined with the Cauchy comination method under different w values.
  #family   'gaussian' or 'binomial' corresponding to continuous or binary response variables. 
  #method   using "lasso" or "ridge" in the hdi package for the de-sparsified LASSO or ridge regression
  # ... : other arguments passed to the hdi package.
  
  #Value
  #powersize  List of power and average size matrix fitting with a linear, nonlinear and mixed interaction model, containing the following components:
     #poweravsizeL  Matrix of power and average size matrix of linear interaction model.
     #poweravsizeNL  Matrix of power and average size matrix of nonlinear interaction model.
     #poweravsizeM  Matrix of power and average size matrix of mixed interaction model
  
  
  #define the column number of interaction effect matrix according to m 
  inter_num<-choose(d,2)
  
  #pvall, pvalg and pvalm are the list of p-value.
  pvall=pvalg=pvalm=list()
  #Test_P is the matrix of p-value obtained from B times of HDI. The p-value lists are composed of them.
  test_Pl=test_Pg=test_Pm=matrix(0,length(B),(d+inter_num))
  #powerl, powerg and powerm are the power and size matrix.
  powerl=powerg=powerm=matrix(0,length(C),(d+inter_num))
  
  for (t in 1:length(C)) # for each C value, run length(B) times replicates
  {
    c=C[t]
    for (i in 1:length(B)) 
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
      set.seed(B[i]+1000)
      x=mvrnorm(n,rep(0,d),sigma)
      colnames(x)<-paste0("x",1:ncol(x))
      
      #Using linear kernel,gaussian kernel and a kernel that is similar to guassian kernel to define the interaction effect, we combine three types of interaction terms with main effects matrix to constuct three types of independent variable matrix. We apply a temporary matrix to place the interaction term obtained from each calculation in case we can only get the last calculation.
      inter_setl<-matrix(0,n,)
      for (m in 1:(ncol(x)-1)){
        for (e in (m+1):(ncol(x))){
          m1<-x[,m]
          m2<-x[,e]
          templ<-cbind(m1,m2)
          #use linear kernel to define the linear interaction
          templ1<-apply(templ,1,function(x) x[1]*x[2])
          templ1<-as.matrix(templ1)
          #name the interaction terms as xd-xe
          colnames(templ1)<-paste0(colnames(x)[m],"-",colnames(x)[e])
          #put each interaction term into the temporaty matrix 
          inter_setl<-cbind(inter_setl,templ1)
        } 
      }
      inter_setl1<-scale(inter_setl[,-1])
      kernell<-scale(cbind(x,inter_setl1))
    
      inter_setg<-matrix(0,n,)
      for (f in 1:(ncol(x)-1)){
        for (g in (f+1):(ncol(x))){
          f1<-x[,f]
          f2<-x[,g]
          tempg<-cbind(f1,f2)
          #use guassian kernel to define the nonlinear interaction
          tempg1<-apply(tempg,1,function(x) exp(-0.5*((x[1]-x[2])^2)))
          tempg1<-as.matrix(tempg1)
          #name the interaction terms as xf-xg
          colnames(tempg1)<-paste0(colnames(x)[f],"-",colnames(x)[g])
          #put each interaction term into the temporaty matrix 
          inter_setg<-cbind(inter_setg,tempg1)
        } 
      }
      inter_setg1<-scale(inter_setg[,-1])
      kernelg<-scale(cbind(x,inter_setg1))
    
      #The reason we construc this kernel is that we will use this kernel as a nonlinear kernel to generate response variable to avoid using the same kernel to generate and analyze the data.
      inter_setg.<-matrix(0,n,)
      for (j in 1:(ncol(x)-1)){
        for (l in (j+1):(ncol(x))){
          j1<-x[,j]
          j2<-x[,l]
          tempg.<-cbind(j1,j2)
          #use this kernel to define the interaction
          tempg.1<-apply(tempg.,1,function(x) exp(-0.5*(abs(x[1]-x[2]))))
          tempg.1<-as.matrix(tempg.1)
          #name the interaction terms as xj-xl
          colnames(tempg.1)<-paste0(colnames(x)[j],"-",colnames(x)[l])
          #put each interaction term into the temporaty matrix 
          inter_setg.<-cbind(inter_setg.,tempg.1)
        } 
      }
      inter_setg.1<-scale(inter_setg.[,-1])
      kernelg.<-scale(cbind(x,inter_setg.1))
      
    #generating the response variable. 
    if (family=="binomial") {
      #Response variable is binary and no error term is needed.
      tau<-NULL
      Y0<-beta0+x%*%alpha+c*(inter_setl1%*%beta)+(1-c)*(inter_setg.1%*%beta)
      Y=rbinom(n,1,1/(1+exp(-Y0)))
    } else {
      #Response variable is continuous.
      mu<-0
      tau<-tau
      error<-rnorm(n,mu,tau)
      Y<-beta0+x%*%alpha+c*(inter_setl1%*%beta)+(1-c)*(inter_setg.1%*%beta)+error 
    }

  ########High-dimensional Inference (HDI)########
    #HDI in linear interaction model
    set.seed(1000)
    if (method == "lasso") fitlinear <- lasso.proj (kernell, Y, family = family) else fitlinear <- ridge.proj(kernell, Y, family = family)
    test_Pl[i,]<-fitlinear$pval
    pvall[[t]]<-test_Pl
    
    #HDI in nonlinear interaction model
    set.seed(1000)
    if (method == "lasso") fitgaussian <- lasso.proj (kernelg, Y, family = family) else fitgaussian <- ridge.proj(kernelg, Y, family = family)
    test_Pg[i,]<-fitgaussian$pval
    pvalg[[t]]<-test_Pg
    
    #HDI in mixed interaction model
    pval=matrix(0,length(w),ncol(kernell))
    for (o in 1:length(w)){
      #establishing the mixed interaction matrix according to each w, we combine it with x to construct the predictor variables to do high-dimenional inference.
      kernelmix=w[o]*inter_setl1+(1-w[o])*inter_setg1
      kernelmix1<-cbind(scale(x),kernelmix)
      set.seed(1000)
      if (method == "lasso") fitmix <- lasso.proj (kernelmix1, Y, family = family) else fitmix <- ridge.proj(kernelmix1, Y, family = family)
      pval[o,]=fitmix$pval
    }
    
    #using Cauchy combined method to combine p-values obtained with HDI under different w values
    CauchyP = function(p)
    {
      cct = sum(tan(pi*(0.5-p))/length(p))
      return(1-pcauchy(cct))
    }
    
    test_Pm[i,]<-as.matrix(apply(pval,2,CauchyP))
    pvalm[[t]]<-test_Pm
    }
    
  ########Power and size calculation########
  #calculate the power and size of the linear interaction effect model. 
  test_Pl1<-test_Pl
  test_Pl1[test_Pl1>=0.05]=2
  test_Pl1[test_Pl1<0.05]=1
  test_Pl1[test_Pl1==2]=0
  powerl[t,]<-apply(test_Pl1,2,function(x) sum(x)/length(x))
  #calculate the power and size of the nonlinear interaction effect model.
  test_Pg1<-test_Pg
  test_Pg1[test_Pg1>=0.05]=2
  test_Pg1[test_Pg1<0.05]=1
  test_Pg1[test_Pg1==2]=0
  powerg[t,]<-apply(test_Pg1,2,function(x) sum(x)/length(x))
  #calculate the power and size of the mixed interaction effect model.
  test_Pm1<-test_Pm
  test_Pm1[test_Pm1>=0.05]=2
  test_Pm1[test_Pm1<0.05]=1
  test_Pm1[test_Pm1==2]=0
  powerm[t,]<-apply(test_Pm1,2,function(x) sum(x)/length(x))
  colnames(powerl)<-colnames(powerg)<-colnames(powerm)<-colnames(kernell)
  rownames(powerl)<-rownames(powerg)<-rownames(powerm)<-c("0","0.2","0.4","0.6","0.8","1.0")
  }
  
  #output the power matrix of the linear interaction model
  mainpowerl<-powerl[,which(alpha!=0)]
  interpowerl<-powerl[,length(alpha)+which(beta!=0)]
  #output the size of main effects and interaction effects of the linear interaction model
  mainsizel<-powerl[,which(alpha==0)]
  intersizel<-powerl[,length(alpha)+which(beta==0)]
  #calculate the average size of main effect and interaction effect
  mainavsizel<-as.matrix(apply(mainsizel,1,mean))
  colnames(mainavsizel)<-"mainsize"
  interavsizel<-as.matrix(apply(intersizel,1,mean))
  colnames(interavsizel)<-"intersize"
  #combine power matrix and average size matrix 
  poweravsizeL<-cbind(mainpowerl,mainavsizel,interpowerl,interavsizel)
  
  #output the power matrix of the nonlinear interaction model
  mainpowerg<-powerg[,which(alpha!=0)]
  interpowerg<-powerg[,length(alpha)+which(beta!=0)]
  #calcuate the size of main effects and interaction effects of the nonlinear interaction model
  mainsizeg<-powerg[,which(alpha==0)]
  intersizeg<-powerg[,length(alpha)+which(beta==0)]
  #calculate the average size of main effects and interaction effects
  mainavsizeg<-as.matrix(apply(mainsizeg,1,mean))
  colnames(mainavsizeg)<-"mainsize"
  interavsizeg<-as.matrix(apply(intersizeg,1,mean))
  colnames(interavsizeg)<-"intersize"
  #combine power matrix and average size matrix 
  poweravsizeNL<-cbind(mainpowerg,mainavsizeg,interpowerg,interavsizeg)
  
  #output the power matrix of the mixed interaction model
  mainpowerm<-powerm[,which(alpha!=0)]
  interpowerm<-powerm[,length(alpha)+which(beta!=0)]
  #calculate the size of main effects and interaction effects of the mixed interaction model
  mainsizem<-powerm[,which(alpha==0)]
  intersizem<-powerm[,length(alpha)+which(beta==0)]
  #calculate the average size of main effects and interaction effects
  mainavsizem<-as.matrix(apply(mainsizem,1,mean))
  colnames(mainavsizem)<-"mainsize"
  interavsizem<-as.matrix(apply(intersizem,1,mean))
  colnames(interavsizem)<-"intersize"
  #combine power matrix and average size matrix 
  poweravsizeM<-cbind(mainpowerm,mainavsizem,interpowerm,interavsizem)
  
  powersize<-list(poweravsizeL=poweravsizeL,poweravsizeNL=poweravsizeNL,poweravsizeM=poweravsizeM)
  
  return(powersize)
}

#Different simulation scenarios
#As in our paper,we set the iterations of all scenarios as 1000 which may take soem time to run depending on your computer speed.

######################################################
## scenario when the response variable is continuous##
######################################################
B=c(1:1000)
C=seq(0,1,by=0.2)
tau=0.8
beta0=0.5
d=10
alpha<-rep(0,d)
alpha[1:5]<-rep(0.2,5)  # the first five main effects for the power and the rest zeros for the average size. 
beta<-rep(0,d*(d-1)/2)
beta[c(2,11,16,19,27)]<-rep(0.2,5)  #non-zero coefficients for the power of interactions and the rest zeros for the average size of interactions
w=seq(0,1,by=0.2)

sim1_d10_n100=mksml(B,n=100,d,C,tau,beta0,alpha,beta,w,family=c("gaussian"),method=c("lasso"))
#The meesage shown in the screen during running the program comes from the HDI package. It is not a warning message.
#sim1_d10_n200=mksml(B,n=200,d,C,tau,beta0,alpha,beta,w,family=c("gaussian"),method=c("lasso"))

#The scenario of d=20 or 50 when response variable is continious are similar to the scenario of n=100, d=10. 
#Just need to modify d and n, the rest settings remain the same.
#when d changes to 20
d=20
alpha<-rep(0,d) 
alpha[1:5]<-rep(0.2,5)   # the first five main effects for the power and the rest zeros for the average size. 
beta<-rep(0,d*(d-1)/2)
beta[c(2,11,16,19,27)]<-rep(0.2,5)  #non-zero coefficients for the power of interactions and the rest zeros for the average size of interactions
#sim1_d20_n100=mksml(B,n=100,d,C,tau,beta0,alpha,beta,w,family=c("gaussian"),method=c("lasso"))
#sim1_d20_n200=mksml(B,n=200,d,C,tau,beta0,alpha,beta,w,family=c("gaussian"),method=c("lasso"))


##################################################
## scenario when the response variable is binary##
################################################## 
#tau=NULL
d=10
beta0=0.05
alpha=alpha<-rep(0,d)
alpha[1:5]<-rep(0.5,5)  # the first five main effects for the power and the rest zeros for the average size. 
beta<-rep(0,d*(d-1)/2)
beta[c(2,11,16,19,27)]<-rep(0.5,5)  #non-zero coefficients for the power of interactions and the rest zeros for the average size of interactions
w=seq(0,1,by=0.2)
sim2_d10_n100=mksml(B,n=100,d,C,tau,beta0,alpha,beta,w,family=c("binomial"),method=c("lasso"))
#sim2_d10_n200=mksml(B,n=200,d,C,tau,beta0,alpha,beta,w,family=c("binomial"),method=c("lasso"))

#The scenario of d=20 or 50 when response variable is binary are similar to the scenario of n=100, d=10. 
#Just need to modify d and n, the rest settings remain the same.
#when d changes to 20
d=20
alpha=alpha<-rep(0,d)
alpha[1:5]<-rep(0.5,5)  # the first five main effects for the power and the rest zeros for the average size. 
beta<-rep(0,d*(d-1)/2)
beta[c(2,11,16,19,27)]<-rep(0.5,5)  #non-zero coefficients for the power of interactions and the rest zeros for the average size of interactions
sim2_d20_n100=mksml(B,n=100,d,C,tau,beta0,alpha,beta,w,family=c("binomial"),method=c("lasso"))
#sim2_d20_n200=mksml(B,n=200,d,C,tau,beta0,alpha,beta,w,family=c("binomial"),method=c("lasso"))

