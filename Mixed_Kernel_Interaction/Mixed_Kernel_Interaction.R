#The mk function is used to test significant main and interaction effects in high-dimensional (generalized) linear model.
#R version: R/3.6.3.

#Description
#This function can be used to perform inference to test significant main and interaction effects in high-dimensional (generalized) linear model. HDI package is needed to be installed before running this function. 

library(hdi)

#Usage
mk<- function (X, Y, cov, w, family = c("gaussian","binomial"), method = c("lasso", "ridge"), ...){	
  
  #Arguments
  #X   main effect matrix that can be gene expression data.
  #Y   response vector that can be either continuous or binary.
  #cov   data.frame or matrix of covariates. 
  #w   vector of weights used to analyze the data with a mixed kernel. p-values are combined with the Cauchy comination method under different w values. 
  #family   either 'gaussian' or 'binomial', relying on the type of response. 
  #method   either "lasso" or "ridge" to estimate the main and interaction effects.
  # ... : other arguments passed to hdi.
  
  #Value 
  #Caup   Vetcor of p-values that are cauchy combination of p-values obtained under different weights.
  
  ########Step-1  define linear and nonlinear interaction#########
  #Choosing the ith and jth X variables, we use a linear kernel to define the linear interaction b/w them. 
  temp_totall<-matrix(0,nrow(X),)
  for (i in 1:(ncol(X)-1)){
    for (j in (i+1):(ncol(X))){
      a1<- X[,i]
      a2<- X[,j]
      templ<- cbind(a1,a2)
      #define linear interaction     
      templ1<- apply(templ,1,function(x) x[1]*x[2])
      templ1<- as.matrix(templ1)
      #name the interaction terms as Xi_Xj
      colnames(templ1)<- paste0(colnames(X)[i],"_",colnames(X)[j])
      temp_totall<- cbind(temp_totall,templ1)
    } 
  }
  #scale the variables
  kernell<- scale(temp_totall[,-1])

  #Choosing the ith and jth X variables, we use a Gaussian kernel to define the non-linear interaction b/w them.
  temp_totalnl<- matrix(0,nrow(X),)
  for (i in 1:(ncol(X)-1)){
    for (j in (i+1):(ncol(X))){
      a1<- X[,i]
      a2<- X[,j]
      tempnl<- cbind(a1,a2)
      #use a Gaussian kernel to define the non-linear interaction
      tempnl1<- apply(tempnl,1,function(x) exp(-0.5*((x[1]-x[2])^2)))
      tempnl1<- as.matrix(tempnl1)
      #name the interaction terms as Xi_Xj
      colnames(tempnl1)<- paste0(colnames(X)[i],"_",colnames(X)[j])
      temp_totalnl<- cbind(temp_totalnl,tempnl1)
    } 
  }
  kernelnl<- scale(temp_totalnl[,-1])
  
  ########Step-2   High-dimensional Inference (HDI)#########
  #Based on the linear and nonlinear interaction matrices, we constitute the mixed interaction matrix under different w.
    if (is.null(cov)) {
    #construct a p-value matrix containing p-values of X and its interaction terms      
      pval=matrix(0, length(w), ncol(X)+ncol(kernell)) 
      pvaladj=matrix(0, length(w), ncol(X)+ncol(kernell)) 
      for (k in 1: length(w)){
      #establishing the mixed interaction matrix according to each w.
      kernelm = w[k]* kernell+ (1-w[k])* kernelnl
      Eff <- cbind (scale(X), kernelm)
      #HDI
      set.seed(2000)
      if (method == "lasso") fit <- lasso.proj (Eff, Y, family = family,standardize=FALSE) else fit <- ridge.proj(Eff, Y, family = family,standardize=FALSE)
      pval[k,]<- fit$pval
      pvaladj[k,]<-pval[k,]*mta(Eff)
      }
      #Name the column of p-value matrix obtained from HDI
      colnames(pvaladj)<-c(colnames(X),colnames(kernell))
      } else {
    #If there are covariates in the model, then main effects, interaction effects and covariates are all fitted.
    pval=matrix(0, length(w), ncol(cov)+ncol(X)+ncol(kernell))
    pvaladj=matrix(0, length(w), ncol(cov)+ncol(X)+ncol(kernell))
    for (k in 1 :length(w)){ 
      kernelm = w[k]* kernell+ (1-w[k])* kernelnl 
      Eff<- cbind (scale(cov), scale(X), kernelm)
      Eff1<-cbind(scale(X), kernelm)
      set.seed(2000)
      if (method == "lasso") fit <- lasso.proj(Eff, Y, family = family,standardize=FALSE) else fit <- ridge.proj(Eff, Y, family = family,standardize=FALSE)
      pval[k,]<- fit$pval
      pvaladj[k,]<-pval[k,]*mta(Eff1)
    }
    #Name the column of p-value matrix obtained from HDI
    colnames(pvaladj)<-c(colnames(cov),colnames(X),colnames(kernell))
  }
  #There are part of p-values larger than 1 after multiple testing adjustment, we let them be 1. 
  pvaladj[which(pvaladj>=1)]=1
  
  
  ########STEP-3   P-value combination#########
  #P-value obtained from HDI according to different w need to be combined for further inference.
  Caup<- as.matrix(apply(pvaladj,2,delp))
  return(Caup)
}

#Multiple testing adjustment

mta<-function(x)
{
  corx<-cor(x)
  ecl<-eigen(corx)
  ecl1<-ecl$values
  ecl11<-abs(ecl1)
  ecl111<-ifelse(ecl11>=1,1,0)
  ecl112<-ecl11-floor(ecl11)
  fx<-ecl111+ecl112
  sum.fx<-sum(fx)
  return(sum.fx)
}

#Cauchy combination
CauchyP = function(p)
{
  cct = sum(tan(pi*(0.5-p))/length(p))
  return(1-pcauchy(cct))
}
#If w p-values are all equal to 1 in a variable, then they are not combined by Cauchy. Their combination is 1. 
#If w p-values are all less than 1, then they are combined by Cauchy. If part of p-values are less then 1 and the others are equal to 1, then we just use Cauchy to combine the ones less than 1.
delp<-function(x){
  if (length(x)==1){
    out<-x
  } else {
    index1<-which(x<1)
    index2<-which(x==1)
    if(length(index1)==length(x)){
      out<-CauchyP(x)
    } else if (length(index2)==length(x)){
      out<-CauchyP(x)
    } else {
      x1<-c(x[-index2])
      out<-CauchyP(x1)
    }
  }
  return(out)
}

#Examples
#The gene expression data of lung cancer was downloaded from the [GEO] repository https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4115 
  #corresponding to the case study 2 in the paper.
#Cancer status, clinical characters and gene expression are included in the data. We separated the dataset into three files named ending, charac and 
  #expression respectively.

#####input the GEO data#####
library(data.table)
ending<-fread("ending.txt",header=T,sep="\t")
#delete the suspicious cancer samples
ending1<-as.matrix(ending[-c(188:192),])
#According to the status, we replace the status with 0 and 1. 
ending1[,2][which(ending1[,2]=="no cancer")]<-0
ending1[,2][which(ending1[,2]=="cancer")]<-1
#When loading ending, the rowname of it is regarded as first column, so we adjust it.
ending2<-as.matrix(ending1[,-1])
rownames(ending2)<-ending1[,1]
expression<-fread("expression.txt",header=T,sep="\t")
#The gene expression data is obtained based on probeset, so we need to map probesets to genes. 
  #The "gene and probeID" file can also be get from the GEO website.
#Only gene and probeID are picked from the file.
pro<-fread("gene and probeID.txt",header=T,sep="\t")
proexp<-cbind(pro[,2],expression[,-1])
#For genes with different probesets , we calculated the average of different probesets as the gene level expression. 
proexp1<-aggregate(x=proexp[,-1],by=list(proexp$GeneSymbol),FUN=mean)
#input clinical characters data. Seven clinical characters containing age, gender, race, smoking status, pack years of smoking, hemopytsis and lymphadenopathy were picked up into analysis.
pheno<-fread("pheno.txt",header=T,sep="\t")
pheno1<-as.matrix(pheno[,-1])
rownames(pheno1)<-as.matrix(pheno[,1])

#####multivariate logistic regression to select signficant covariates to be included in the model #####
#samples with missing values were deleted.
phenona<-is.na(pheno1[,1])
indexna<-which(phenona==TRUE)
pheno2<-pheno1[-indexna,]
state1<-as.numeric(ending2[,1])[-indexna]
age1<-as.numeric(pheno2[,1])
packyear1<-as.numeric(pheno2[,5])
logistmulti<-glm(state1~age1+pheno2[,2]+pheno2[,3]+pheno2[,4]+packyear1+pheno2[,6]+pheno2[,7],family="binomial")
summary(logistmulti)
#Only age and lymphadenopathy are remained as meaningful covariates. 
lymphadenopathy<-pheno2[,7]
xx<-data.frame(age1,lymphadenopathy)
xxx<-model.matrix(~lymphadenopathy,xx)
cov<-cbind(age1,xxx[,2])

#####mapping the KEGG pathway genes and GEO genes#####
genename<-fread("hsa05222.txt",header=T,sep="\t") #other pathways can also be changed here to do the following analysis.
#mapping genes to the KEGG pathway.
genename<-as.matrix(genename)
int<-intersect(proexp1[,1],genename[,2])
rownames(proexp1)<-proexp1[,1]
proexp2<-proexp1[,-1]
#there are 64 genes mapped to the pathway.
expressiongene<-proexp2[int,]
#The rows represent genes which need to transposed. 
expressiongene1<-t(expressiongene[,-c(188:192)])
#samples with missing values were deleted.
expressiongene2<-expressiongene1[-indexna,]

w=seq(0,1,by=0.2)

hsa05222<-mk(expressiongene2,state1,cov,w,family="binomial",method="lasso")
