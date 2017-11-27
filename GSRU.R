#Gaussian Seidel Residual Updater - Multiple Markers
#Rversion

source("~/Documents/OneDrive - University of Edinburgh/MME_Iterators/GSRU/GSRU_Functions.R")

###########################################################

#### OPTIONS ####

marker_class = "RANDOM" #markers can be treated as fixed or random
genomic_prediction = "NO" #if YES, then validation genotypes read in and GEBVs created for validation set
lamda1=0.1 #shrinkage value applied to model when SNPs fitted as random effects. Not used if marker_class is set to FIXED
CONV=1 #set starting convergence value, 1 is suggested
CONV_TOL=0.00000000000001 #Set convergence threshold, 10^-9 is suggested

###########################################################

### DATA ENTRY ###

###########################################################

#change to directory holding your data files
setwd("~/Documents/OneDrive - University of Edinburgh/MME_Iterators/GSRU/Data_Files/MultipleMarkers")

#read in individuals (Data/No Data Score)
x<-read.csv("X_JH.csv",sep=",",header=F)
x<-as.matrix(x)
class(x)<-"numeric"
ndata=nrow(x)

#read in genotypes
z<-read.table("Genmult_JH.csv",sep=",",header=F)
z<-as.matrix(z)
class(z)<-"numeric"

#read in residuals (required to genertate differing phenotypes)
e<-read.table("E_JH.csv",sep=",",header=F)
e<-as.matrix(e)
class(e)<-"numeric"
average_e<-sum(e)/ndata

#generate phenotypes
allele_effect=read.table("Allele_Effects_JH.csv",sep=",",header=F); allele_effect = as.matrix(allele_effect)
#allele_effect=c(0.5,-1,1,-0.5)
phenotype=rep(0,ndata)
for (i in 1:ndata){
  phenotype[i] = (sum(z[i,]*allele_effect)) + e[i]
}

#############################################################################################################

#calculate allele frequencies & variances
allele_freq_stats=calc_allele_freq_per_locus(z,allele_effect,phenotype)

### MATRIX MANIPULATION ###

#assign phenotypes to y matrix
y=phenotype

#xpy<-as.matrix(t_mult_dupl(x,y))
#zpy<-as.matrix(t_mult_dupl(z,y))

#combine fixed effects and random effects into a single matrix
X=cbind(x,z)

#X<-z
if (marker_class=="FIXED"){
  #Count the number of expected solutions
  neq=ncol(X)
  #Count the number of expected of fixed effect means
  nmeans=ncol(X)
  ######################### SNP EFFECTS = FIXED ##############################
}
if (marker_class=="RANDOM") {
  #Count the number of expected solutions
  neq=ncol(X)
  #Count the number of expected of fixed effect means
  nmeans=ncol(x)

  ######################### SNP EFFECTS = RANDOM ##############################
}

#provide phenotypes as the starting values for residual updater
e=as.matrix(y)

#provide empty matrices to store Sum Squares for each Iteration
ssqd_on=as.matrix(rep(0,neq))
ssq_n=as.matrix(rep(0,neq))
#create empty "X prime X" matrix
xpx=rep(0,neq)

start=1
#for each effect in the model
for (j in 1:neq){
  #Update effect counts & generate "X prime X"
  xpx[j]<-as.matrix(t_mult_dupl(X[,j],X[,j]))  
}

#create intital empty vector to store newly generated solutions at each iteration
sol=rep(0,neq)
#create intial empty vector to store previously generated solutions at each iteration
old_sol=rep(0,neq)

#######################################################################################

#################################

####  Guass-Seidel Iterator  ####

#################################

#######################################################################################

#set counter for iterator
count=1
#until convergence is reached, iterate process
repeat{
  #for each effect in the model (Fixed + Random)
  for (j in start:neq){
    #If iterating over fixed effects
    if (j==nmeans){
      #form lhs
      lhs=xpx[j]
    }
    #If iterating over random effects
    else{
      #form lhs by adding shrinkage value to diagonal of "X prime X"
      lhs=xpx[j]+lamda1
    }
    #form rhs with y corrected by other effects (formula 1)
    rhs=(X[,j]%*%e) + (xpx[j]*sol[j])
    #do Gauss Seidel
    val=rhs/lhs
    #MCMC sample solution from its conditional (commented out here)
    #val_num=(rhs/lhs)
    #val_den=(1/lhs)
    #val=rnorm(1,(rhs/lhs),(1/lhs))
    #update e with current estimate (formula 2)
    e = e - (X[,j]*(val-sol[j]))
    #store old solution
    old_sol[j]=sol[j]
    #update sol
    sol[j]=val
    #Calcualte convergence value using the SumSquare of old vs new solutions
    ssqd_on[j]<-ssqd(old_sol[j],sol[j])
    ssq_n[j]<-sum(sol[j]^2)
    CONV<-ssqd_on[j]/ssq_n[j]
  }
  if (CONV<=CONV_TOL) break
  count=round(count+1/neq,1)
  print(count)
  ######################### ITERATOR CONVERGED ##############################
}

print(paste0("SNP EFFECTS ",marker_class," took ",count," iterations to solve",sep=" "))

sim_effects=allele_effect
plot(sol[-(nmeans)],sim_effects,main="SNP Effect Estimates vs True Values")

ebv<-X%*%sol
gv<-z%*%allele_effect
cor.test(gv,ebv)
plot(gv,ebv,main="Genomic Estimated Breeding Values (GEBVS) vs Genetic Values : Training Set")

######################### TRAINING SET EBV ESTIMATION -  FINISHED ##############################

##########################

### Genomic Prediction ###

##########################

if (genomic_prediction=="YES"){
  #read in validation set individuals (No Phen)
  z_pred=read.csv("Z_Validation.csv",sep=",",header=F);z_pred<-as.matrix(z_pred)
  class(z_pred)<-"numeric"
  ndata=nrow(z_pred)
  effects_estimate<-sol[-(nmeans)]
  
  gebv<-z_pred%*%effects_estimate
  gv = z_pred%*%sim_effects
  cor.test(gv,gebv)
  plot(gv,gebv,main="Genomic Estimated Breeding Values (GEBVS) vs Genetic Values : Validation Set")
  
  ######################### GENOMIC PREDICTION - FINISHED ##############################
}

######################################################################################

### Estimate Variance Components using Bayes C (π=0)

#ss=sum(sol_random^2)+length(e)*Sa
#vara=ss/chi(length(sol_random)+ncol(z))
#ss=sum(e**2)+nue*Se
#vare=ss/chi(nue+ndata) 

######################################################################################
