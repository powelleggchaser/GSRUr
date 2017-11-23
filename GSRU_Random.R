#Gaussian Seidel Residual Updater - Multiple Markers
#Rversion

###########################################################

##  FUNCTIONS

###########################################################

t_mult<-function(matrix){
  t_matrix<-t(matrix)
  prod_matrix<-t_matrix%*%matrix
  return(prod_matrix)
}

t_mult_dupl<-function(matrix1,matrix2){
  t_matrix<-t(matrix1)
  prod_matrix<-t_matrix%*%matrix2
  return(prod_matrix)
}

t_mult_diag<-function(matrix){
  t_matrix<-t(matrix)
  prod_matrix<-t_matrix%*%matrix
  diag_prod_matrix<-diag(prod_matrix)
  return(diag_prod_matrix)
}

dupl_mult_diag<-function(matrix1,matrix2){
  prod_matrix<-matrix2%*%matrix1
  diag_matrix<-diag(prod_matrix)
  return(diag_matrix)
}

lamda<-function(vare,vara){
  vare<-as.numeric(vare)
  vara<-as.numeric(vara)
  lamda=vare/vara
  return(lamda)
}

normalise <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

ssqd <- function(x,y) {
  ssqd <-sum((x-y)^2)
  return(ssqd)
}

calc_allele_freq_per_locus <- function(geno,allele_effect,pheno) {
  p=(colSums(geno))/(2*ndata)
  q=1-p
  pq2=2*p*q
  var_a=pq2*(allele_effect^2)
  var_p=var(phenotype)
  var_e=var_p-var_a
  lamda=lamda(var_e,var_a)
  output = cbind(p,q,pq2,var_a,var_p,var_e,lamda)
  #output = cbind(p,q,pq2,var_a,var_p,var_e)
  return(output)
}

###########################################################
###########################################################


### DATA ENTRY ###


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

#calculate allele frequencies & variances
allele_freq_stats=calc_allele_freq_per_locus(z,allele_effect,phenotype)

lamda1=0.1

### MATRIX MANIPULATION ###

#assign phenotypes to y matrix
y=phenotype

#xpy<-as.matrix(t_mult_dupl(x,y))
#zpy<-as.matrix(t_mult_dupl(z,y))

#combine fixed effects and random effects into a single matrix
X=cbind(x,z)

#X<-z

#Count the number of expected solutions
neq=ncol(X)
#Count the number of expected of fixed effect means
nmeans=ncol(x)

#provide phenotypes as the starting values for residual updater
e=as.matrix(y)
#set starting convergence value
CONV=1
#Set convergence threshold
CONV_TOL=0.00000000000001

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


###############################################################

### FIT SNPS AS RANDOM EFFECTS ################################

###############################################################

#create intital empty vector to store newly generated solutions at each iteration
sol_random=rep(0,neq)
#create intial empty vector to store previously generated solutions at each iteration
old_sol_random=rep(0,neq)

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
    rhs=(X[,j]%*%e) + (xpx[j]*sol_random[j])
    #do Gauss Seidel
    val=rhs/lhs
    #MCMC sample solution from its conditional (commented out here)
    #val_num=(rhs/lhs)
    #val_den=(1/lhs)
    #val=rnorm(1,(rhs/lhs),(1/lhs))
    #update e with current estimate (formula 2)
    e = e - (X[,j]*(val-sol_random[j]))
    #store old solution
    old_sol_random[j]=sol_random[j]
    #update sol
    sol_random[j]=val
    #Calcualte convergence value using the SumSquare of old vs new solutions
    ssqd_on[j]<-ssqd(old_sol_random[j],sol_random[j])
    ssq_n[j]<-sum(sol_random[j]^2)
    CONV<-ssqd_on[j]/ssq_n[j]
  }
  if (CONV<=CONV_TOL) break
  count=round(count+1/neq,1)
  print(count)
}

print(paste0("SNP EFFECTS (RANDOM) took ",count," iterations to solve",sep=" "))

sim_effects=allele_effect
plot(sol_random[-(nmeans)],sim_effects,main="SNP Effect Estimates vs True Values")

ebv<-X%*%sol_random
gv<-z%*%allele_effect
cor.test(gv,ebv)
plot(gv,ebv,main="Genomic Estimated Breeding Values (GEBVS) vs Genetic Values : Training Set")
######################################################################################

##########################

### Genomic Prediction ###

##########################


#read in validation set individuals (No Phen)
z_pred=read.csv("Z_Validation.csv",sep=",",header=F);z_pred<-as.matrix(z_pred)
class(z_pred)<-"numeric"
ndata=nrow(z_pred)
effects_estimate<-sol_random[-(nmeans)]

gebv<-z_pred%*%effects_estimate
gv = z_pred%*%sim_effects
cor.test(gv,gebv)
plot(gv,gebv,main="Genomic Estimated Breeding Values (GEBVS) vs Genetic Values : Validation Set")

###########################

### Estimate Variance Components using Bayes C (π=0)

#ss=sum(sol_random^2)+length(e)*Sa
#vara=ss/chi(length(sol_random)+ncol(z))
#ss=sum(e**2)+nue*Se
#vare=ss/chi(nue+ndata) 
