#Gauss-Seidel Residual Updater
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

