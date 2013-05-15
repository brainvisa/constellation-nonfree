#faire un p4 add
library(cluster)
coucoukmedoids <- function(x, k, display){
  cl <- pam(x, k)
  summary(cl)
  if (display) {
  plot(cl)
  }
  return(cl)
}

kmedoidsondist <- function(x, k, init_medoids){
  cl <- pam(x, k, diss=TRUE, medoids=init_medoids)
  summary(cl)
  
  return(cl)
}

kmedoidsondist2 <- function(x, k){
  cl <- pam(x, k, diss=TRUE)
  summary(cl)
  
  return(cl)
}

kmedoidsWithDistMatrix <- function(x, k){
  cl <- pam(x, k, diss=TRUE)
  summary(cl)
  return(cl)
}