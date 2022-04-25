# Required function

# Matricization (Hoff (2015))
mat<-function(A,j)
{
  Aj<-t(apply(A,j,"c"))
  if(nrow(Aj)!=dim(A)[j])  { Aj<-t(Aj) }
  Aj
}

# Array-matrix product (Hoff (2015))
#' @param A a real valued array
#' @param M a real matrix
#' @param k an integer, a mode of \code{A}

amprod<-function(A,M,k)
{
  K<-length(dim(A))
  AM<-M%*%mat(A,k)
  AMA<-array(AM, dim=c(dim(M)[1],dim(A)[-k]) )
  aperm(AMA,  match(1:K,c(k,(1:K)[-k]) ) )
}

# For a tensor T \in R^{d_1, d_2, d_3, p} and a vector v \in R^p, compute mode-4 tensor vector product T \times_{4} v (Sun and Li (2017))
Tv = function(T, v){

  DD = dim(T)
  v = as.numeric(v)
  p = length(v)
  if(DD[length(DD)] != p){
    error("Dimensions of T and v do not fit")
  }else{

    tmp = array(0,DD[1:(length(DD) - 1)])
    for(i in 1:p){
      if(length(DD) == 3){
        tmp = tmp + T[,,i] * v[i]
      }else if(length(DD) == 4){
        tmp = tmp + T[,,,i] * v[i]
      }
    }
  }
  tmp
}

# For a tensor T \in R^{d_1, d_2, d_3, p} and a matrix M \in R^{p * n} , compute mode-4 tensor vector product T \times_{4} M_i for each i = 1, ...,n (Sun and Li (2017))
TM = function(T, M){
  n = ncol(M)
  out = list()
  for(i in 1:n){
    out[[i]] = Tv(T,M[,i])
  }
  out
}

# TM, then vectorize
TMv = function(T, M){
  n = ncol(M)
  out = list()
  for(i in 1:n){
    out[[i]] = as.vector(Tv(T,M[,i]))
  }
  out
}






