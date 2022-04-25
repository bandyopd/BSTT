#' Fit Bayesian Skewed Tensor-t (BSTT) model with tensor spike-and-slab lasso (TSSL) prior to GAAD data
#' Lee et al. (2021+) A New Class of Skewed Tensor Distributions
#'
#' @param Y t x s x b x n array of skewed tensor response
#' @param X p x n matrix of covariates
#' @param vecy tsb x n matrix: vectorized tensor
#' @param n.burn burn-in period
#' @param n.save number of posterior samples
#' @param thin thinning size
#'
#' @return Returns a list with the following components:
#' \item rho: posterior samples of correlation of each mode of tensor response (3 x n.save)
#' \item sigma.sq: posterior samples of variance parameter (b x n.save)
#' \item lam.est: posterior samples of skewness parameters (b x n.save)
#' \item B.est: posterior samples of tensor regression coefficients (t x s x b x p x n.save)
#' \item nu: posterior samples of degree of freedom (1 x n.save)
#' @export
#'

BSTT_TN <- function(Y,X,vecy, n.burn = 10, n.save = 100, thin = 1){

  t = dim(Y)[1]; s = dim(Y)[2]; b = dim(Y)[3]; p = dim(X)[1]; n = dim(X)[2];

  # Store MCMC out

  rho.save <- matrix(NA, 2, n.save)
  B.est.save <- array(NA, c(p, t*s, n.save))
  sigma.sq.save <- matrix(NA, 1, n.save)
  lam.est.save <- matrix(NA, 1, n.save)
  nu.save <- matrix(NA, 1, n.save)

  source("./functions_tensor.R")

  #load required packages
  library(expm); library(Matrix); library(matrixcalc); library(LaplacesDemon)
  library(MASS); library(msm); library(truncnorm); library(abind)
  library(magrittr); library(doParallel); library(TruncatedNormal)
  registerDoParallel(cores=2)


  #  Updating functions
  #---------------------------------------------------------------------------------

  # Update rho1 & rho2 & rho3 (MH algorithm) using closed form of equicorrelation assumption for R1 & R2 & R3

  rho.update <- function(t,s,b,rho,vecy,X,B.est,lam.est,W,sigma.sq, R1, R2,delta){

    kron<-function(...)
    {
      M<-list(...)
      if(is.list(M[[1]])) { M<-M[[1]] }
      JM<-1 ; for(j in 1:length(M)){ JM<-kronecker(JM,M[[j]]) }
      JM
    }

    rho.prop = c(rbeta(1,2.5,2), rbeta(1,2.5,2), rbeta(1,2.5,2)) #(rho1.prop, rho2.prop, rho3.prop)

    det.curr = ((sigma.sq[1]*sigma.sq[2]*((1-rho[3])^(b-1))*(1 + (b-1)*rho[3]))^(t*s))*((((1-rho[2])^(s-1))*(1 + (s-1)*rho[2]))^(t*b))*((((1-rho[1])^(t-1))*(1 + (t-1)*rho[1]))^(s*b))
    det.prop = ((sigma.sq[1]*sigma.sq[2]*((1-rho.prop[3])^(b-1))*(1 + (b-1)*rho.prop[3]))^(t*s))*((((1-rho.prop[2])^(s-1))*(1 + (s-1)*rho.prop[2]))^(t*b))*((((1-rho.prop[1])^(t-1))*(1 + (t-1)*rho.prop[1]))^(s*b))

    inv.curr = kron(diag(c(1/sqrt(sigma.sq[1]),1/sqrt(sigma.sq[2])))%*%((1/(1-rho[3]))*(diag(b) - (rho[3]/ (1 + (b-1)*rho[3]))*matrix(1,b,b)))%*%diag(c(1/sqrt(sigma.sq[1]),1/sqrt(sigma.sq[2]))), ((1/(1-rho[2]))*(diag(s) - (rho[2]/ (1 + (s-1)*rho[2]))*matrix(1,s,s))), ((1/(1-rho[1]))*(diag(t) - (rho[1]/ (1 + (t-1)*rho[1]))*matrix(1,t,t)))) # \{D_{\sigma}^{-1}R_{\rho_3}^{-1}D_{\sigma}^{-1} \otimes R_{\rho_2}^{-1} \otimes R_{\rho_1}^{-1}\}
    inv.prop = kron(diag(c(1/sqrt(sigma.sq[1]),1/sqrt(sigma.sq[2])))%*%((1/(1-rho.prop[3]))*(diag(b) - (rho.prop[3]/ (1 + (b-1)*rho.prop[3]))*matrix(1,b,b)))%*%diag(c(1/sqrt(sigma.sq[1]),1/sqrt(sigma.sq[2]))), ((1/(1-rho.prop[2]))*(diag(s) - (rho.prop[2]/ (1 + (s-1)*rho.prop[2]))*matrix(1,s,s))), ((1/(1-rho.prop[1]))*(diag(t) - (rho.prop[1]/ (1 + (t-1)*rho.prop[1]))*matrix(1,t,t))))

    inv.R21.curr = kron( ((1/(1-rho[2]))*(diag(s) - (rho[2]/ (1 + (s-1)*rho[2]))*matrix(1,s,s))), ((1/(1-rho[1]))*(diag(t) - (rho[1]/ (1 + (t-1)*rho[1]))*matrix(1,t,t))) )
    inv.R21.prop = kron( ((1/(1-rho.prop[2]))*(diag(s) - (rho.prop[2]/ (1 + (s-1)*rho.prop[2]))*matrix(1,s,s))), ((1/(1-rho.prop[1]))*(diag(t) - (rho.prop[1]/ (1 + (t-1)*rho.prop[1]))*matrix(1,t,t))) )

    inv.R3.curr = ((1/(1-rho[3]))*(diag(b) - (rho[3]/ (1 + (b-1)*rho[3]))*matrix(1,b,b)))
    inv.R3.prop = ((1/(1-rho.prop[3]))*(diag(b) - (rho.prop[3]/ (1 + (b-1)*rho.prop[3]))*matrix(1,b,b)))

    logdens.curr = -0.5*n*log(det.curr) - 0.5*sum( (t(vecy) - t(X)%*%B.est - t(W)%*%(kron(diag(lam.est,2,2),diag(t*s))) )%*%inv.curr%*%(vecy - t(B.est)%*%X - (kron(diag(lam.est,2,2),diag(t*s))%*%W) ) )
    logdens.prop = -0.5*n*log(det.prop) - 0.5*sum( (t(vecy) - t(X)%*%B.est - t(W)%*%(kron(diag(lam.est,2,2),diag(t*s))) )%*%inv.prop%*%(vecy - t(B.est)%*%X - (kron(diag(lam.est,2,2),diag(t*s))%*%W) ) )

    logdens.curr = - 0.5*sum( (t(vecy) - t(X)%*%B.est - t(W)%*%(kron(diag(lam.est,2,2),diag(t*s))) )%*%inv.curr%*%(vecy - t(B.est)%*%X - (kron(diag(lam.est,2,2),diag(t*s))%*%W) ) )
    logdens.prop = - 0.5*sum( (t(vecy) - t(X)%*%B.est - t(W)%*%(kron(diag(lam.est,2,2),diag(t*s))) )%*%inv.prop%*%(vecy - t(B.est)%*%X - (kron(diag(lam.est,2,2),diag(t*s))%*%W) ) )

    logratio = logdens.prop - logdens.curr
    if(log(runif(1)) > logratio) {rho = rho} else {rho = rho.prop}
    if(log(runif(1)) > logratio) {inv.Sigma = inv.curr} else {inv.Sigma = inv.prop}
    if(log(runif(1)) > logratio) {inv.R21 = inv.R21.curr} else {inv.R21 = inv.R21.prop}
    if(log(runif(1)) > logratio) {inv.R3 = inv.R3.curr} else {inv.R3 = inv.R3.prop}
    if(log(runif(1)) > logratio) {det = det.curr} else {det = det.prop}

    return(list(rho, inv.Sigma, inv.R21, inv.R3, det))


  }

  #---------------------------------------------------------------------------------
  # Update sigma.sq

  sigma.sq.update <- function (t,s,b,n,vecy,Y,X,B.est,lam.est,W,inv.R,sigma.sq,rho,delta,inv.R21,inv.R3){

    kron<-function(...)
    {
      M<-list(...)
      if(is.list(M[[1]])) { M<-M[[1]] }
      JM<-1 ; for(j in 1:length(M)){ JM<-kronecker(JM,M[[j]]) }
      JM
    }

    g1 = 2; g2 = 2; #prior distribution for inv.sigma.sq ~ Ga(2,2)
    Sww <-  foreach(l = 1:n,.combine='+') %dopar% {crossprod(W[1:(t*s),l])}
    S <- foreach(l = 1:n,.combine='+') %dopar% { ( c(Y[,,1,l]) - t(X)[l,]%*%B.est[,1:(t*s)] - lam.est[1]*t(W[1:(t*s),l]) )%*%kron(inv.R3,inv.R21)[1:(t*s),1:(t*s)]%*%t( c(Y[,,1,l]) - t(X)[l,]%*%B.est[,1:(t*s)] - lam.est[1]*t(W[1:(t*s),l]) ) }
    inv.sigma.sq <- rgamma(1, g1 + n*t*s, 0.5*S + 0.5*Sww + g2)
    sigma1.sq <- 1/inv.sigma.sq

    Sww2 <- foreach(l = 1:n,.combine='+') %dopar% {crossprod(W[(t*s + 1):(2*t*s),l])}
    S2 <- as.numeric(foreach(l = 1:n,.combine='+') %dopar% { ( c(Y[,,2,l]) - t(X)[l,]%*%B.est[,(t*s + 1):(2*t*s)] - lam.est[2]*t(W[(t*s + 1):(2*t*s),l]) )%*%kron(inv.R3,inv.R21)[(t*s + 1):(2*t*s),(t*s + 1):(2*t*s)]%*%t( c(Y[,,2,l]) - t(X)[l,]%*%B.est[,(t*s + 1):(2*t*s)] - lam.est[2]*t(W[(t*s + 1):(2*t*s),l]) ) })
    inv.sigma.sq2 <- rgamma(1, g1 + n*t*s, 0.5*S2 + 0.5*Sww2 + g2)
    sigma2.sq <- 1/inv.sigma.sq2

    return(c(sigma1.sq,sigma2.sq))
  }

  #---------------------------------------------------------------------------------
  # Update skewness parameter lam.est

  lam.est.update <- function(W,inv.Sigma,B.est,Y,X,n,vecy,lam.est,t,s,b,delta){

    kron<-function(...)
    {
      M<-list(...)
      if(is.list(M[[1]])) { M<-M[[1]] }
      JM<-1 ; for(j in 1:length(M)){ JM<-kronecker(JM,M[[j]]) }
      JM
    }

    Swinvw <- foreach(l = 1:n,.combine='+') %dopar% {t(W[1:(t*s),l])%*%inv.Sigma[1:(t*s),1:(t*s)]%*%W[1:(t*s),l]}
    A.lam <- (b^2)*Swinvw + 1
    B.lam <- foreach(l = 1:n,.combine='+') %dopar% {( c(Y[,,1,l]) - t(X)[l,]%*%B.est[,1:(t*s)])%*%inv.Sigma[1:(t*s),1:(t*s)]%*%W[1:(t*s),l]*(b^2) + 4 }
    lam1.est <- rnorm(1, mean = B.lam/A.lam, sd = 1/(sqrt(A.lam)))

    Swinvw2 <- foreach(l = 1:n,.combine='+') %dopar% {t(W[(t*s + 1):(2*t*s),l])%*%inv.Sigma[(t*s + 1):(2*t*s),(t*s + 1):(2*t*s)]%*%W[(t*s + 1):(2*t*s),l]}
    A.lam2 <- (b^2)*Swinvw2 + 1
    B.lam2 <- foreach(l = 1:n,.combine='+') %dopar% {( c(Y[,,2,l]) - t(X)[l,]%*%B.est[,(t*s + 1):(2*t*s)])%*%inv.Sigma[(t*s + 1):(2*t*s),(t*s + 1):(2*t*s)]%*%W[(t*s + 1):(2*t*s),l]*(b^2) + 4}
    lam2.est <- rnorm(1, mean = B.lam2/A.lam2, sd = 1/(sqrt(A.lam2)))

    return(c(lam1.est,lam2.est))

  }

  #---------------------------------------------------------------------------------
  # Update W = abs(Z_2) : Note that the each component is sampled from univariate truncated normal distribution

  W.update <- function(t,s,n, lam.est, inv.Sigma, B.est, vecy, X, W, nu){

    kron<-function(...)
    {
      M<-list(...)
      if(is.list(M[[1]])) { M<-M[[1]] }
      JM<-1 ; for(j in 1:length(M)){ JM<-kronecker(JM,M[[j]]) }
      JM
    }

    nu_W <- t*s + nu

    for (N in 1:n){
      W[1:(t*s),N] <- rtmvt(n = 1, mu = (solve( (lam.est[1]^2)*inv.Sigma[1:(t*s),1:(t*s)]/delta[N] + diag(t*s))%*%((lam.est[1]*inv.Sigma[1:(t*s),1:(t*s)]/delta[N])%*%(vecy[1:(t*s),] - t(B.est)[1:(t*s),]%*%X)) )[,N] , sigma = c(( solve((lam.est[1]^2)*inv.Sigma[1:(t*s),1:(t*s)]/delta[N] + diag(t*s)) + nu )/nu_W)*(diag(t*s) + (lam.est[1]^2)*inv.Sigma[1:(t*s),1:(t*s)]*delta[N]), df = nu_W, lb = rep(0,t*s), ub = rep(Inf,t*s))
    }

    for (N in 1:n){
      W[(t*s + 1):(2*t*s),N] <- rtmvt(n = 1, mu = (solve( (lam.est[2]^2)*inv.Sigma[(t*s+1):(2*t*s),(t*s+1):(2*t*s)]/delta[N] + diag(t*s))%*%((lam.est[2]*inv.Sigma[(t*s+1):(2*t*s),(t*s+1):(2*t*s)]/delta[N])%*%(vecy[(t*s + 1):(2*t*s),] - t(B.est)[(t*s + 1):(2*t*s),]%*%X)) )[,N] , sigma = c(( solve((lam.est[2]^2)*inv.Sigma[(t*s + 1),(2*t*s)]/delta[N] + diag(t*s)) + nu )/nu_W)*(diag(t*s) + (lam.est[2]^2)*inv.Sigma[(t*s+1):(2*t*s),(t*s+1):(2*t*s)]*delta[N]), df = nu_W, lb = rep(0,t*s), ub = rep(Inf,t*s))
    }

    return(W)

  }


  #---------------------------------------------------------------------------------
  # Update beta
  B.est.update <- function(X,vecy,W,lam.est, inv.Sigma, inv.R1, inv.R2, t,s,p, B.est){

    kron<-function(...)
    {
      M<-list(...)
      if(is.list(M[[1]])) { M<-M[[1]] }
      JM<-1 ; for(j in 1:length(M)){ JM<-kronecker(JM,M[[j]]) }
      JM
    }

    Sxy1 <- foreach(l = 1:n,.combine='+') %dopar% { t(t(X[,l]))%*%t(as.matrix(t(vecy)[l,1:(t*s)]))%*%inv.Sigma[1:(t*s),1:(t*s)] }
    Sxw1 <- foreach(l = 1:n,.combine='+') %dopar% { t(t(X[,l]))%*%t(as.matrix(t(W)[l,1:(t*s)]))%*%inv.Sigma[1:(t*s),1:(t*s)]*lam.est[1] }
    Sxy2 <- foreach(l = 1:n,.combine='+') %dopar% { t(t(X[,l]))%*%t(as.matrix(t(vecy)[l,(t*s+1):(2*t*s)]))%*%inv.Sigma[(t*s+1):(2*t*s),(t*s+1):(2*t*s)] }
    Sxw2 <- foreach(l = 1:n,.combine='+') %dopar% { t(t(X[,l]))%*%t(as.matrix(t(W)[l,(t*s+1):(2*t*s)]))%*%inv.Sigma[(t*s+1):(2*t*s),(t*s+1):(2*t*s)]*lam.est[2] }
    Sxx <- foreach(l = 1:n,.combine='+') %dopar% { t(t(X[,l]))%*%t(t(X)[l,]) }
    A.inv1 <- solve(kron(Sxx,inv.Sigma[1:(t*s),1:(t*s)]) + 0.01*diag(t*s*p))
    A.inv2 <- solve(kron(Sxx,inv.Sigma[(t*s+1):(2*t*s),(t*s+1):(2*t*s)]) + 0.01*diag(t*s*p))
    vecB1 <- mvrnorm(1, mu = A.inv1%*%(as.vector(t(Sxy1)) - as.vector(t(Sxw1))), Sigma = A.inv1, tol = 1e-3)
    vecB2 <- mvrnorm(1, mu = A.inv2%*%(as.vector(t(Sxy2)) - as.vector(t(Sxw2))), Sigma = A.inv1, tol = 1e-3)
    B.est <- cbind(mat(array(vecB1, dim = c(t,s,p)),3), mat(array(vecB2, dim = c(t,s,p)),3))

  }

  # Update nu.est # Q_fun (nu,delta)
  nu.est.update <- function(n,nu,delta,sigma.sq){

    source("functions_tensor.R")
    log.nu.curr = (0.5*n*nu)*log(0.5*nu) - n*log(gamma(0.5*nu)) + (0.5*nu - 1)*sum(log(delta)) - (0.5*nu)*sum(delta) - log(sigma.sq) + 0.5*log( (gamma(0.5*nu) * digamma(0.5*nu) * digamma(0.5*nu) + gamma(0.5*nu) * trigamma(0.5*nu))/gamma(0.5*nu) - gamma(0.5*nu) * digamma(0.5*nu) * (gamma(0.5*nu) * digamma(0.5*nu))/gamma(0.5*nu)^2 - ( (gamma(0.5*(nu+1)) * digamma(0.5*(nu+1)) * digamma(0.5*(nu+1)) + gamma(0.5*(nu+1)) * trigamma(0.5*(nu+1)))/gamma(0.5*(nu+1)) - gamma(0.5*(nu+1)) * digamma(0.5*(nu+1)) * (gamma(0.5*(nu+1)) * digamma(0.5*(nu+1)))/gamma(0.5*(nu+1))^2 ) - (2*(nu+3))/(nu*(nu+1)^2) )

    q_1st <- q_1st_part(nu) + 0.5*log(sum(delta)) - 0.5*sum(delta)
    q_2nd <- q_2nd(nu)
    mu_nu.prop =  nu - q_1st/q_2nd
    var_nu.prop = -1/q_2nd
    nu.prop = rtruncnorm(n=1, a = 2, b = 300, mean = mu_nu.prop, sd = sqrt(var_nu.prop))
    log.nu.prop = (0.5*n*nu.prop)*log(0.5*nu.prop) - n*log(gamma(0.5*nu.prop)) + (0.5*nu.prop - 1)*sum(log(delta)) - (0.5*nu)*sum(delta) - log(sigma.sq) + 0.5*log( (gamma(0.5*nu.prop) * digamma(0.5*nu.prop) * digamma(0.5*nu.prop) + gamma(0.5*nu.prop) * trigamma(0.5*nu.prop))/gamma(0.5*nu.prop) - gamma(0.5*nu.prop) * digamma(0.5*nu.prop) * (gamma(0.5*nu.prop) * digamma(0.5*nu.prop))/gamma(0.5*nu.prop)^2 - ( (gamma(0.5*(nu.prop+1)) * digamma(0.5*(nu.prop+1)) * digamma(0.5*(nu.prop+1)) + gamma(0.5*(nu.prop+1)) * trigamma(0.5*(nu.prop+1)))/gamma(0.5*(nu.prop+1)) - gamma(0.5*(nu.prop+1)) * digamma(0.5*(nu.prop+1)) * (gamma(0.5*(nu.prop+1)) * digamma(0.5*(nu.prop+1)))/gamma(0.5*(nu.prop+1))^2 ) - (2*(nu.prop+3))/(nu.prop*(nu.prop+1)^2) )

    logratio = log.nu.curr - log.nu.prop
    if(log(runif(1)) > logratio[1]) {nu = nu} else {nu = nu.prop}

    return(nu)

  }

  #initial values

  B.est <- matrix(rep(c(1,0.01,0.2,0.2,0.2,0.2),t*s*b), p , t*s*b)
  W <- matrix(0.2, t*s*b,n)
  rho <- matrix(c(0.6,0.6,0.6), 3, 1)
  sigma.sq <- c(1,1)
  lam.est <- c(1,1)
  nu = 4
  delta <- c(rep(1,n));b_delta <- c(rep(1,n));

  R1 <- matrix(0,t,t)
  for(j in 1:t){
    for(k in 1:t){
      if (j != k ){
        R1[j,k] = 0.6
      }
      else{
        R1[j,k] <- 1
      }
    }
  }

  R2 <- matrix(0,s,s)
  for(j in 1:s){
    for(k in 1:s){
      if (j != k ){
        R2[j,k] = 0.6
      }
      else{
        R2[j,k] <- 1
      }
    }
  }

  R3 <- matrix(0,b,b)
  for(j in 1:b){
    for(k in 1:b){
      if (j != k ){
        R3[j,k] = 0.6
      }
      else{
        R3[j,k] <- 1
      }
    }
  }

  kron<-function(...)
  {
    M<-list(...)
    if(is.list(M[[1]])) { M<-M[[1]] }
    JM<-1 ; for(j in 1:length(M)){ JM<-kronecker(JM,M[[j]]) }
    JM
  }

  # fill missing responses with the new missing values
  for (N in 1:n){
    vecy[,N] = ifelse (is.na(vecy[,N]), MASS::mvrnorm(length(delta_p[,N][(delta_p[,N]==1)]), mu = t(B.est)%*%X[,N] + kron(diag(lam.est),diag(t*s))%*%W[,N], Sigma = kron(R3,R2,R1) ) , vecy[,N])
  }
  Y <- array(vecy, dim = c(t,s,b,n))

  begin_sampler <- proc.time()[3]

  # MCMC iterations
  for (i in 1:(n.burn + n.save*thin)) { #}

    mat<-function(A,j)
    {
      Aj<-t(apply(A,j,"c"))
      if(nrow(Aj)!=dim(A)[j])  { Aj<-t(Aj) }
      Aj
    }

    # Kronecker product
    kron<-function(...)
    {
      M<-list(...)
      if(is.list(M[[1]])) { M<-M[[1]] }
      JM<-1 ; for(j in 1:length(M)){ JM<-kronecker(JM,M[[j]]) }
      JM
    }

    # Integrating missing values out
    for (N in 1:n){
      vecy[,N] = ifelse (is.na(vecy[,N]), mvrnorm(length(kappa[,N][(kappa[,N]==1)]), mu = t(B.est)%*%X[,N] + diag( c(rep(lam.est[1,1], t*s), rep(lam.est[2,2],t*s)),t*s*b)%*%W[,N], Sigma = kron(R3,R2,R1)) , vecy[,N])
    }

    rho.result = rho.update(t,s,b,rho,vecy,X,B.est,lam.est,W,sigma.sq, R1, R2,delta)
    rho = rho.result[[1]]; inv.Sigma = rho.result[[2]]; inv.R21 = rho.result[[3]]; inv.R3 = rho.result[[4]]; det = rho.result[[5]];

    R1 <- matrix(0,t,t)
    for(j in 1:t){
      for(k in 1:t){
        if (j != k ){
          R1[j,k] = rho[1]
        }
        else{
          R1[j,k] <- 1
        }
      }
    }

    R2 <- matrix(0,s,s)
    for(j in 1:s){
      for(k in 1:s){
        if (j != k ){
          R2[j,k] = rho[2]
        }
        else{
          R2[j,k] <- 1
        }
      }
    }

    R3 <- matrix(0,b,b)
    for(j in 1:b){
      for(k in 1:b){
        if (j != k ){
          R3[j,k] = rho[3]
        }
        else{
          R3[j,k] <- 1
        }
      }
    }

    sigma.sq = sigma.sq.update(t,s,b,n,vecy,Y,X,B.est,lam.est,W,inv.R,sigma.sq,rho,delta,inv.R21,inv.R3)
    lam.est = lam.est.update(W,inv.Sigma,B.est,Y,X,n,vecy,lam.est,t,s,b,delta)
    W = W.update(t,s,n, lam.est, inv.Sigma, B.est, vecy, X, W, nu)
    B.est =  B.est.update(X,vecy,W,lam.est, inv.Sigma, inv.R1, inv.R2, t,s,p, B.est)
    a_delta <- nu/2 + t*s*b;
    for (l in 1:290){
      b_delta[l] <- nu/2 + 0.5* t(vecy - t(B.est)%*%X - kron(diag(lam.est),diag(t*s))%*%W)[l,]%*%inv.Sigma%*%( vecy - t(B.est)%*%X -kron(diag(lam.est),diag(t*s))%*%W )[,l]
      delta[l] <- rgamma(1, a_delta, b_delta[l])
    }

    nu = nu.est.update(n,nu,delta,sigma.sq)

    ##Keep track of MCMC output:
    if(i > n.burn & (i - n.burn)%%thin==0){
      ii = (i - n.burn)/thin

      rho.save[,ii] <- rho
      sigma.sq.save[,ii] <- sigma.sq
      lam.est.save[,ii] <- lam.est
      B.est.save[,,ii] <- B.est
      nu.save[,ii] <- nu

    }

    if(i %% ceiling((n.save)/10) == 0){
      cat(
        paste0(
          "##### ",
          Sys.time(),
          " Iteration # ",i," of ", (n.save*thin + n.burn),
          " #####\n",
          "##### Time elapsed: ",proc.time()[3] - begin_sampler, " seconds #####\n",
          "##### Time per each iteration: ",(proc.time()[3] - begin_sampler)/i, " seconds #####\n"
        )
      )
    }

  } # end sampler


  result = list(rho.save, sigma.sq.save, lam.est.save, B.est.save,nu.save)
  names(result) = c('rho', 'sigma.sq','lam.est','B.est','nu')

  return(result)

} # end BSTT_TN function


