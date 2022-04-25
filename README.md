# A New Class of Skewed Tensor Distributions

This repository includes code for BSTT (Bayesian Skewed Tensor-t) model with sparsity inducing prior described in the paper "A New Class of Skewed Tensor Distributions", by Lee et al (2021+).

## Installation

The BSTT package can be installed with the `devtools` package:
  
  ```{r, eval = FALSE}
library(devtools) 
devtools::install_github(repo='InkooLee/BSTT')
```
## Documentation

Vignette is available at [/InkooLee/BSTT/blob/master/vignettes/BSTT_vignette.html](http://htmlpreview.github.io/?https://github.com/InkooLee/BSTT/blob/master/vignettes/BSTT_vignette.html)

## Usage

#### Bayesian Skewed Tensor-t with tensor spike-and-slab lasso/tensor normal prior

The `BSTT` package provides posterior sampling algorithm for skewed tensor-t model. 

  
  ```{r, eval = FALSE}
## Load library and required package
library(BSTT)
library(abind)

## load GAAD data
setwd('')
data <- read.csv("PDdata.csv")

## load required functions
source("./functions_tensor.R")
source("./BSTT_TSSL.R")
source("./BSTT_TN.R")

#'@ Inputs
#' Y \in R^{t x s x b x n} three way tensor of responses with n subjects
#' X \in R^{p} vector of subject-level covariates (with intercept, (p+1) x 1)

t = 28 # maximum number of teeth per one subject
s = 6 # maximum number of surfaces per each tooth
b = 2 # two bio-markers (i.e., PPD and CAL)
n = 290 # number of subjects

# Three way tensor response
YP <- array(data$PPD, dim = c(s,t,n)) # PPD
YC <- array(data$CAL, dim = c(s,t,n)) # CAL
YP <- aperm(YP, c(2,1,3)) # t x s matrices for PPD
YC <- aperm(YC, c(2,1,3)) # t x s matrices for CAL

Y <- abind(YP, YC) # t x s x nl
Y <- array(Y, dim = c(t,s,b,n)) # t x s x b x n

#'@PPD/CAL
#'Y[,,,1] # 290 subjects for PPD
#'Y[,,,2] # 290 subjects for CAL

## covariates
X <- data[,7:11]
X <- X[seq(1, nrow(X), 168), ]
rownames(X) <- NULL
X <- t(as.matrix(X)) # matrix of covariates (Age,Gender,BMI,Smoker,HbA1c) x number of subjects
X <- rbind(t(matrix(1,n)),X) # including intercept term
p <- dim(X)[1]

# Missing Responses
vecy <- t(mat(Y,4))
vecym <- matrix(NA, t*s*b, n)
vecym[which(!is.na(vecy))] <- 0; vecym[which(is.na(vecy))] <- 1; delta_p <- vecym; # delta_p is t*s*b by n matrix with missing = 1, observed = 0
```

We can fit the model as 

```{r, eval = FALSE}
## Fit the model
result_TSSL <- BSTT_TSSL(Y,X,vecy, n.burn = 100, n.save = 1000, thin = 5)
result_TN <- BSTT_TN(Y,X,vecy, n.burn = 100, n.save = 1000, thin = 5)
```


