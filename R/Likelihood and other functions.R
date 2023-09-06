###############################################################################
### Functions used in simulation ##############################################

#### Helper functions used in simulation study presented in:
#### Sollmann, R., "Estimating the temporal scale of lagged 
#### responses in species abundance and occurrence", Ecosphere (accepted).
#### For background/details, please see publication


## Function to generate correlated variable - 
## adapted from Chandler & Hepinstall-Cymerman, 2016,
## Supplementary Material 1; doi.org/10.1007/s10980-016-0380-z

spcov2 <- function(D, alpha=2, standardize=TRUE) {
  V <- exp(-D/alpha)
  cov1 <- t(chol(V)) %*% rnorm(nrow(D))
  if(standardize)
    z <- as.numeric((cov1 - mean(cov1)) / sd(cov1))
  else
    z <- as.numeric(cov1)
  return(z)
}


###############################################################################
##### LIKELIHOOD FUNCTIONS FOR TEMPORAL SCALE MODEL ###########################

##in all likelihood functions:

## Input:
## y. = observations (vectors for Poisson and logistic regression,
##      matrices for occupancy, N-mixture models)
## Xs. = static predictor, vector
## Xc. = repeated measures predictor, matrix

## Parameters to be estimated:
## sig=temporal scale parameter; here, of half-normal weights function
## beta0=intercept, on link scale, of response variable
## beta1=coefficient of static predictor variable 
## beta2=coefficient of composite predictor (ie, with repeatd measures over time)

## For occupancy, N-mixture models:
## alpha0 = detection probability (on logit scale)


###Poisson regression (perfectly observed count data)

nll <- function(params, y., Xs., Xc.) {
  sig <- exp(params[1]) 
  beta0 <- params[2]
  beta1 <- params[3]
  beta2 <- params[4]
  
  cov.w <- apply(Xc., 1, function(x) {
    w0 <- exp(-d^2 / (2*sig^2))
    w <- w0/sum(w0)
    sum(x * w)
  })
  
  lam <- exp(beta0+beta1*Xs.+beta2*cov.w)
  loglik <- dpois(y., lam, log=TRUE)
  -sum(loglik)
}


##Logistic regression (perfectly observed binary data)
nll.p <- function(params, y., Xs., Xc.) {
  sig <- exp(params[1])
  beta0 <- params[2]
  beta1 <- params[3]
  beta2 <- params[4]
  
  cov.w <- apply(Xc., 1, function(x) {
    w0 <- exp(-d^2 / (2*sig^2))
    w <- w0/sum(w0)
    sum(x * w)
  })
  
  psi <- plogis(beta0+beta1*Xs.+beta2*cov.w)
  loglik <- dbinom(y., 1, psi, log=TRUE)
  -sum(loglik)
}


##N-mixture model (imperfectly observed abundance)
## Function adapted from R package unmarked pcount() 
## (github.com/rbchan/unmarked/blob/master/R/pcount.R)
## Note that variable naming is different between this function and pcount()

nll.nm <- function(params, y., Xs., Xc.) {

  sig <- exp(params[1])
  beta0 <- params[2]
  beta1 <- params[3]
  beta2 <- params[4]
  alpha0 <- params[5] 
  
  cov.w <- apply(Xc., 1, function(x) {
    w0 <- exp(-d^2 / (2*sig^2))
    w <- w0/sum(w0)
    sum(x * w)
  })
  
  ##move prep into nll function
  KK<-max(y.)+100 #upper limit for Poisson integration
  kk <- 0:KK 
  J <- nrow(y.) #number of sites
  K <- ncol(y.) #number of occasions
  k.ik <- rep(kk, J) #abundance states, repeated for each site
  k.ijk <- rep(kk, J * K) #abundance states, repeated for each site and occasion
  y.ij <- as.numeric(t(y.))
  y.ijk <- rep(y.ij, each = KK + 1) #repeat data for each possible abundance state
  ijk <- expand.grid(kk = 0:KK, j = 1:K, i = 1:J)
  ijk.to.ikj <- with(ijk, order(i, kk, j))
  
  
  theta.i <- exp(beta0+beta1*Xs.+beta2*cov.w) #site-length vector of expected abundance
  p.ij <- plogis(rep(alpha0, J*K))        #site by occasion length vector of p
  
  theta.ik <- rep(theta.i, each = KK + 1) #exp. abundance, repeated for each abundance state
  p.ijk <- rep(p.ij, each = KK + 1)       #p, repeated for each abundance state
  bin.ijk <- dbinom(y.ijk, k.ijk, p.ijk)  #p(data) for all possible abundance states, some will be 0
  bin.ijk[which(is.na(bin.ijk))] <- 1     #does not apply here
  bin.ik.mat <- matrix(bin.ijk[ijk.to.ikj], J * (KK + 
                                                   1), K, byrow = TRUE)                #p(data) in matrix with ncol=number of occasions
  g.ik <- apply(bin.ik.mat, 1, prod)      #rowProds(bin.ik.mat)   #total p(data) across all visits, some will be 0 (when k<n.detected)
  #but for each site, across all abundance states, the sums should never be 0
  
  f.ik <- dpois(k.ik, theta.ik)  #p(abundance state) given exp. abundance, some could be 0 (extremes), but should not sum to 0
  #for a given site across all possible abundance states, ie, some should be possible
  
  dens.i.mat <- matrix(f.ik * g.ik, J, KK + 1, byrow = TRUE)
  dens.i <- apply(dens.i.mat, 1, sum)
  -sum(log(dens.i))
}



## Occupancy model (imperfectly observed abundance)
## Function adapted from R package unmarked occu() 
## (github.com/rbchan/unmarked/blob/master/R/occu.R)
nll.occ <- function(params, y., Xs., Xc.) {
  
  sig <- exp(params[1])
  beta0 <- params[2]
  beta1 <- params[3]
  beta2 <- params[4]
  alpha0 <- params[5] 
  
  nd <- ifelse(rowSums(y., na.rm = TRUE) == 0, 1, 0)
  yvec <- as.numeric(t(y.)) 
  J<-length(nd)
  K<-length(yvec)/J
  
  cov.w <- apply(Xc., 1, function(x) {
    w0 <- exp(-d^2 / (2*sig^2))
    w <- w0/sum(w0)
    sum(x * w)
  })
  
  psi <- plogis(beta0+beta1*Xs.+beta2*cov.w)
  pvec <- rep(plogis(alpha0), J*K)
  cp <- (pvec^yvec) * ((1 - pvec)^(1 - yvec))
  cpmat <- matrix(cp, J, K, byrow = TRUE)
  loglik <- log(apply(cpmat,1,prod) * psi + nd * (1 - psi))
  -sum(loglik)
}


