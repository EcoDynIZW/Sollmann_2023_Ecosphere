################################################################################
#### Simulating and analyzing data under temporal scale model ##################

#### Script to repeat the simulation study presented in: 
#### Sollmann, R., "Estimating the temporal scale of lagged 
#### responses in species abundance and occurrence", Ecosphere (in press).
#### For background/details on parameters, please see publication

##load likelihood, helper functions
source('Likelihood and other functions.R')

### Setting scenarios based on sigma, J, beta.c

##temporal scale parameter
sigma<-c(1, 2.5, 4) #smaller sigma - closer years matter more

##sample size
J<- c(100, 175, 250) 

##effect of composite variable
beta.c<- c(0.5, 0.9, 1.3)

## Create data frame with all scenario combinations of parameters
tr<-expand.grid(sigma, J, beta.c)
sig.in<-tr$Var1
J.in<-tr$Var2
b.in<-tr$Var3

n.scen<-dim(tr)[1]

#intercept abundance, occurrence
beta0 <- 0.7 

## effect of additional static covariate
beta1 <- 0.5 

##correlation of repeated measures of predictor across time
rho<-1 #higher rho --> higher correlation 

tt<-10 #years before focal year
d<-tt:1 #time since focal year, assuming predictor is arranged chronologically
#ie, first column is earliest year, thus farthest away from focal year
# can also be started at 0 instead of 1

##matrix of distances between years (temp distances)
D<-matrix(NA, tt,tt)
for (i in 1:tt){
  if (i ==1) {D[i,]<-(i-1):(tt-1)} else{
    D[i,]<-c((i-1):1, 0:(tt-i))}
}

##used in occu and Nmix only
alpha0<-0 #intercept p, logit scale
K<-3 #secondary occasions

#number of simulations per scenario
nsim=100

# manually set model type for which to run simulation, out of
# 'Poisson', 'LogReg', 'Occu', 'Nmix'
# Note: Nmix is fairly slow 

mod<-'Occu'

##set number of parameters estimated depending on model type
npar<-ifelse(mod %in% c('Poisson', 'LogReg'),4,5)

#objects to hold param estimates and convergence info
out<-array(NA, c(nsim, npar, n.scen))
se.out<-array(NA, c(nsim, npar, n.scen))
conv<-matrix(NA, nsim, n.scen)


## run simulations for all scenarios
for (k in 1:n.scen){ ##scenario loop
  
  #generate static covariate, always the same for a given scenario
  Xs<-rnorm(J.in[k],0,1) 
  
  ##matrix/array to hold all generated data
  if(mod %in% c('Occu', 'Nmix')){
  ally<-array(NA, c(nsim, J.in[k], K))} else{
    ally<-matrix(NA, nsim, J.in[k])
  }
  
  ##array to hold all generated composite variables
  allX<-array(NA, c(nsim, J.in[k], tt+1))
  
  for (i in 1:nsim){  #sim loop
    
    ##simulate temporally correlated covariate
    Xc<-matrix(NA, J.in[k], tt)
    for(j in 1:J.in[k]){
      cov <- spcov2(D, alpha=rho, standardize = F)#
      Xc[j,]<-cov 
    }
    
    ##generate weights based on input sigma
    w0 <- exp(-d^2 / (2*sig.in[k]^2))
    w<-w0/sum(w0)
    
    ##calculate composite variable
    Xcomp<-apply(Xc, 1, function(x){sum(x*w)})
    
    ##generate data, fit model, separate for each model type
    if (mod=='Poisson'){
      ##calculate exp. abundance, generate abundances
      lam<-exp(beta0+beta1*Xs+b.in[k]*Xcomp)
      N<-rpois(J.in[k],lam)
      y<-N
      
      ##fit model
      o1 <- try(optim(c(log(sig.in[k]), 0, 0, 0), nll, y.=y, Xs.=Xs, Xc.=Xc, hessian=TRUE,
                      method='Nelder-Mead', #for Poisson only
                      control = list(maxit=1000)))
    }
    
    if (mod=='LogReg'){
      
      ##calculate occu prob, generate occurrences
      psi<-plogis(beta0+beta1*Xs+b.in[k]*Xcomp)
      Z<-rbinom(J.in[k],1,psi)
      y<-Z
      
      ##fit model
      o1 <- try(optim(c(log(sig.in[k]), 0, 0, 0), nll.p, 
                      y.=y, Xs.=Xs, Xc.=Xc, hessian=TRUE,
                      method='BFGS', 
                      control = list(maxit=1000)))          
    }
    
    if (mod=='Nmix'){
      ##calculate exp. abundance, generate abundances
      lam<-exp(beta0+beta1*Xs+b.in[k]*Xcomp)
      N<-rpois(J.in[k],lam)
      
      ##generate detection data
      y<-matrix(rbinom(J.in[k]*K, rep(N, each=K), plogis(alpha0)),
                J.in[k], K, byrow=T)
      
      ##fit model
      o1 <- try(optim(c(log(sig.in[k]), 0, 0, 0,0), nll.nm,
                      y.=y, Xs.=Xs, Xc.=Xc,
                      hessian=TRUE,
                      method='BFGS',  
                      control = list(maxit=1000)))
    } 
    
    if (mod=='Occu'){
      
      ##calculate occu prob, generate occurrences
      psi<-plogis(beta0+beta1*Xs+b.in[k]*Xcomp)
      Z<-rbinom(J.in[k],1,psi)
      
      ##generate detection data
      y<-matrix(rbinom(J.in[k]*K, rep(Z, each=K), plogis(alpha0)),
                J.in[k], K, byrow=T)
      
      ##fit model
      o1 <- try(optim(c(log(sig.in[k]), 0, 0, 0, 0), nll.occ, 
                      y.=y, Xs.=Xs, Xc.=Xc,
                      hessian=TRUE,
                      method='BFGS', 
                      control = list(maxit=1000)))          
    }
    
    ##ensure rest of code runs if optim failed
    if(identical(class(o1), "try-error")) {
      o1 <- list()
      o1$par <- rep(NA, npar)
      o1$convergence <- -99
      se<-rep(NA, npar)
    } else {
      
      ##calculate SE, set to NA if cannot be calculated and mark as 
      ##not converged
      se<-try(sqrt(diag(solve(o1$hessian))))
      if(identical(class(se), "try-error")) {
        se<-rep(NA, npar)
        o1$convergence <- -98
      }
      
      
      ## if at least one SE cannot be calculates, set as not converged
      if(sum(is.na(se))>0) o1$convergence<- -98 
      ## If sigma is estimated too high (>15) for 10-yr time series, 
      ## set as not converged
      if(o1$par[1]>=log(15)) o1$convergence <- -97
      ## If beta.c is estimated way too high, set as not converged
      if(o1$par[4]>=10) o1$convergence <- -97
      ## If any SE is extreme, set to not converged
      if((sum(is.na(se)))==0 & (any(se>10))) o1$convergence <- -96
    }
    
    ##write results into matrix
    out[i,,k]<-o1$par
    se.out[i,,k]<-se
    conv[i,k]<-o1$convergence
    
    ##write data into matrix 
    if (mod %in% c('Occu', 'Nmix')) {ally[i,,]<-y} else{
      ally[i,]<-y
    }
    allX[i,,1:tt]<-Xc
    allX[i,,tt+1]<-Xcomp
    
  } ##sim loop
  ##save all data per scenario
  saveRDS(list(Xs=Xs, ally=ally, allX=allX), 
          paste(mod, '_Data_sig.', sig.in[k], '_J.', J.in[k],'_beta.', b.in[k], 
                '.rds', sep=''))
  print(k) ##keep track of where you are
} ##scenario loop

# ##write out results across all scenarios for given model type
saveRDS(list(out=out, se.out=se.out, conv=conv),
        paste(mod, '_Final_results.rds', sep=''))

