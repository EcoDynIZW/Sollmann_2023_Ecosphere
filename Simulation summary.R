###############################################################################
### summarize simulation results for time scale model #########################

## packages for plots only
library(reshape2)
library(ggplot2)

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

#number of simulations per scenario, ran up to 100
nsim=100

# manually set model type for which to summarize simulation, out of
# 'Poisson', 'LogReg', 'Occu', 'Nmix'

mod<-'Occu'

##set number of parameters estimated depending on model type
npar<-ifelse(mod %in% c('Poisson', 'LogReg'),4,5)

## read in simulation results
res<-readRDS(paste(mod, '_Final_results.rds', sep=''))

###############################################################################

## sum non-converged iterations across iterations for each scenario
non.c<-apply(res$conv, 2, function(x){sum(x!=0)})/nsim


## plot non-convergence rates against levels of input parameters
#pdf(paste('Plotconvergence.', mod, '.pdf', sep=''))
par(mfrow=c(2,2))
##against beta - strongest
plot(b.in, non.c, ylab='Rate of non-convergence',
        xlab='Beta(composite X)', pch=19, axes=F, ylim=c(0, 1))
axis(side=2)
axis(side=1, at=beta.c)
box()
mean.c<-apply(matrix(non.c, nrow=9, ncol=3), 2, median)
points(beta.c, mean.c, pch=19, cex=1.3, col='blue')
points(beta.c, mean.c, type='l', col='blue')


##against sigma
plot(sig.in, non.c, ylab='Rate of non-convergence',
     xlab='Sigma', pch=19, axes=F, ylim=c(0, 1))
axis(side=2)
axis(side=1, at=sigma)
box()
mean.c<-apply(matrix(non.c, nrow=9, ncol=3, byrow=T), 2, median)
points(sigma, mean.c, pch=19, cex=1.3, col='blue')
points(sigma, mean.c, type='l', col='blue')


##against J - weakest
plot(J.in, non.c, ylab='Rate of non-convergence',
     xlab='Sample size', pch=19, axes=F, ylim=c(0, 1))
axis(side=2)
axis(side=1, at=J)
box()
mean.c<-c(median(non.c[J.in==100]), median(non.c[J.in==175]), median(non.c[J.in==250]))
points(J, mean.c, pch=19, cex=1.3, col='blue')
points(J, mean.c, type='l', col='blue')
#dev.off()


################################################################################
## calculate summary info for parameters: absolute bias, SE
## as point summary, use median, because of skewed distributions

resmat<-matrix(NA, 0, 7)
colnames(resmat)<-c('Scenario', 'MeanEst', 'Avg.SE', 'Avg.Bias', 'RMSE', 'Truth',
                    'N.conv')

for (k in 1:n.scen){
  if(npar==5){
  truth<-c(log(sig.in[k]), beta0, beta1, b.in[k], alpha0)}else {
    truth<-c(log(sig.in[k]), beta0, beta1, b.in[k]) } 
  sub.c<-res$conv[,k]
  sub<-res$out[sub.c==0,,k]
  sub.e<-res$se.out[sub.c==0,,k]
  ressub<-cbind(rep(k, npar),
                apply(sub, 2, median),
                apply(sub.e, 2, median),
                apply(sub-matrix(truth, nrow(sub), ncol(sub), byrow=T), 2, median),
                sqrt(apply((sub-matrix(truth, nrow(sub), ncol(sub), byrow=T))^2,2,mean)),
                truth,
                rep(dim(sub)[1], npar) )
  resmat<-rbind(resmat, ressub)
}


## add parameter names, turn into dataframe
if (npar==5){
pps<-c('log(sig)', 'beta0', 'beta1', 'beta.c', 'alpha0')
} else {pps<-c('log(sig)', 'beta0', 'beta1', 'beta.c')}

results<-as.data.frame(resmat)
results$Parameter<-rep(pps, 27)
results<-results[,c(1,8,6,2,3,4,5,7)]
results[,3:7]<-round(results[,3:7], dig=3)

##write out as csv
##Note: formatted (dec, sep) for German computer
write.table(results, file=paste('Results.',mod,'.csv', sep=''), sep=';', col.names = TRUE,
            dec=',', row.names = F)

##write out as R file
saveRDS(results, paste('Results.',mod,'.rds', sep=''))


################################################################################
####  make main figures ########################################################

#### Requires running simulations for all model types first, or 
#### adjust to only model types of interest


modvec<-c('Poisson', 'LogReg', 'Nmix', 'Occu') 

##first, 'flatten' results

out.list<-list()

for (ii in 1:length(modvec)){
  if (modvec[ii] %in% c('Poisson', 'LogReg')) npar<-4
  else npar<-5
  ## read in results
  res<-readRDS(paste(modvec[ii], '_Final_results.rds', sep=''))
  
  ##set non-converged to NA so they are not plotted
  res.conv<-res$out
  for (i in 1:nsim){
    for (j in 1:n.scen){
      if(res$conv[i,j]==0)next
      res.conv[i,,j]<-NA
    }
  }
  
  ##parameter order: sigma is 1, beta.c is 4
  out.flat<-melt(res.conv, varnames = c('Iteration', 'Parameter', 'Scenario'))
  out.flat$sigma<-rep(sig.in, each=npar*nsim) 
  out.flat$beta.c<-rep(b.in, each=npar*nsim)
  out.flat$J<-rep(J.in, each=npar*nsim)
  out.flat$Model<-rep(mod[ii], dim(out.flat)[1])
  
  
  out.list[[ii]]<-out.flat
}

all.out<-do.call(rbind, out.list)


################ plot results for sigma ########################################

##truth, sigma
v_line <- data.frame(
  yintercept = rep( log(rep(c(1, 2.5, 4), each=3)), 4),
  beta.c = rep( rep(beta.c, 3), 4),
  sigma=rep( rep(sigma, each=3), 4),
  Model=rep( rep(mod, each = 3), 3)
)



##row and column names for panel plot
sig_names <- c(
  `1` = 'sigma = 1',
  `2.5` = 'sigma = 2.5',
  `4` = 'sigma = 4'
)

b_names <- c(
  `0.5` = 'b.c = 0.5',
  `0.9` = 'b.c = 0.9',
  `1.3` = 'b.c = 1.3'
)

m_names <- c(
  `LogReg` = 'Logistic',
  `Nmix` = 'N-mixture',
  `Occu` = 'Occupancy',
  `Poisson`= 'Poisson'
)

all.out2<-all.out
all.out2$Model<-factor(all.out2$Model,
                       levels=c('Poisson', 'Nmix', 'LogReg', 'Occu'))

##remove extreme values per response variable, otherwise yaxis makes plot useless
ab<-all.out[all.out$Model %in% c('Poisson', 'Nmix') & all.out$Parameter == 1,]
ab<-ab[!is.na(ab$value),]
mima<-quantile(ab$value, p=c(0.005, 0.995), na.rm=T)
ab1<-ab[ab$value>mima[1] & ab$value<mima[2],]

oc<-all.out[all.out$Model %in% c('LogReg', 'Occu') & all.out$Parameter == 1,]
oc<-oc[!is.na(oc$value),]
mima<-quantile(oc$value, p=c(0.005, 0.995), na.rm=T)
oc1<-oc[oc$value>mima[1] & oc$value<mima[2],]

##put abundance and occurrence results back together
aboc.or<-rbind(ab, oc)
aboc<-rbind(ab1, oc1)

##This corresponds to Figure 1
ggplot(aboc, aes(x = J, y = value)) +
  geom_violin(color = "black",  aes(group=J), fill='grey97') + # alpha = .3,
  labs(x = "J", y = "Estimate log(sigma)") +
  facet_nested(sigma ~ factor(Model, levels = c('Poisson', 'Nmix', 'LogReg', 'Occu'))+beta.c,
               labeller = labeller(sigma  = as_labeller(sig_names),
                                   beta.c = as_labeller(b_names),
                                   Model = as_labeller(m_names))
  ) +
  geom_hline(data=v_line,aes(yintercept = yintercept), 
             linetype='dashed', color = 'orangered', size=0.5)+
  scale_x_continuous(breaks=c(100, 175, 250)) +
  stat_summary(data= aboc.or,
               fun = "median",
               geom = "point",
               color = "black",
               size = 0.75)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(panel.spacing.x	= unit(0.1, "lines"))


############# plot results for beta ############################################

##truth, beta
v_line2 <- data.frame(
  yintercept = rep(rep(beta.c, 3),2),
  beta.c = rep(rep(beta.c, 3),2),
  sigma=rep(rep(sigma, each=3),2),
  Model=rep( rep(c('Poisson', 'Nmix'), each = 3), 3)
)

v_line2b <- data.frame(
  yintercept = rep(rep(beta.c, 3),2),
  beta.c = rep(rep(beta.c, 3),2),
  sigma=rep(rep(sigma, each=3),2),
  Model=rep( rep(c('LogReg', 'Occu'), each = 3), 3)
)

## needs two plots with different y scales
## plot estimates of beta.c
## for better visualization, cut off top and bottom 0.5% of data points
sub1<-all.out[all.out$Parameter==4 & 
                all.out$Model %in% c('Poisson', 'Nmix'),]
sub1<-sub1[!is.na(sub1$value),]
mima<-quantile(sub1$value, p=c(0.005, 0.995), na.rm=T)
sub1.sub<-sub1[sub1$value>mima[1] & sub1$value<mima[2],]


g1<-ggplot(sub1.sub, aes(x = J, y = value)) +
  geom_violin(color = "black",  aes(group=J), fill='grey97') + # alpha = .3,
  labs(x = "J", y = "Estimate b.c") +
  facet_nested(sigma ~ factor(Model, levels = c('Poisson', 'Nmix'))+beta.c,
               labeller = labeller(sigma  = as_labeller(sig_names),
                                   beta.c = as_labeller(b_names),
                                   Model = as_labeller(m_names))
  ) +
  geom_hline(data=v_line2,aes(yintercept = yintercept), 
             linetype='dashed', color = 'orangered', size=0.5)+
  scale_x_continuous(breaks=c(100, 175, 250)) +
  stat_summary(data = sub1,
               fun = "median",
               geom = "point",
               color = "black",
               size = 0.75)+
  theme_bw()+
  theme(panel.spacing.x	= unit(0.1, "lines"))+
  theme(plot.margin = unit(c(5.5,.2,5.5,5.5), "points"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 


sub2<-all.out[all.out$Parameter==4 & 
                all.out$Model %in% c('LogReg', 'Occu'),]
sub2<-sub2[!is.na(sub2$value),]
mima<-quantile(sub2$value, p=c(0.005, 0.995), na.rm=T)
sub2.sub<-sub2[sub2$value>mima[1] & sub2$value<mima[2],]

g2<-ggplot(sub2.sub, aes(x = J, y = value)) +
  geom_violin(color = "black", aes(group=J), fill='grey97') + # alpha = .3, 
  labs(x = "J", y='') +
  facet_nested(sigma ~ factor(Model, levels = c('LogReg', 'Occu'))+beta.c,
               labeller = labeller(sigma  = as_labeller(sig_names),
                                   beta.c = as_labeller(b_names),
                                   Model = as_labeller(m_names))
  ) +
  
  geom_hline(data=v_line2b,aes(yintercept = yintercept), 
             linetype='dashed', color = 'orangered', size=0.5)+
  scale_x_continuous(breaks=c(100, 175, 250)) +
  stat_summary(data = sub2,
               fun = "median",
               geom = "point",
               color = "black",
               size = 0.75)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(panel.spacing.x	= unit(0.1, "lines"))+
  theme(plot.margin = unit(c(5.5,5.5,5.5,0.1), "points"))

##Corresponds to Figure 2
gall<-grid.arrange(g1,g2, nrow = 1, ncol=2)
