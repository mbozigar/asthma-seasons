
rm(list = ls()) 
getwd()
#setwd()

options(scipen = 999) # disables scientific notation in favor of decimal format

library(mcmcplots)
library(nimble)
library(coda)


# model running info#########################################################################################################
niter=50000; nburnin=45000; thin=5; nchains=1;
sample.total=((niter-nburnin)/thin)*nchains; sample.total
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


# Model definition #########################################################################################################
model.code<-nimbleCode({
  for (i in 1:m) # index ED visit records
  {	
    for (j in 1:4) # index case day (1) and control days (3)
    {
      y[i,j] ~ dbern(pi[i,j]) # data model for binary outcome ED visit yes/no
      logit(pi[i,j]) <- a0 + 
        #b1*MCOc[i,j] + 
        b2*MNOxc[i,j] + 
        b3*MO3c[i,j] + 
        b4*MPMxc[i,j] + 
        b5*MPM25c[i,j] +
        b6*MSO2c[i,j] + 
        b7*Mtempc[i,j] + b8*Mdewpc[i,j] + b9*Mtedec[i,j] +
        c[combos[i]] 

      LL[i,j] <- y[i,j]*log(pi[i,j]) + (1-y[i,j])*log(1-pi[i,j])
      dev[i,j] <- -2*LL[i,j]
    }
    sdevc[i] <- sum(dev[i,1:4]) # sum of deviances for each case/control day
  }
  sdevo <- sum(sdevc[1:m]) # sum of deviances (overall deviance) for each record
  
  for (l in 1:q)
  {
    c[l] ~ dnorm(0,tauc)
  }
  
  a0~dnorm(0,taua)  #prior
  taua<-1/pow(sda,2)
  sda~dunif(0,4)  #hyperprior
  #b1~dnorm(0,taub1); taub1~dgamma(2,1) # CO
  b2~dnorm(0,taub2); taub2~dgamma(2,1) # NOx
  b3~dnorm(0,taub7); taub3~dgamma(2,1) # O3
  b4~dnorm(0,taub4); taub4~dgamma(2,1) # PMx
  b5~dnorm(0,taub5); taub5~dgamma(2,1) # PM2.5
  b6~dnorm(0,taub6); taub6~dgamma(2,1) # SO2
  b7~dnorm(0,taub7); taub7~dgamma(2,1) # temp
  b8~dnorm(0,taub8); taub8~dgamma(2,1) # dewp
  b9~dnorm(0,taub9); taub9~dgamma(2,1) # temp*dewp
  tauc ~ dgamma(2,0.5)
})
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


# load ##################################################################################################################
load(file = "path/constants_SC_winter.rda") # constants/data

load(file = "path/outcome_SC_winter.rda") # outcome

load(file = "path/inits_SC_winter.rda") # initial values
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


# nimble model code, compile, onfigure, set monitors, and build ##########################################################
# create project
begin <- Sys.time()
model.n<-nimbleModel(code=model.code,name="model.n", constants=Rdata, data=Routcome, inits=inits) 
# compile
model.n.comp<-compileNimble(model.n) 
# configure and set monitors ()
model.n.conf<-configureMCMC(model.n,print=TRUE)
# add monitors
model.params <- c("a0", "sda", 
                  #"b1", "taub1", # CO
                  "b2", "taub2", # NOx
                  "b3", "taub3", # O3
                  "b4", "taub4", # PMx
                  "b5", "taub5", # PM2.5
                  "b6", "taub6", # SO2
                  "b7", "b8","b9",
                  "taub7", "taub8","taub9",
                  "tauc",
                  "sdevo") 
model.n.conf$addMonitors(model.params)
# build
model.n.build<-buildMCMC(model.n.conf) 
# recompile 
model.n.build.re<-compileNimble(model.n.build,project=model.n.comp)
end <- Sys.time(); elapsed<-end-begin; elapsed
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


# run model #############################################################################################
set.seed(1)
begin <- Sys.time()
# run the model (USE RECOMPILED FOR WAIC)
samples<-runMCMC(model.n.build.re, niter=niter, nburnin=nburnin, nchains=nchains, summary=TRUE, samplesAsCodaMCMC=TRUE, thin=thin) # WAIC = TRUE, nburnin=100000, thin=50
end <- Sys.time(); elapsed<-end-begin; elapsed
# &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


# plots, diagnostics, and results #################################################################################
library(Rcpp)
library(plyr)
library(dplyr)
library(tidyverse)
library(ggplot2)

variablenames<-colnames(samples$samples)
newvector <- variablenames[1:length(variablenames)]
sampledata <- list() 
sampledata <- samples$samples[,newvector]
# sampledata[[1]] <- samples$samples[[1]][,newvector] # chain 1
# sampledata[[2]] <- samples$samples[[2]][,newvector] # chain 2
sampledata <- as.mcmc(sampledata)
# sampledata <- as.mcmc.list(sampledata)

options(max.print=8000) # For big lists
begin <- Sys.time()
results <- summary(sampledata) # create a list from the summary output
end <- Sys.time(); elapsed<-end-begin; elapsed
results[] <- lapply(results,round,5) # round

# estimates
summary(sampledata)
# convergence
geweke.diag(sampledata[,model.params])
denplot(sampledata, parms=model.params)

# trace plot
traplot(sampledata, parms=model.params)




