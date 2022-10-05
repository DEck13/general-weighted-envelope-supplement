rm(list = ls())

library(tidyverse)
library(TRES)
library(MASS)
library(parallel)
library(xtable)
library(faraway)

RNGkind("L'Ecuyer-CMRG")

source("source.R")
numCores

# envelope basis
p <- 6; u <- 3
init <- rep(1, 6)
v1 <- matrix(init, nrow = p)
O <-qr.Q(qr(v1), complete = TRUE)
Gamma <- O[, 1:u]
Gamma0 <- O[, (u+1):p]

# material and immaterial variation
Omega <- diag(exp(c(0,1,2))) 
Omega0 <- diag(exp(c(-4,-3,-1)))

# envelope setup
SigmaX <- Gamma%*%Omega%*%t(Gamma) + Gamma0%*%Omega0%*%t(Gamma0)
eig <- eigen(SigmaX)
SigmaX.half <- eig$vec %*% diag(sqrt(eig$val)) %*% t(eig$vec)
beta <- -(5*Gamma[, 1] + 5*Gamma[, 2] + 2*Gamma[, 3]) / sqrt(54)
as.matrix(beta, ncol = 1)



# envelope simulation
set.seed(13)
nMC <- 100
ns <- c(200, 400, 600, 800)
nboot <- 1e3
system.time({
	psims_MC_B <- lapply(ns, function(j){
		do.call(rbind, lapply(lapply(1:nMC, function(x){  
			unpack(sim_function(n = j, p = p, nboot = nboot, 
													numCores = numCores, SigmaX.half, beta, u = u, 
													family = "poisson"), beta = beta)
		}), function(xx) unlist(xx)))
	})
})

save(psims_MC_B, file = "psimsB.RData")

