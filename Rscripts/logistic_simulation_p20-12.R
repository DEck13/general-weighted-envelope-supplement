
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


p <- 15; u <- 5
foo <- matrix(0, nrow = p, ncol = u)

for(j in 1:u){
	foo[(2*(j-1)+1):(2*j),j ] <- 1
}


O <-qr.Q(qr(foo), complete = TRUE)
Gamma <- O[, 1:u]
Gamma0 <- O[, (u+1):p]
beta <- rowSums(Gamma) / 7
as.matrix(beta, ncol = 1)

Omega <- diag( 2 * 1:5 ) 
Omega0 <- diag(exp( seq(from = -7, to = -1, length = 10) ))
SigmaX <- Gamma%*%Omega%*%t(Gamma) + Gamma0%*%Omega0%*%t(Gamma0)
eig <- eigen(SigmaX)
SigmaX.half <- eig$vec %*% diag(sqrt(eig$val)) %*% t(eig$vec)


# envelope simulation
set.seed(24)
nMC <- 10 # Monte Carlo sample size
ns <- c(400, 525, 675, 800) # sample sizes
nboot <- 500 # bootstrap sample size

system.time({
	lsims_MC_12 <- lapply(ns, function(j){
		do.call(rbind, lapply(lapply(1:nMC, function(x){  
			a <- proc.time()[3]
			out <- unpack(sim_function(
				n = j, 
				p = p, 
				nboot = nboot, 
				numCores = numCores, 
				SigmaX.half, 
				beta, 
				u = u, 
				family = "binomial"), 
				beta = beta)
			b <- proc.time()[3]
			cat("iter: ", x, " time: ", b-a, "\n")
			out
		}), function(xx) unlist(xx)))
	})
})


save(lsims_MC_12, file = "lsims20-12.RData")

