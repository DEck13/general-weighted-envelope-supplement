
rm(list = ls())
library(tidyverse)
library(TRES)
library(MASS)
library(parallel)
library(xtable)
library(faraway)
library(ggplot2)

numCores <- 10
RNGkind("L'Ecuyer-CMRG")

##################################################
#         1D optimization function               #
##################################################
fun1D <- function(W, M, U){
	f <- log(t(W) %*% M %*% W) + log(t(W) %*% chol2inv(chol(M+U)) %*% W)
	df <- 2*(M %*% W/(as.numeric(t(W) %*% M %*% W))+
					 	solve(M+U) %*% W/(as.numeric(t(W) %*% chol2inv(chol(M+U)) %*% W)))
	list(F = f, G = df)
}

##################################################
#    get initial value for 1D algorithm          #
##################################################
get_ini1D <- function(M, U){
	p <- dim(U)[2]
	v1 <- eigen(M)$vectors
	v2 <- eigen(M+U)$vectors
	v <- cbind(v1, v2)
	index <- v[1,] < 0
	v[,index] <- -v[,index]
	W0 <- Re(v[, 1]) ## Ensure the real number
	Fw0 <- fun1D(W0, M, U)$F
	for (i in 2:(2*p)) {
		W <- Re(v[, i])
		Fw <- fun1D(W, M, U)$F
		if (Fw < Fw0) {
			W0 <- W
			Fw0 <- Fw
		}
	}
	W0
}


################################################################
#         1D solver for individual objective function          #
#         using function OptManiMulitBallGBB                   #
################################################################
ballGBB1D <- function(M, U, ...) {
	
	# Options for function OptManiMulitBallGBB
	opts <- list(...)
	W0 <- get_ini1D(M, U)
	if (is.null(opts$xtol) || opts$xtol < 0 || opts$xtol > 1) opts$xtol <- 1e-8
	if (is.null(opts$gtol) || opts$gtol < 0 || opts$gtol > 1) opts$gtol <- 1e-8
	if (is.null(opts$ftol) || opts$ftol < 0 || opts$ftol > 1) opts$ftol <- 1e-12
	# parameters for control the linear approximation in line search
	if (is.null(opts$rho) || opts$rho < 0 || opts$rho > 1) opts$rho <- 1e-04
	# factor for decreasing the step size in the backtracking line search
	if (is.null(opts$eta))
		opts$eta <- 0.2 else if (opts$eta < 0 || opts$eta > 1)
			opts$eta <- 0.1
	# parameters for updating C by HongChao, Zhang
	if (is.null(opts$gamma) || opts$gamma < 0 || opts$gamma > 1) opts$gamma <- 0.85
	if (is.null(opts$tau) || opts$tau < 0 || opts$tau > 1) opts$tau <- 1e-03
	# parameters for the  nonmontone line search by Raydan
	if (is.null(opts$m)) opts$m <- 10
	if (is.null(opts$STPEPS)) opts$STPEPS <- 1e-10
	if (is.null(opts$maxiter) || opts$maxiter < 0 || opts$maxiter > 2^20) opts$maxiter <- 800
	if (is.null(opts$nt) || opts$nt < 0 || opts$nt > 100) opts$nt <- 5
	# if (is.null(opts$record)) opts$record <- 0
	if (is.null(opts$eps)) opts$eps <- 1e-14
	
	fit <- OptManiMulitBallGBB(W0, opts, fun1D, M, U)
	list(X = fit$X, out = fit$out)
}

###########################################################################
#         Line search algorithm for optimization on manifold              #
###########################################################################
OptManiMulitBallGBB <- function(X, opts=NULL, fun, ...) {
	# Line search algorithm for optimization on manifold:
	#
	#    min f(X), s.t., ||X_i||_2 = 1, where X \in R^{n,p}
	#        g(X) = grad f(X)
	#    X = [X_1, X_2, ..., X_p]
	#    each column of X lies on a unit sphere
	
	# Input:
	#    X --- initialization. ||X_i||_2 = 1, each column of X lies on a unit sphere
	#    opts --- option structure with fields:
	#       maxiter     max number of iterations
	#       xtol        stop control for ||X_k - X_{k-1}||
	#       gtol        stop control for the projected gradient
	#       ftol        stop control for |F_k - F_{k-1}|/(1+|F_{k-1}|)
	#       usually, max{xtol, gtol} > ftol
	#    fun --- objective function and its gradient:
	#    [F, G] = fun(X, data1, data2)
	#
	#    F, G are the objective function value and gradient, repectively
	#    data1, data2 are addtional data, and can be more.
	#
	# Calling syntax:
	#    OptManiMulitBallGBB(X0, opts, fun, data1, data2);
	#
	# Output:
	#     x --- solution
	#     g --- gradient of x
	#     Out --- output information
	#
	# For example, consider the maxcut SDP:
	#     X is n by n matrix
	#     max Tr(C*X), s.t., X_ii = 1, X is PSD.
	#
	#     low rank model is:
	#         X = V'*V, V = [V_1, ..., V_n], V is a p by n matrix
	#         max Tr(C*V'*V), s.t., ||V_i|| = 1,
	#
	#     # Define the function returning objective and the gradient:
	#     maxcut_quad <- function(V, C){
	#       g = 2*(V*C)
	#       f = sum(dot(g,V))/2
	#       return(list(F = f, G = g))
	#     }
	#
	#     # Call function OptManiMulitBallGBB
	#     OptManiMulitBallGBB(x0, opts, maxcut_quad, C);
	#
	#
	#
	#     Reference: Z. Wen and W. Yin (2013), A feasible method for optimization
	#       with orthogonality constraints
	#
	
	##Size information
	X <- as.matrix(X)
	if (length(X) == 0)
		print("input X is an empty") else {
			n <- dim(X)[1]; k <- dim(X)[2]}
	
	if (is.null(opts$xtol) || opts$xtol < 0 || opts$xtol > 1) opts$xtol <- 1e-8
	if (is.null(opts$gtol) || opts$gtol < 0 || opts$gtol > 1) opts$gtol <- 1e-8
	if (is.null(opts$ftol) || opts$ftol < 0 || opts$ftol > 1) opts$ftol <- 1e-12
	# parameters for control the linear approximation in line search
	if (is.null(opts$rho) || opts$rho < 0 || opts$rho > 1) opts$rho <- 1e-04
	# factor for decreasing the step size in the backtracking line search
	if (is.null(opts$eta))
		opts$eta <- 0.2 else if (opts$eta < 0 || opts$eta > 1)
			opts$eta <- 0.1
		# parameters for updating C by HongChao, Zhang
		if (is.null(opts$gamma) || opts$gamma < 0 || opts$gamma > 1) opts$gamma <- 0.85
		if (is.null(opts$tau) || opts$tau < 0 || opts$tau > 1) opts$tau <- 1e-03
		# parameters for the  nonmontone line search by Raydan
		if (is.null(opts$m)) opts$m <- 10
		if (is.null(opts$STPEPS)) opts$STPEPS <- 1e-10
		if (is.null(opts$maxiter) || opts$maxiter < 0 || opts$maxiter > 2^20) opts$maxiter <- 800
		if (is.null(opts$nt) || opts$nt < 0 || opts$nt > 100) opts$nt <- 5
		if (is.null(opts$eps)) opts$eps <- 1e-14
		if (is.null(opts$record)) opts$record <- 0
		
		# copy parameters
		xtol <- opts$xtol
		gtol <- opts$gtol
		ftol <- opts$ftol
		rho  <- opts$rho
		m <- opts$m
		STPEPS <- opts$STPEPS
		eta <- opts$eta
		gamma <- opts$gamma
		eps <- opts$eps
		nt <- opts$nt
		crit <- matrix(1, opts$maxiter, 3)
		record <- opts$record
		
		# normalize x so that ||x||_2 = 1
		nrmX <- apply(X*X, 2, sum)
		nrmX <- matrix(nrmX, 1, k)
		if (sqrt(sum((nrmX-1)^2)) > 1e-8) {
			X <- sweep(X, 2, sqrt(nrmX),"/")
		}
		args = list(X, ...)
		eva <- do.call(fun, args)
		f <- eva$F; g <- as.matrix(eva$G)
		out <- c()
		out$nfe <- 1
		
		
		Xtg <- apply(X*g, 2, sum)
		Xtg <- matrix(Xtg, 1, k)
		gg <- apply(g*g, 2, sum)
		gg <- matrix(gg, 1, k)
		XX <- apply(X*X, 2, sum)
		XX <- matrix(XX, 1, k)
		XXgg <- XX*gg
		temp <- sweep(X, 2, Xtg, "*")
		dtX <- matrix(temp, n, k) - g
		nrmG <- sqrt(sum((dtX)^2))
		
		Q <- 1; Cval <- f; tau <- opts$tau
		
		## print iteration header if record >= 1
		if (record >= 1){
			cat(paste("\n", '------ Gradient Method with Line search ----- ',"\n"),
					sprintf("%4s %8s %10s %9s %9s %3s", 'Iter', 'tau', 'F(X)', 'nrmG', 'XDiff', 'nls'), "\n")
		}
		if (record == 10) out$fvec = f
		
		# ##
		# X_list <- vector("list", opts$maxiter)
		# XDiff_list <- vector("list", opts$maxiter)
		# FDiff_list <- vector("list", opts$maxiter)
		# Q_list <- vector("list", opts$maxiter)
		# Cval_list <- vector("list", opts$maxiter)
		# ##
		
		##main iteration
		for (itr in 1:opts$maxiter) {
			Xp <- X; fp <- f; gp <- g; dtXP <- dtX
			
			nls <- 1; deriv = rho*nrmG^2
			
			while (TRUE) {
				## calculate g, f
				tau2 <- tau/2
				beta <- (1 + (tau2^2)*(-(Xtg^2) + XXgg))
				a1 <- ((1 + tau2*Xtg)^2 - (tau2^2)*XXgg)/beta
				a2 <- -tau*XX/beta
				X <- sweep(Xp, 2, a1, "*") + sweep(gp, 2, a2, "*")
				
				
				args = list(X, ...)
				eva <- do.call(fun, args)
				f <- eva$F; g <- as.matrix(eva$G)
				out$nfe <- out$nfe + 1
				
				if (f <= Cval - tau*deriv || nls >= 5)
					break
				tau <- eta * tau
				nls <- nls + 1
			}
			
			if (record == 10)
				out$fvec <- rbind(out$fvec,f)
			
			Xtg <- apply(X*g, 2, sum)
			Xtg <- matrix(Xtg, 1, k)
			gg <- apply(g*g, 2, sum)
			gg <- matrix(gg, 1, k)
			XX <- apply(X*X, 2, sum)
			XX <- matrix(XX, 1, k)
			XXgg <- XX*gg
			temp <- sweep(X, 2, Xtg, "*")
			dtX <- matrix(temp, n, k) - g
			
			nrmG <- sqrt(sum(dtX^2))
			s <- X - Xp
			XDiff <- sqrt(sum(s^2))/sqrt(n)
			FDiff <- abs(fp - f)/(abs(fp) + 1)
			
			if (record >= 1)
				# cat(paste('------ Gradient Method with Line search ----- ',"\n"),
				#     sprintf("%4s %8s %8s %10s %10s %10s", 'Iter', 'tau', 'F(X)', 'nrmG', 'XDiff','nls'))
				cat(sprintf('%4d  %3.2e  %3.2e  %3.2e  %3.2e %2d',
										itr, tau, f, nrmG, XDiff, nls), "\n")
			
			crit[itr,] <- cbind(nrmG, XDiff, FDiff)
			r <- length((itr - min(nt, itr)+1):itr)
			temp1 <- matrix(crit[(itr - min(nt, itr)+1):itr,], r, 3)
			mcrit <- apply(temp1, 2, mean)
			
			
			if ((XDiff < xtol && FDiff < ftol) || nrmG < gtol
					|| all(mcrit[2:3] < 10 * c(xtol, ftol))) {
				out$msg = "converge"
				break
			}
			
			y <- dtX - dtXP
			sy <- sum(s*y)
			tau <- opts$tau
			sy <- abs(sy)
			
			if (sy > 0) {
				if (itr %% 2 == 0)
					tau <- sum(s*s)/sy else { tau <- sy/sum(y*y)}
				
				tau <- max(min(tau, 1e+20), 1e-20)
			}
			Qp <- Q; Q <- gamma*Qp + 1
			Cval <- (gamma*Qp*Cval + f)/Q
		}
		
		if (itr >= opts$maxiter)
			out$msg <- "exceed max iteration"
		
		Xn <- apply(X*X, 2, sum)
		Xn <- matrix(Xn, 1, k)
		out$feasi <- svd(Xn - 1)$d[1]
		
		if (out$feasi > eps) {
			nrmX <- apply(X*X, 2, sum)
			X <- sweep(X, 2, sqrt(nrmX),"/")
			args = list(X, ...)
			eva <- do.call(fun, args)
			f <- eva$F; g <- as.matrix(eva$G)
			out$nfe <- out$nfe + 1
			nrmX.n <- apply(X*X, 2, sum)
			out$feasi <- svd(nrmX.n - 1)$d[1]
		}
		
		out$nrmG <- nrmG
		out$fval <- f
		out$itr <- itr
		
		return(list(X = X, g = g, out = out))
}


## new oneD_bic function which performs manifold1D 
## and oneD_bic from TRES
oneD_bic_mod <- function(M=M,U=U,C=1,n=n,maxdim=10){
	p <- dim(M)[2]
	Mnew <- M
	Unew <- U
	G <- matrix(0, p, maxdim)
	G0 <- diag(p)
	phi <- rep(0, p)
	for (k in 1:maxdim) {
		if (k == p) 
			break
		gk <- ballGBB1D(Mnew, Unew)$X
		phi[k] <- n * (log(t(gk) %*% Mnew %*% gk) + 
									 	log(t(gk) %*% chol2inv(chol(Mnew + Unew)) %*% gk)) + 
			log(n) * C
		G[, k] <- G0 %*% gk
		G0 <- qr.Q(qr(G[, 1:k]), complete = TRUE)[, (k + 1):p]
		Mnew <- t(G0) %*% M %*% G0
		Unew <- t(G0) %*% U %*% G0
	}
	bicval <- rep(0, maxdim)
	for (k in 1:maxdim) {
		bicval[k] <- sum(phi[1:k])
	}
	u <- which.min(bicval)
	Gamma_full <- G
	return(list(u = u, bicval = bicval, Gamma = Gamma_full))
}

oneD_bic <- function(M, U, n, C=1, maxdim=10) {
	p <- dim(M)[2]
	Mnew <- M
	Unew <- U
	G <- matrix(0, p, maxdim)
	G0 <- diag(p)
	phi <- rep(0, p)
	for (k in 1:maxdim){
		if (k == p) break
		gk <- ballGBB1D(Mnew, Unew)$X
		phi[k] <- n*(log(t(gk) %*% Mnew %*% gk)+ log(t(gk) %*% chol2inv(chol(Mnew+Unew)) %*% gk))+
			log(n)*C
		G[, k] <- G0 %*% gk
		G0 <- qr.Q(qr(G[, 1:k]), complete=TRUE)[, (k+1):p]
		Mnew <- t(G0) %*% M %*% G0
		Unew <- t(G0) %*% U %*% G0
	}
	bicval <- rep(0, maxdim)
	for (k in 1:maxdim) {
		bicval[k] <- sum(phi[1:k])
	}
	u <- which.min(bicval)
	Gamma <- G[,1:u]
	return(list(u = u, bicval = bicval, Gamma = Gamma))
}

## compute M and U with normal predictors for logistic
## regression model with normal predictors
Logistic_cov <- function(Y,X,a,b){
	n <- nrow(X); p <- ncol(X)
	theta <- a + X%*%b
	wts <- as.numeric(exp(theta)/((1 + exp(theta))^2))
	wts <- wts/mean(wts)
	Exw = t(wts) %*% X / n 
	Sxw <- t(X-do.call(rbind, replicate(n, Exw, simplify=FALSE))) %*% 
		diag(wts) %*% (X-do.call(rbind, replicate(n, Exw, simplify=FALSE))) / n
	Ys <- theta + (Y-exp(theta)/(1+exp(theta)))/wts
	Eyw <- t(wts) %*% Ys / n
	Sxyw <- t(X-do.call(rbind, replicate(n, Exw, simplify=FALSE))) %*% 
		diag(wts) %*% (Ys-do.call(rbind, replicate(n, Eyw, simplify=FALSE))) / n
	Syw = t(Ys-do.call(rbind, replicate(n, Eyw, simplify=FALSE))) %*% 
		diag(wts) %*% (Ys-do.call(rbind, replicate(n, Eyw, simplify=FALSE))) / n
	M <- Sxw
	U <-  Sxyw%*%t(Sxyw) / as.numeric(Syw)
	out = list(M = M, U = U)
	return(out)
}

## compute M and U with normal predictors for poisson 
## regression model with normal predictors
Poisson_cov <- function(Y,X,a,b){
	n <- nrow(X); p <- ncol(X)
	theta <- a + X%*%b
	wts <- as.numeric(exp(theta))
	wts <- wts/mean(wts)
	Exw = t(wts) %*% X / n 
	Sxw <- t(X-do.call(rbind, replicate(n, Exw, simplify=FALSE))) %*% 
		diag(wts) %*% (X-do.call(rbind, replicate(n, Exw, simplify=FALSE))) / n
	Ys <- theta + (Y-exp(theta))/wts
	Eyw <- t(wts) %*% Ys / n
	Sxyw <- t(X-do.call(rbind, replicate(n, Exw, simplify=FALSE))) %*% 
		diag(wts) %*% (Ys-do.call(rbind, replicate(n, Eyw, simplify=FALSE))) / n
	Syw = t(Ys-do.call(rbind, replicate(n, Eyw, simplify=FALSE))) %*% 
		diag(wts) %*% (Ys-do.call(rbind, replicate(n, Eyw, simplify=FALSE))) / n
	M <- Sxw
	U <-  Sxyw%*%t(Sxyw) / as.numeric(Syw)
	out = list(M = M, U = U)
	return(out)
}

## weighted envelope estimation using the 1D algorithm
weighted_env1D <- function(M, U, n, betahat){
	p <- length(betahat)
	oneD <- oneD_bic_mod(M = M, U = U, n = n, C = 1, maxdim = p)
	bic_val <- oneD$bicval
	min_bic <- min(bic_val)
	w <- exp(min_bic - bic_val) / sum(exp(min_bic - bic_val))
	rowSums(do.call(cbind, lapply(1:p, function(k){
		if(k < p){
			w[k] * tcrossprod(oneD$Gamma[, 1:k]) %*% betahat	
		}
		else{
			w[p] * betahat
		}
	})))
}

## manifoldFG
mod <- Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
mani.params <- get.manifold.params(IsCheckParams = TRUE)

FGfun_mfd <- function(M, U, u) {
	n <- dim(M)[2]
	mw <- function(w) { matrix(w, n, u) }
	f <- function(w) { W <- mw(w); log(det(t(W) %*% M %*% W)) + log(det(t(W) %*% chol2inv(chol(M+U)) %*% W))  }
	g <- function(w) { W <- mw(w); 2*(M %*% W %*% chol2inv(chol(t(W) %*% M %*% W))+ chol2inv(chol(M+U)) %*% W %*% chol2inv(chol((t(W) %*% chol2inv(chol(M + U)) %*% W))) ) }
	
	prob <- methods::new(mod$RProblem, f, g)
	mani.defn <- ManifoldOptim::get.stiefel.defn(n, u)
	
	return(list(prob=prob, mani.defn=mani.defn))
}

manifoldFG <- function(M, U, u, Gamma_init = NULL){
	if (is.null(Gamma_init)) {
		if (missing(u)) 
			stop("u is required since Gamma_init is NULL.")
		Gamma_init <- manifold1D(M, U, u)
	}
	p <- dim(Gamma_init)[1]
	u <- dim(Gamma_init)[2]
	opts <- list()
	if (is.null(opts$maxiter)) 
		opts$maxiter = 500
	else if (opts$maxiter < 0 || opts$maxiter > 2^20) 
		opts$maxiter = 500
	if (is.null(opts$tol)) 
		opts$tol = 1e-08
	else if (opts$tol < 0 || opts$tol > +1) 
		opts$tol = 1e-08
	if (is.null(opts$method)) 
		opts$method = "RBFGS"
	if (is.null(opts$check)) 
		opts$check = FALSE
	if (dim(U)[1] != dim(U)[2]) {
		U = tcrossprod(U)
	}
	mani.params <- get.manifold.params(IsCheckParams = opts$check)
	solver.params <- ManifoldOptim::get.solver.params(Max_Iteration = opts$maxiter, 
																										Tolerance = opts$tol, IsCheckParams = opts$check)
	if (u < p) {
		res <- FGfun_mfd(M, U, u)
		prob <- res$prob
		mani.defn <- res$mani.defn
		gamma <- ManifoldOptim::manifold.optim(prob, mani.defn, 
																					 method = opts$method, mani.params = mani.params, 
																					 solver.params = solver.params, x0 = Gamma_init)
		Gamma <- matrix(gamma$xopt, p, u)
	}
	else {
		Gamma <- diag(p)
	}
	Gamma
}


## FG BIC function 
FG_bic <- function(M, U, n){
	which.min(sapply(1:p, function(k) {
		FGfun(manifoldFG(M=M, U=U, u=k), M=M, U=U)$F + k*log(n)/n
	}))
}

## weighted envelope estimation using FG optimization
weighted_envFG <- function(M, U, n, betahat){
	p <- length(betahat)
	Gammas <- sapply(1:p, function(k) {
		manifoldFG(M=M, U=U, u=k)
	})
	bic_val <- n * sapply(1:p, function(k) {
		FGfun(Gammas[[k]], M=M, U=U)$F + k*log(n)/n
	})
	min_bic <- min(bic_val)
	w <- exp(min_bic - bic_val) / sum(exp(min_bic - bic_val))
	rowSums(do.call(cbind, lapply(1:p, function(k){
		w[k] * tcrossprod(Gammas[[k]]) %*% betahat 
	})))
}

## compute the geometric mean of a vector
geoMean <- function(x){
	prod(x)^(1/length(x))
}

## model boot function
model_boot <- function(model, nboot, numCores, intercept = TRUE){
	fam <- model$family$family
	X <- model$x; Y <- model$y
	n <- nrow(X); p <- ncol(X)
	dat <- model$data
	MU_fit <- NULL
	if(intercept){
		if(fam == "binomial"){
			MU_fit <- Logistic_cov(Y=Y, X=X[, -1], 
														 a=coef(model)[1], b=coef(model)[-1])        
		}
		if(fam == "poisson"){
			MU_fit <- Poisson_cov(Y=Y, X=X[, -1], 
														a=coef(model)[1], b=coef(model)[-1])       
		}
	}
	if(!intercept){
		if(fam == "binomial"){
			MU_fit <- Logistic_cov(Y=Y, X=X, a=0, 
														 b=coef(model))      
		}
		if(fam == "poisson"){
			MU_fit <- Poisson_cov(Y=Y, X=X, a=0, 
														b=coef(model))      
		}    
	}
	M <- MU_fit$M; U <- MU_fit$U
	
	# estimated u
	u_fit1D <- oneD_bic(M = M, U = U, n = n, maxdim = p)$u
	u_fitFG <- FG_bic(M = M, U = U, n = n)
	
	MU_fit_star <- a_star <- b_star <- NULL
	mat <- do.call(rbind, mclapply(1:nboot, mc.cores = numCores, function(j){
		
		dat_star <- dat[sample(1:n, replace = TRUE), ]
		m1_star <- update(model, data = dat_star)
		
		# MLE
		betahat_star <- coef(m1_star)
		if(intercept){
			a_star <- betahat_star[1]
			b_star <- betahat_star[-1]
			if(fam == "binomial"){
				MU_fit_star <- Logistic_cov(Y=m1_star$y, 
																		X=m1_star$x[, -1], a=a_star, b=b_star)        
			}
			if(fam == "poisson"){
				MU_fit_star <- Poisson_cov(Y=m1_star$y, 
																	 X=m1_star$x[, -1], a=a_star, b=b_star)                
			}
		}
		if(!intercept){
			a_star <- 0
			b_star <- betahat_star
			if(fam == "binomial"){      
				MU_fit_star <- Logistic_cov(Y=m1_star$y, 
																		X=m1_star$x[, -1], a=a_star, b=b_star)
			}
			if(fam == "poisson"){      
				MU_fit_star <- Poisson_cov(Y=m1_star$y, 
																	 X=m1_star$x[, -1], a=a_star, b=b_star)
			}      
		}
		
		M_star <- MU_fit_star$M; U_star <- MU_fit_star$U
		
		# env fixedu
		Gamma_fixedu_1D_star <- manifold1D(M = M_star, 
																			 U = U_star, u = u_fit1D)
		betaenv_fixedu_1D_star <- as.numeric(
			Gamma_fixedu_1D_star %*% t(Gamma_fixedu_1D_star) %*% 
				b_star)
		betaenv_fixedu_FG_star <- as.numeric(
			tcrossprod(manifoldFG(M = M_star, 
														U = U_star, u = u_fitFG)) %*% b_star)
		
		# env varu
		u1D_star <- oneD_bic(M = M_star, U = U_star, 
												 n = n, maxdim = p)$u
		Gamma_1D_ustar <- manifold1D(M = M_star, 
																 U = U_star, u = u1D_star)
		betaenv_varu_1D_star <- as.numeric(
			Gamma_1D_ustar %*% t(Gamma_1D_ustar) %*% b_star)
		uFG_star <- FG_bic(M = M_star, U = U_star, n = n)
		betaenv_varu_FG_star <- as.numeric(
			tcrossprod(manifoldFG(M = M_star, 
														U = U_star, u = uFG_star)) %*% b_star)
		
		# weighted envelope
		w_betaenv_1D_star <- weighted_env1D(M = M_star, 
																				U = U_star, n = n, betahat = b_star)
		w_betaenv_FG_star <- weighted_envFG(M = M_star, 
																				U = U_star, n = n, betahat = b_star)		
		
		# estimators
		t(c(b_star, betaenv_fixedu_1D_star, 
				betaenv_varu_1D_star, w_betaenv_1D_star, 
				u1D_star, betaenv_fixedu_FG_star, 
				betaenv_varu_FG_star, w_betaenv_FG_star, 
				uFG_star))
	}))
	
	mat
}


## simulation function
sim_function <- function(n, p = p, nboot, numCores, 
												 SigmaX.half, beta, u = u, family){
	
	# setup
	X <- matrix(rnorm(n*p), nrow = n) %*% SigmaX.half 
	Y <- rep(0, n); m1 <- NULL
	data_sim <- NULL 
	if(family == "binomial"){
		Y <- rbinom(n, size = 1, prob = 
									1 / (1 + exp(-as.numeric(X %*% beta))))    
		data_sim <- as.data.frame(cbind(Y, X))
		m1 <- glm(Y ~ -1 + ., family = "binomial", 
							data = data_sim, x =TRUE, y = TRUE)
	}
	if(family == "poisson"){
		Y <- rpois(n, lambda = exp(as.numeric(X %*% beta)))
		data_sim <- as.data.frame(cbind(Y, X))
		m1 <- glm(Y ~ -1 + ., family = "poisson", 
							data = data_sim, x =TRUE, y = TRUE)
	}
	
	# estimates of beta
	betahat <- coef(m1)
	
	# true u 
	M <- U <- NULL
	if(family == "binomial"){
		MUfit <- Logistic_cov(Y=Y, X=X, a=0, b=betahat)
		M <- MUfit$M; U <- MUfit$U    
	}
	if(family == "poisson"){
		MUfit <- Poisson_cov(Y=Y, X=X, a=0, b=betahat)
		M <- MUfit$M; U <- MUfit$U    
	}
	
	betaenv_true_1D <- 
		as.numeric(tcrossprod(manifold1D(M = M, U = U, u = u)) %*% betahat)
	betaenv_true_FG <- 
		as.numeric(tcrossprod(manifoldFG(M = M, U = U, u = u)) %*% betahat)
	
	# estimated u
	u_fit1D <- oneD_bic(M = M, U = U, n = n, maxdim = p)$u
	betaenv_fixedu_1D <- betaenv_varu_1D <- 
		as.numeric(tcrossprod(manifold1D(M = M, 
																		 U = U, u = u_fit1D)) %*% betahat)
	u_fitFG <- FG_bic(M = M, U = U, n = n)
	betaenv_fixedu_FG <- betaenv_varu_FG <- 
		as.numeric(tcrossprod(manifoldFG(M = M, 
																		 U = U, u = u_fitFG)) %*% betahat)  
	
	# weighted envelope estimates
	w_betaenv_1D <- weighted_env1D(M = M, U = U, n = n, 
																 betahat = betahat)
	w_betaenv_FG <- weighted_envFG(M = M, U = U, n = n, 
																 betahat = betahat)
	
	## estimates of beta
	betaest <- cbind(betahat, betaenv_true_1D, betaenv_fixedu_1D, 
									 betaenv_varu_1D, w_betaenv_1D, betaenv_true_FG, 
									 betaenv_fixedu_FG, betaenv_varu_FG, w_betaenv_FG)
	
	## bootstrap
	M_star <- M; U_star <- U
	mat <- do.call(rbind, mclapply(mc.cores = numCores, 
																 1:nboot, function(j){
																 	
																 	# MLE                                   
																 	m1_star <- update(m1, 
																 										data = data_sim[sample(1:n, replace = TRUE), ])
																 	betahat_star <- coef(m1_star)
																 	
																 	# true u                            
																 	if(family == "binomial"){
																 		MUfit <- Logistic_cov(Y=m1_star$y, X=m1_star$x, 
																 													a=0, b=betahat_star)
																 		M_star <- MUfit$M; U_star <- MUfit$U    
																 	}
																 	if(family == "poisson"){
																 		MUfit <- Poisson_cov(Y=m1_star$y, X=m1_star$x, 
																 												 a=0, b=betahat_star)
																 		M_star <- MUfit$M; U_star <- MUfit$U    
																 	}
																 	betaenv_true_1D_star <- as.numeric(
																 		tcrossprod(manifold1D(M = M_star, 
																 													U = U_star, u = u)) %*% betahat_star)
																 	betaenv_true_FG_star <- as.numeric(
																 		tcrossprod(manifoldFG(M = M_star, 
																 													U = U_star, u = u)) %*% betahat_star)
																 	
																 	# env fixedu
																 	betaenv_fixedu_1D_star <- as.numeric(
																 		tcrossprod(manifold1D(M = M_star, 
																 													U = U_star, u = u_fit1D)) %*% betahat_star)
																 	betaenv_fixedu_FG_star <- as.numeric(
																 		tcrossprod(manifoldFG(M = M_star, 
																 													U = U_star, u = u_fitFG)) %*% betahat_star)      
																 	
																 	# env varu
																 	u_star1D <- oneD_bic(M = M_star, 
																 											 U = U_star, n = n, maxdim = p)$u
																 	betaenv_varu_1D_star <- as.numeric(
																 		tcrossprod(manifold1D(M = M_star, 
																 													U = U_star, u = u_star1D)) %*% betahat_star)
																 	u_starFG <- FG_bic(M = M_star, 
																 										 U = U_star, n = n)
																 	betaenv_varu_FG_star <- as.numeric(
																 		tcrossprod(manifoldFG(M = M_star, 
																 													U = U_star, u = u_starFG)) %*% betahat_star)                                   
																 	
																 	# weighted envelope
																 	w_betaenv_1D_star <- weighted_env1D(M = M_star, 
																 																			U = U_star, n = n, betahat = betahat_star)
																 	w_betaenv_FG_star <- weighted_envFG(M = M_star, 
																 																			U = U_star, n = n, betahat = betahat_star)
																 	
																 	# output                             
																 	t(c(betahat_star, betaenv_true_1D_star, betaenv_fixedu_1D_star, 
																 			betaenv_varu_1D_star, w_betaenv_1D_star, u_star1D, 
																 			betaenv_true_FG_star, betaenv_fixedu_FG_star, 
																 			betaenv_varu_FG_star, w_betaenv_FG_star, u_starFG,
																 			m1_star$iter))
																 }))
	
	## output
	out <- list(betaest = betaest, 
							betahat_star = mat[, 1:p], 
							betaenv_true_1D_star = mat[, (p+1):(2*p)], 
							betaenv_fixedu_1D_star = mat[, (2*p+1):(3*p)], 
							betaenv_var_1D_star = mat[, (3*p+1):(4*p)], 
							w_betaenv_1D_star = mat[, (4*p+1):(5*p)], 
							u_star1D = mat[, 5*p+1],
							betaenv_true_FG_star = mat[, (5*p+2):(6*p+1)], 
							betaenv_fixedu_FG_star = mat[, (6*p+2):(7*p+1)], 
							betaenv_var_FG_star = mat[, (7*p+2):(8*p+1)], 
							w_betaenv_FG_star = mat[, (8*p+2):(9*p+1)],
							u_starFG = mat[, 9*p+2])
	out
}


## unpack function for inference
unpack <- function(xlist, beta){
	
	# intermediary quantities
	#betahat_star, betaenv_true_1D_star, betaenv_fixedu_1D_star, 
	#betaenv_varu_1D_star, w_betaenv_1D_star, u_star1D, 
	#betaenv_true_FG_star, betaenv_fixedu_FG_star, 
	#betaenv_varu_FG_star, w_betaenv_FG_star, u_starFG
	names_u <- c("u_star1D", "u_starFG")
	names_star <- c("betahat_star", "betaenv_true_1D_star", 
									"betaenv_fixedu_1D_star", "betaenv_var_1D_star", 
									"w_betaenv_1D_star", "betaenv_true_FG_star", 
									"betaenv_fixedu_FG_star", "betaenv_var_FG_star", 
									"w_betaenv_FG_star")
	betaest <- xlist$betaest
	p <- nrow(betaest)
	names_est <-  colnames(betaest)
	
	# means
	means <- lapply(names_star, function(xx){ 
		colMeans(xlist[[xx]], na.rm = TRUE) 
	})
	names(means) <- paste("means_", names_star, sep = "")
	
	# get variances
	variances <- lapply(names_star, function(xx){ 
		var(xlist[[xx]], na.rm = TRUE) 
	})
	names(variances) <- paste("variances_", names_star, sep = "")
	
	# get multivariate coverage
	n <- nrow(xlist[[6]])
	multi <- lapply(1:ncol(betaest), function(j){
		as.numeric( t(beta -  betaest[, j]) %*% 
									solve(variances[[j]]) %*% (beta -  betaest[, j])
								< p * (n - 1) / (n - p) * qf(0.05, p, n - p, lower = F)          
		)
	})
	names(multi)  <- paste("chisquared_", names_star, sep = "")
	
	# get determinants
	dets <- lapply(variances, function(xx){ det(xx)^(1/p) })
	names(dets) <- paste("pthrootdets_", names_star, sep = "")
	
	# get standard deviations
	sds <- lapply(variances, function(xx){ sqrt(diag(xx)) })
	names(sds) <- paste("sds_", names_star, sep = "")
	
	# null probability
	nulls <- lapply(1:ncol(betaest), function(j){ 
		as.numeric(abs(betaest[, j] / sds[[j]]) < qnorm(0.975))
	})
	names(nulls) <- paste("nulls_", names_star, sep = "")
	
	# coverage probability
	coverages <- lapply(1:ncol(betaest), function(j){ 
		as.numeric(
			(beta > betaest[, j] - qnorm(0.975) * sds[[j]]) &
				(beta < betaest[, j] + qnorm(0.975) * sds[[j]])
		)
	})
	names(coverages) <- paste("coverages_", names_star, sep = "")
	
	# get ratios
	ratios <- lapply(names(sds)[-1], function(xx){ 
		sds[[1]] / sds[[xx]]
	})
	names(ratios) <- paste("ratio_", names_star[-1], sep = "")
	
	# get bias 
	bias <- lapply(names_star, function(xx){ 
		colMeans(xlist[[xx]], na.rm = TRUE) - beta
	})
	names(bias) <- paste("bias_", names_star, sep = "")
	
	# get MSE
	MSE <- lapply(names_star, function(xx){ 
		#print(dim(xlist[[xx]] - beta))
		#print(dim(tcrossprod(xlist[[xx]] - beta)))
		#sum(eigen(tcrossprod(xlist[[xx]] - beta))$val)
		#mean(apply(xlist[[xx]] - beta, 1, function(xx) sum(eigen(tcrossprod(xx))$val)))
		foo <- xlist[[xx]] - beta
		nboot <- nrow(foo)
		bar <- lapply(1:nboot, function(j) tcrossprod(foo[j, ]))
		baz <- bar[[1]] / nboot
		for(j in 2:nboot) baz <- baz + bar[[j]] / nboot
		det(baz)^(1/p)
		
	})
	names(MSE) <- paste("MSE_", names_star, sep = "")
	
	# get distribution of u
	u_boot <- lapply(names_u, function(xx){
		out <- rep(0, p)
		tab <- table(xlist[[xx]]) / nboot
		loc <- as.numeric(names(tab))
		out[loc] <- tab
		out
	})
	names(u_boot) <- c("u_star1D", "u_starFG")
	
	out <- c(ratios, coverages, nulls, dets, multi, bias, MSE, u_boot)
	out
}



