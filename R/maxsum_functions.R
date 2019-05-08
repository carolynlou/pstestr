#' Collection of helper functions for main simulation
#'
#' @import CompQuadForm
#' @import splines
#' @import aSPU
#' @import geigen
#' @import AssotesteR
#' @import abind
#' @import quantreg
#' @import doParallel
#'


# FUNCTIONS FOR MAXKAT TEST
#library(pracma) # for repmat

# library(CompQuadForm) # for LC of chi-squareds
# library(splines)
# library(aSPU)
# library(geigen)
# library(AssotesteR)
# library(abind)
# library(quantreg)

#library(rARPACK)
#library(microbenchmark)
#library(SKAT)


# Rayleigh Quotient
# function to estimate alpha vector
# This optimizes a Rayleigh quotient which is the solution
# to a generalized eigenvalue problem. alpha is the first
# eigenvector.

rayleigh.opt <- function(u=NULL, v=NULL, par=NULL) {
  (t(par) %*% u %*% par)/(t(par) %*% v %*% par)
}

# Sum of Projection Matrix
# sums the p^2 nXn sub-blocks of an
# npXnp projection matrix
sumProj <- function(P, n, p){
  out = matrix(0, n, n)
  for(j in 1:p){
    for(k in 1:p){
      # ijth block of P
      out = out + P[ (n*(j-1)+1):(n*j), (n*(k-1)+1):(n*k) ]
    }
  }
  out
}


# Sum of Orthogonal Basis Vectors
# sums the p^2 nXn sub-blocks of an
# npXnp projection matrix by summing the
# outer product of the underlying orthogonal
# basis vectors. This is the slowest step.
# sums the p^2 nXn sub-blocks of an
# npXnp projection matrix by summing the
# outer product of the underlying orthogonal
# basis vectors.
sumOrth <- function(Psi, n, p){
  out = matrix(0, n, n)
  for(j in 1:p){
    for(k in 1:p){
      # ijth block of P
      out = out + Psi[ (n*(j-1)+1):(n*j), ] %*% t(Psi[(n*(k-1)+1):(n*k),] )
    }
  }
  out
}


# Sum of Orthogonal Basis Vectors, faster
# sums the p^2 nXn sub-blocks of an
# npXnp projection matrix by summing the
# outer product of the underlying orthogonal
# basis vectors. This is the slowest step.
# sums the p^2 nXn sub-blocks of an
# npXnp projection matrix by summing the
# outer product of the underlying orthogonal
# basis vectors. This is roughly twice as fast as sumOrth().

sumOrth2 <- function(Psi, n, p){
  out1 <- out <- matrix(0, n, n)
  for(j in 1:p){
    out1 = out1 + Psi[ (n*(j-1)+1):(n*j), ] %*% t(Psi[(n*(j-1)+1):(n*j),] )
    k=j+1
    while(k <= p){
      # ijth block of P
      out = out + Psi[ (n*(j-1)+1):(n*j), ] %*% t(Psi[(n*(k-1)+1):(n*k),] )
      k = k+1
    }
  }
  out1 + out + t(out)
}
# test = microbenchmark(sumOrth2(Psi, n, p), sumOrth(Psi, n, p), times=50L )


# Sum of Orthogonal Basis Vectors, different
# sums the p^2 nXn sub-blocks of an
# npXnp projection matrix by summing the
# outer product of the underlying orthogonal
# basis vectors. This is the slowest step.
# sums the p^2 nXn sub-blocks of an
# npXnp projection matrix by summing the
# outer product of the underlying orthogonal
# basis vectors. This functions slightly differently.

sumOrth3 <- function(Psi, i, j, n, p){
  out <- matrix(0, n, n)
  for(ell in 1:p){
    out = out + Psi[ (n*(ell -1)+1):(n* ell), i, drop=FALSE] %*% t(Psi[(n*(ell-1)+1):(n*ell), j, drop=FALSE] )
  }
  out + t(out)
}


# Estimates the numerator of Q0
Q0num = function(Y, contrasts, G){
  cons = apply(contrasts, 2, function(con) t(Y) %*% G %*% diag(con) )
  t(cons) %*% cons
}


# Estimates the denominator of Q0
Q0denom = function(contrast, dGtG){
  cons = apply(contrast, 2, function(con) dGtG %*% diag(con) )
  Vmat = t(cons) %*% cons
}

# Estimates the numerator of R0
R0num = function(Y, contrasts, G){
  cons = t(Y) %*% G %*% contrasts
  t(cons) %*% cons
}

# Estimates the denominator of Q0
R1denom = function(contrast, G){
  cons = G %*% contrast
  Vmat = t(cons) %*% cons
}


# perform 1 permutation of the data
Q0_nullperm = function(Y, contrasts, G, Vmat){
  yperm = sample(Y)
  Umat = Q0num(yperm, contrasts, G)
  optim(par=rep(1, ncol(contrasts)), fn=rayleigh.opt, u=Umat, v=Vmat, control = list(fnscale=-1))$value
}

# compute the pvalue for the ratio of quadratic forms
pratioqform = function(nminusm, r, R){
  ds = c(rep(nminusm, r), rep(0, nminusm-r)) - R
  imhof(q=0, lambda=ds/sum(abs(ds) ))$Qq
}

MaxSumEstimate = function(Y, covY, GQ){
  Vinv = solve(t(GQ) %*% covY %*% GQ)
  coefs = c( Vinv %*% t(GQ) %*% Y)
  Rmaxsum = c(t(Y) %*% GQ %*% coefs)
  list(Rmaxsum = Rmaxsum, bcoefs=coefs)
}

buildbasis = function(bdim, p){
	inds = 1:p
	vec = seq(0, p, length.out=bdim+1)
	svd(sapply(1:bdim, function(x) ifelse(inds>vec[x] & inds <= vec[x+1], 1, 0) ), nv=0)$u
}

# get a single value which gives the 1-alpha treshold for
# the unnormed weight vector.
normalweightthreshold = function(G, W, alpha=0.05 ){
	V = svd(G %*% W, nu=0)
	d = V$d
	V = V$v
	WV = W %*% V
	covar = rowSums(t(WV) )^2 * 1/d^2
	qmvnorm(1-alpha, sigma=diag(covar), tail='both.tails')$quantile
}


### inference on the estimated weights
#Y=Y, GQ=GQs2[[x]], covY=covY, Q=contrasts2[[x]], nsimint=nsimint)$pvals
maxsumInf = function(Y, GQ, covY, Q, nsimint=2000, alpha=0.05){
	r = ncol(Q)
	# Things for the variance
	#GtQY = t(GQ) %*% Y
	#Vinv = solve( t(GQ) %*% covY %*% GQ - GtQY %*% t(GtQY) )
	Vinv = solve( t(GQ) %*% covY %*% GQ )
	Vsqrtinv = svd(Vinv, nv=0)
	Vsqrtinv = Vsqrtinv$u %*% diag(sqrt(Vsqrtinv$d), r) %*% t(Vsqrtinv$u)
	QVsqrt = Q %*% solve(Vsqrtinv)
	# sqrt of diagonal of Q V Q^T (covariance of the projected U vector)
	diagcov = sqrt(rowSums(QVsqrt^2))
	diagcov = ifelse(diagcov==0, 1, diagcov)


	# Projected scores
	PU = Q %*% t(GQ) %*% Y
	# projected scores are standardized here.
	zetastar = PU/diagcov

	pinftynorm = function(z){
	  max(abs(QVsqrt %*% z)/diagcov) # divide by variances
	}

	# Monte carlo integration -- or just simulation under the null
	cat('performing simulations\n')
	simnormals = matrix(rnorm(nsimint*r), ncol=r, nrow=nsimint)
	siminftynorms = apply(simnormals, 1, pinftynorm)
	Fn = ecdf(siminftynorms)
	# define a quantile function and compute voxelwise pvalues (and return the result)
	list(PUstd = zetastar, rejreg = quantile(siminftynorms, probs=1-alpha), pvals = 1-Fn(abs(zetastar)), Fn=Fn )
}

buildbasis = function(r, labelvec){
	basis = cut(labelvec, c(-1, seq(0, max(labelvec), length.out=r+1)) )
	# first column is no label
	basis = model.matrix(~ basis)[,-1]
	# don't need to square for norm b/c they're indicators
	basis = apply(basis, 2, function(x) x/sqrt(sum(x)) )
}

# This function residualizes to covariates and performs SVD of the imaging data
svdfunc = function(Gtilde=NULL, form=NULL, pcinds=NULL, data=NULL){
	# design matrix for covariates
	X = model.matrix(form, data=data)
	# dof assumes that X is full rank
	A = svd((diag(nrow(X)) - X %*% solve(t(X) %*% X) %*% t(X) ), nv=0, nu=nrow(X) - ncol(X) )$u
	G = t(A) %*% Gtilde
	PCs = svd(G, nu=0, nv=max(pcinds) )
	ds = PCs$d
	PCs = PCs$v
	#ds = svd( G, nv=0, nu=0 )$d
	return(list(ds=ds, PCs=PCs))
}

# computes covariance of Y for bernoulli data
# This is not really covY yet
covYfunc = function(Y0, X){
	null.mod = glm(as.factor(Y0) ~ X, family='binomial')
	Yhat = predict(null.mod, type='response')

	covY = diag(Yhat*(1-Yhat))
	XtcovY = t(X) %*% covY
	# now it is
	covY = covY - t(XtcovY) %*% solve( XtcovY %*% X) %*% XtcovY
	list(covY=covY, Y = (Y0-Yhat))
}

# Code by JP Fortin for the CAT plots
getOverlap <- function(ranks1, ranks2){
	kk <- 1:length(ranks1)
	r <- factor(pmax(ranks1, ranks2), levels=kk)
	tab <- table(r)
	cs <- cumsum(tab)
	cs/kk
}


# Computes AUC of CAT plot
catauc = function(r, Y1, Y2, G1tilde, G2tilde, covY1, covY2, PCs=NULL){
   # change to PCs1 PCs2 for treating X as random
   GQ1 = G1tilde %*% PCs[,1:r, drop=FALSE]
   GQ2 = G2tilde %*% PCs[,1:r, drop=FALSE]
   # get standardized statistics for CAT plot
   # change to PCs1 PCs2 for treating X as random
   invisible(capture.output(PUstd1 <- maxsumInf(Y=Y1, GQ=GQ1, covY=covY1, Q=PCs[,1:r, drop=FALSE], nsimint=1)$PUstd ))
   invisible(capture.output(PUstd2 <- maxsumInf(Y=Y2, GQ=GQ2, covY=covY2, Q=PCs[,1:r, drop=FALSE], nsimint=1)$PUstd ))

   ### CAT PLOTS ###
   # Code from JP Fortin
   ranks1 <- rank(-abs(PUstd1))
   ranks2 <- rank(-abs(PUstd2))
   o <- getOverlap(ranks1, ranks2)
   mean(o)
}

