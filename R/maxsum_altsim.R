#' Performs power analysis for the projected score test under a single alternative hypothesis
#' simulation under a few alternative hypotheses
#'
#' @description This is an indep function which calls setup from maxsum_altsim_setup.R,
#' which sets up all of the variables needed for one value of kperc, percentage of independent variables with nonzero signal
#' and one value of mbeta, mean coefficient for non zero effects uniform around that value
#' That script loops through the variable k and beta
#'
#' @param n defaults to 100
#' @param p defaults to 1000
#' @param model can be specified as 'normal' (default) for linear regression, otherwise does logistic regression
#' @param sigma defaults to 1
#' @param nsim defaults to 500
#' @param alpha significance level, defaults to 0.05
#' @param seed set a seed for the power calculation, defaults to 2019
#' @param mbeta mean coefficient for nonzero effects, defaults to 0
#' @param kperc percentage of independent variables with nonzero signal, defaults to 40
#' @param rho spatial correlation in G parameter, AR1 structure, defaults to 0.9
#' @param betasp indicator of presence of spatial information, defaults to TRUE
#' @param rs investigator-specified set of "contrasts" of G, defaults to c(10, 20, 50)
#' @param mc.cores number of cores to run simulation on
#' @import parallel
#' @import foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom stats ecdf glm lm model.matrix optim pchisq predict quantile rbinom resid rnorm
#' @importFrom utils capture.output
#' @importFrom mvtnorm qmvnorm
#' @importFrom iterators icount
#' @return A list of the entire matrix of simulation results, the power results, and "out"
#' @keywords projected score test
#' @export

pstest = function(n = 100, p = 1000, model = 'normal',
                  sigma = 1, nsim = 500, alpha = 0.05, seed = 2019,
                  mbeta = 0, #mean bet
                  kperc = 40, #current k
                  rho = 0.9,
                  betasp = TRUE,
                  rs = c(10, 20, 50),
                  mc.cores = 1){

  #set up simulation
  sobj = sim_setup(n = n, p = p, model = model, sigma = sigma, nsim = nsim,
                   alpha = alpha, seed = seed,
                   rho = rho,
                   betasp = betasp,
                   rs = rs) #create setup object with params specified in main function

  Gprime = sobj$Gprime
  GQs = sobj$GQs
  GQs2 = sobj$GQs2
  R1nams = sobj$R1nams
  R2nams = sobj$R2nams
  nams = sobj$nams
  simresults = sobj$simresults
  powresults = sobj$powresults
  H1 = sobj$H1
  A = sobj$A
  G = sobj$G
  linkatlambda = sobj$linkatlambda

  betas = seq(0, mbeta, length.out=kperc/2+1)[-1]

  # spatial information
  if(betasp){
  	beta = c(rep(0,kperc), -1*betas, rep(0, kperc*5), betas,  rep(0, p-kperc*7) )
  # no spatial information (besides correlation amongst indep variables)
  } else {
  	beta = rep(0, p)
  	beta[round(seq(10, 990, length.out=kperc))] = c(betas, betas) * c(-1, -1, 1, 1, -1,  1)
  }


  #parallelize the simulations
  cl = makeCluster(mc.cores)
  registerDoParallel()
  r = foreach(icount(nsim), .combine = rbind) %dopar% {
      if(tolower(model)=='normal'){
      # Simulate Y under no covariates
      Y = Gprime %*% beta + rnorm(n, 0, sigma)
      # regress out mean
      null.mod = lm(Y ~ 1)
      X = model.matrix(null.mod)
      H = X %*% solve(t(X) %*% X) %*% t(X)
      Y0 = resid(null.mod)
      # covariance structure for Y after removal of covariates.
      # Based on fisher information.
      sigmasqhat = sum(Y0^2)/(n-ncol(X))
      Vtheta0 = sigmasqhat*(diag(n)-H)

    # else assume logistic regression
    } else {
      Y = rbinom(n, size=1, prob=1 - (1/(1 + exp(-Gprime %*% beta)) ) )
      null.mod = glm(as.factor(Y) ~ 1, family='binomial')
      X = model.matrix(null.mod)
      Yhat = predict(null.mod, type='response')

      Y0 = Y - Yhat
      Vtheta0 = diag(Yhat*(1-Yhat))
      XtVtheta0 = t(X) %*% Vtheta0
      Vtheta0 = Vtheta0 - t(XtVtheta0) %*% solve( XtVtheta0 %*% X) %*% XtVtheta0

      A = svd(Vtheta0, nu=n, nv=0)
      linkatlambda = svd( diag(sqrt(A$d)) %*% t(A$u) %*% Gprime, nu=0, nv=0)$d^2
    }

    ### Estimation of R_1
    R1s = sapply(GQs, function(gq) MaxSumEstimate(Y0, Vtheta0, gq)$Rmaxsum )
    R1s2 = sapply(GQs2, function(gq) MaxSumEstimate(Y0, Vtheta0, gq)$Rmaxsum )

    ### Estimation of SKAT
    # same for continuous and categorical
    #SKAT =  t(Y0) %*% Gprime %*% t(Gprime) %*% Y0

    # adaptive SPU (Sum of Powered score (U?))
    # I think you just have to scale Y maybe.
    if(tolower(model) == 'normal'){
      out = aSPU(Y, Gprime, cov=X, resample='perm', model='gaussian')
    } else {
      out = aSPU(Y, Gprime, cov=X, resample='perm', model='binomial')
    }



    if(tolower(model) == 'normal'){
    # For normal data
    # G is the rotated design after removing effects of covariates
    pvalueR1 = sapply(1:length(rs), function(r) pratioqform(nrow(G), rs[r], R1s[r]) )
    pvalueR12 = sapply(1:length(rs), function(r) pratioqform(nrow(G), rs[r], R1s2[r]) )
    } else {
      # asymptotic
      pvalueR1 = 1-pchisq(R1s, df=rs)
      pvalueR12 = 1-pchisq(R1s2, df=rs)
    }

    # unknown variance - needs to be scaled to smaller numbers for imhof to work
    # hence dividing by sum(linkatlambda)

    # value of adaptive test is not meaningful I think
    simresults = data.frame(NA)
    simresults[, R1nams] = c(R1s, pvalueR1)
    simresults[, R2nams] = c(R1s2, pvalueR12)
    simresults[,c('aSPU', 'aSPU_pvalue')] = c(out$Ts['aSPU'], out$pvs['aSPU'])
    simresults[,c('SKAT', 'SKAT_pvalue')] = c(out$Ts['SPU2'], out$pvs['SPU2'])
    simresults[,c('Sum', 'Sum_pvalue')] = c(out$Ts['SPU1'], out$pvs['SPU1'])

    simresults[,-1]

    #if(! i %% 200) cat('done', i, '\n')
  }
  registerDoSEQ()

  simresults = as.data.frame(r)
  powresults[1,] = c(colMeans(simresults[,grep('_pvalue$', nams)]<=alpha, na.rm = T), kperc, mbeta)
  powresults = as.data.frame(powresults)


  return(list("simresults" = simresults, "powresults" = powresults))

}



