#' Performs power analysis for the projected score test under a single alternative hypothesis
#'


# simulation under a few alternative hypotheses
# This is an indep function which calls setup from maxsum_altsim_setup.R,
# which sets up all of the variables needed for one value of k, percentage of independent variables with nonzero signal
# and one value of mbeta, mean coefficient for non zero effects uniform around that value
# That script loops through the variable k and beta


pstest = function(n = 100, p = 1000, model = 'normal', sigma = 1, nsim = 200, alpha = 0.05, seed = 2019,
                  mbetas = c(0, 0.05, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.06,  0.07, 0.08, 0.09, 0.1),
                  ks = c(40, 60, 70, 80, 100),
                  rho = 0.9,
                  betasp = 1,
                  rs = c(10, 20, 50)){

  sobj = setup() #create setup file with params specified in main function

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



  mbeta = mbetas[]
  dir.create('power_results', showWarnings = FALSE)
  powout = sub('setup_', '', sub('.rdata', paste('_k', curk, '_mbeta', mbeta, '_betasp', betasp, '.rdata', sep=''), sub('power_files', 'power_results', setupfile) ) )

  betas = seq(0, mbeta, length.out=curk/2+1)[-1]

  # spatial information
  if(betasp){
  	beta = c(rep(0,curk), -1*betas, rep(0, curk*5), betas,  rep(0, p-curk*7) )
  # no spatial information (besides correlation amongst indep variables)
  } else {
  	beta = rep(0, p)
  	beta[round(seq(10, 990, length.out=curk))] = c(betas, betas) * c(-1, -1, 1, 1, -1,  1)
  }


  for(i in 1:nsim){
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
    # resample='sim' gives incorrect results
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

      # also perform categorical tests
      #nperms = 100
      #assu.out = ASSU.Ord(Y, Gprime, perm=nperms)
      #ssu.out = SSU(Y, Gprime, perm=nperms)
      #asum.out = ASUM.Ord(Y, Gprime, perm=nperms)
      #sum.out = SUM(Y, Gprime, perm=nperms)
      #simresults[, c('aSSU', 'aSSU_pvalue')] = c(assu.out$assu.stat, assu.out$perm.pval)
      #simresults[, c('SSU', 'SSU_pvalue')] = c(ssu.out$ssu.stat, ssu.out$perm.pval)
      #simresults[, c('aSum', 'aSum_pvalue')] = c(asum.out$asum.stat, asum.out$perm.pval)
      #simresults[i, c('Sum', 'Sum_pvalue')] = c(sum.out$asum.stat, sum.out$perm.pval)
    }

    # unknown variance - needs to be scaled to smaller numbers for imhof to work
    # hence dividing by sum(linkatlambda)
    #pvalueSKAT = imhof(q=SKAT/sum(linkatlambda), lambda=linkatlambda/sum(linkatlambda))$Qq



    # output results
    simresults[i, R1nams] = c(R1s, pvalueR1)
    simresults[i, R2nams] = c(R1s2, pvalueR12)
    #simresults[i,c('SKAT', 'SKAT_pvalue')] = c(SKAT, pvalueSKAT)
    # value of adaptive test is not meaningful I think
    simresults[i,c('aSPU', 'aSPU_pvalue')] = c(out$Ts['aSPU'], out$pvs['aSPU'])
    simresults[i,c('SKAT', 'SKAT_pvalue')] = c(out$Ts['SPU2'], out$pvs['SPU2'])
    simresults[i,c('Sum', 'Sum_pvalue')] = c(out$Ts['SPU1'], out$pvs['SPU1'])

    if(! i %% 200) cat('done', i, '\n')
  }

  simresults = as.data.frame(simresults)
  powresults[1,] = c(colMeans(simresults[,grep('_pvalue$', nams)]<=alpha, na.rm = T), curk, mbeta)
  powresults = as.data.frame(powresults)
  return(list=c('simresults', 'powresults', 'out'))

}



