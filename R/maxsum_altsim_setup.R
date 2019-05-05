# maxkat_nullsim
# simulation under a few alternative hypotheses
rm(list=ls())
setwd('/Users/louc/Box/Penn/Penn_Y2/BSTA670/simproject/clou_simonvandekar-pst/power_analyses/')
#setwd('~/Documents/work/maxsum/power_analyses')
logdir = './logdir'
dir.create(logdir, showWarnings=FALSE)
# removes old log files
unlink(file.path(logdir, '*'))

#### THINGS THAT ARE FIXED ####
# fixed design
# might have to center
# consider using spatially correlated variables
setup = function(n = 100, p = 1000, model = 'normal', sigma = 1, nsim = 200, alpha = 0.05, seed = 2019,
                 mbetas = c(0, 0.05, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.06,  0.07, 0.08, 0.09, 0.1),
                 ks = c(40, 60, 70, 80, 100),
                 rho = 0.9,
                 betasp = 1,
                 rs = c(10, 20, 50)){

  set.seed(seed)
  Gprime = matrix(rnorm(n*p), n, p)

  #### SPATIAL CORRELATION IN G ####
  # AR1 structure
  if(rho != 0){
  	V = round(rho^abs(outer(1:p, 1:p, "-")), 10)
  	svd.V = svd(V, nu=p, nv=0)
  	Gprime = Gprime %*% diag(sqrt(svd.V$d)) %*% t(svd.V$u)
  	rm(svd.V, V)
  }

  # Investigator specified set of "contrasts" of G.
  # Use a spline basis to capture "spatially" flexible effects
  # INSTEAD OF SPLINE BASIS USE SOMETHING BASED ON A SPATIAL CORRELATION STRUCTURE
  contrasts = lapply(rs, function(r) svd(bs(1:p, df=r, intercept=TRUE ), nu=r, nv=0)$u )
  #contrasts = lapply(rs, buildbasis, p=p)
  O = svd(Gprime, nu=0, nv=n)$v
  contrasts2 = lapply(rs, function(r) O[,1:r]  )

  # These are nxr matrices used for both models
  GQs = lapply(contrasts, function(con) Gprime %*% con)
  GQs2 = lapply(contrasts2, function(con) Gprime %*% con)
  # don't need these anymore
  rm(contrasts, contrasts2)

  # save out objects we will use
  dir.create('power_files', showWarnings = FALSE)
  setupfile = paste('power_files/setup_n', n, '_p', p, '_model', model, '_sigma', sigma, '_nsim', nsim, '_alpha', alpha, '_rho', rho, '_seed', seed, '.rdata', sep='' )

    #### empty simulation results matrix ####
    R1nams = paste(paste('Rspline', rs, sep='_'), rep(c('', '_pvalue'), each=length(rs)), sep='' )
    R2nams = paste(paste('Rcov', rs, sep='_'), rep(c('', '_pvalue'), each=length(rs)), sep='' )
    nams = c(R1nams, R2nams, 'aSPU', 'aSPU_pvalue', 'SKAT', 'SKAT_pvalue', 'Sum', 'Sum_pvalue')
    # 'SKAT', 'SKAT_pvalue', # same as SPU2
    simresults = matrix(NA, nrow=nsim, ncol=length(nams), dimnames=list(NULL, nams))


    #### Empty power results ####
    powresults = matrix(NA, nrow=1, ncol=(2+sum(grepl('pvalue', nams) )), dimnames=list(NULL, c( sub('_pvalue', '', grep('_pvalue$', nams, value=TRUE) ), 'k', 'mbeta') ) )


    #### parameters for distribution of SKAT statistic ####
    # should make this hat matrix more general to work with other covariates
    H1 = matrix(1, n, 1) %*% matrix(1, 1, n)/n
    # sqrt of variance of Y (idempotent)
    A = svd(diag(rep(1, n)) - H1, nv=0)
    # reducing dimensions of G by rank of H1
    A = A$u[ ,round(A$d, 10)>0]
    G = t(A) %*% Gprime
    linkatlambda = svd(G, nu=0, nv=0)$d^2
    # drop small eigenvalues
    linkatlambda = linkatlambda[ round(linkatlambda, 10) >0]

  return(list(Gprime = Gprime, GQs = GQs, GQs2 = GQs2, R1nams = R1nams,
              R2nams = R2nams, nams = nams, simresults = simresults,
              powresults = powresults, H1 = H1, A = A, G = G, linkatlambda = linkatlambda))


  # #### loop through k ####
  # for(curk in ks){
  #   for(beta in mbetas){
  #
  #      system( paste( 'R --file=./maxsum_altsim.R --args', getwd(), setupfile, curk, mbeta) )
  #
  #   }
  # }
}
