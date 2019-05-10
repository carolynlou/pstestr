#' Sets up necessary parameters for power calculation simulation
#' @description Support function for pstest(), which runs the power calculations that are wrapped in pst_sim(). Sets
#' up objects needed to do power calculation.
#'
#' @param nsim defaults to 500
#' @param seed set a seed for the generation of matrices needed for the various tests, defaults to 2019
#' @param n defaults to 100
#' @param p defaults to 1000
#' @param model can be specified as 'normal' (default) for linear regression, otherwise does logistic regression
#' @param sigma defaults to 1
#' @param rho spatial correlation in G parameter, AR1 structure, defaults to 0.9
#' @param rs investigator-specified set of "contrasts" of G, defaults to c(10, 20, 50)
#' @return A list including spatial correlation parameters, empty dataframes for simresults
#'  and powresults, parameters for distribution of SKAT statistic
#'
#' @export


sim_setup = function(nsim = 500, seed = 2019,
                     n = 100, p = 1000,
                     model = 'normal',
                     sigma = 1,
                     rho = 0.9,
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
  contrasts = lapply(rs, function(r) svd( splines::bs(1:p, df=r, intercept=TRUE ), nu=r, nv=0)$u )
  #contrasts = lapply(rs, buildbasis, p=p)
  O = svd(Gprime, nu=0, nv=n)$v
  contrasts2 = lapply(rs, function(r) O[,1:r]  )

  # These are nxr matrices used for both models
  GQs = lapply(contrasts, function(con) Gprime %*% con)
  GQs2 = lapply(contrasts2, function(con) Gprime %*% con)
  # don't need these anymore
  rm(contrasts, contrasts2)

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

}
