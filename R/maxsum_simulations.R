#' Power analysis for projected score test
#'
#' @description Conducts power analysis for PST along with several other methods
#' for range of user-specified mbetas, ks, and ns. User uses this function to conduct
#' simulation study. Offers option for parallelization. Wraps and combines results from sim_setup() and pstest()
#'
#' @param nsim number of simulations to conduct to assess power, defaults to 500
#' @param seed chosen seed for simulations, defaults to 2019
#' @param mbetas vector of mean beta values
#' @param ks vector of percentage of independent variables with nonzero signal
#' @param n number of observations, defaults to 100
#' @param p number of betas, defaults to 1000
#' @param model can be specified as 'normal' (default) for linear regression, otherwise does logistic regression
#' @param sigma defaults to 1
#' @param alpha significance level, defaults to 0.05
#' @param rho spatial correlation in G parameter, AR1 structure, efaults to 0.9
#' @param betasp indicator of presence of spatial information, defaults to TRUE
#' @param rs investigator-specified set of "contrasts" of G, defaults to c(10, 20, 50)
#' @param mc.cores number of cores to run on, defaults to 1
#' @param doplot if TRUE, makes plots; if FALSE, does not. This produces one power plot per k across a range of mbetas
#' @importFrom graphics legend plot points
#' @return A data frame of power values for PST as well as aSPU, SKAT, and Sum for a range of mbetas and ks. Also plots the power curves.
#'
#' @export

pst_sim = function(nsim = 500,
                   seed = 2019,
                   mbetas = c(0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.06,  0.07, 0.08, 0.09, 0.1),
                   ks = c(40, 60, 70, 80, 100),
                   n = 100,
                   p = 1000,
                   model = 'normal',
                   sigma = 1,
                   alpha = 0.05,
                   rho = 0.9,
                   betasp = TRUE,
                   rs = c(10, 20, 50),
                   mc.cores = 1,
                   doplot = TRUE){

  sobj = sim_setup(n = n, p = p, model = model, sigma = sigma, nsim = nsim,
                   seed = seed, rho = rho, rs = rs)
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


  results = #collecting power results
    lapply(mbetas, function(mbeta) {
        ks.out = lapply(ks, function(kperc) {
          #print(mbeta)
          pstest(nsim = nsim, seed = seed, mbeta = mbeta, kperc = kperc,
                 n = n, p = p,
                 model = model, sigma = sigma, alpha = alpha, betasp = betasp,
                 rs = rs, mc.cores = mc.cores,
                 Gprime = Gprime, GQs = GQs,
                 GQs2 = GQs2, R1nams = R1nams, R2nams = R2nams,
                 nams = nams, H1 = H1,
                 A = A, G = G, linkatlambda = linkatlambda,
                 simresults = simresults,
                 powresults = powresults)
        })

        ks.pow = data.frame(matrix(nrow = length(ks), ncol = 11))
        for(i in 1:nrow(ks.pow)) {ks.pow[i,] = ks.out[[i]]$powresults} #combine all pow results for this mbeta into one dataframe
        names(ks.pow) = names(ks.out[[1]]$powresults) #add names
        return(ks.pow)
    })

  fullpow = do.call("rbind", results) #collapsing power results

  #makes plots for each k
  if(doplot == TRUE){
    #pdffile = paste('n', n, '_p', p, '_model', model, '_sigma', sigma, '_nsim', nsim, '_alpha', alpha, '_rho', rho, '_seed', seed, '.pdf', sep='' )
    shapes = cbind(c(rep(1:2, each=length(rs)), rep(3, ncol(fullpow) - length(rs)*2) ), c(rep(1:length(rs), 2), 1:(ncol(fullpow) - length(rs)*2)) )
    #pdf(pdffile, height=6, width=8)
    for (i in 1:length(ks)){

      k = ks[i]
      subpow = fullpow[ fullpow$k == k, ]
      subpow = subpow[ ! duplicated(subpow$mbeta) ,]
      subpow = subpow[ order(subpow$mbeta),]

      plot(y=subpow$aSPU, x=subpow$mbeta, xlab='Maximum Beta', ylab='Power',
         		type='n', ylim=c(0,1), main=paste('n =',n, 'p =', p, '2k =', k))
     	sapply(names(fullpow)[1:(ncol(fullpow)-2)], function(x) points(subpow$mbeta, subpow[,x],
     		col=shapes[which(names(fullpow)==x),2], pch=shapes[which(names(fullpow)==x),1], type='b') )
       legend('bottomright', legend=names(fullpow)[1:(ncol(fullpow)-2)], col=shapes[1:(ncol(fullpow)-2),2], pch=shapes[1:(ncol(fullpow)-2), 1],
              y.intersp=0.9, bg='white')

    }
  }

  return(fullpow)

}
