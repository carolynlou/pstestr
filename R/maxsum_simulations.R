#' Conducts power analysis for PST along with several other methods
#' for range of user-specified mbetas, ks, and ns. User uses this function to conduct
#' simulation study. Offers option for parallelization
#'
#' @param nsim Number of simulations to conduct to assess power. Defaults to 500
#' @param mbetas Vector of mean beta values
#' @param ks Vector of percentage of independent variables with nonzero signal
#' @param n Defaults to 100
#' @param p Defaults to 1000
#' @param model Can be specified as 'normal' (default) for linear regression, otherwise does logistic regression
#' @param sigma Defaults to 1
#' @param alpha significance level, Defaults to 0.05
#' @param seed set a seed for the power calculation. defaults to 2019
#' @param rho spatial correlation in G parameter, AR1 structure. defaults to 0.9
#' @param betasp indicator of presence of spatial information. defaults to TRUE
#' @param rs Investigator specified set of "contrasts" of G. defaults to c(10, 20, 50)
#' @param mc.cores Number of cores to run on. Defaults to 1
#'
#' @export

pst_sim = function(nsim = 500,
                   mbetas = c(0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.06,  0.07, 0.08, 0.09, 0.1),
                   ks = c(40, 60, 70, 80, 100),
                   n = 100,
                   p = 1000,
                   model = 'normal',
                   sigma = 1,
                   alpha = 0.05,
                   seed = 2019,
                   rho = 0.9,
                   betasp = TRUE,
                   rs = c(10, 20, 50),
                   mc.cores = 1){

  results = #collecting power results
    lapply(mbetas, function(mbeta) {
        ks.out = lapply(ks, function(curk) {
          pstest(nsim = nsim,
                 mbeta = mbeta,
                 kperc = curk,
                 n = n,
                 p = p,
                 model = model,
                 sigma = sigma,
                 alpha = alpha,
                 seed = seed,
                 rho = rho,
                 betasp = betasp,
                 rs = rs,
                 mc.cores = mc.cores)
        })

        ks.pow = data.frame(matrix(nrow = length(ks), ncol = 11))
        for(i in 1:nrow(ks.pow)) {ks.pow[i,] = ks.out[[i]]$powresults} #combine all pow results for this mbeta into one dataframe
        names(ks.pow) = names(ks.out[[1]]$powresults) #add names
        return(ks.pow)
    })

  fullpow = do.call("rbind", results) #collapsing power results

  #makes plots for each k
  shapes = cbind(c(rep(1:2, each=length(rs)), rep(3, ncol(fullpow) - length(rs)*2) ), c(rep(1:length(rs), 2), 1:(ncol(fullpow) - length(rs)*2)) )
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

  return(results)

}
