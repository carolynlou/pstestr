#' Conducts power analysis for PST along with several other methods
#' for range of user-specified mbetas, ks, and ns. User uses this function to conduct
#' simulation study. Offers option for parallelization
#'
#' @param nsim Number of simulations to conduct to assess power. Defaults to 500
#' @param mbetas Vector of mean beta values
#' @param ks Vector of percentage of independent variables with nonzero signal
#' @param mc.cores Number of cores to run on. Defaults to 1

pst_sim = function(nsim = 500,
                   mbetas = c(0, 0.05, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.06,  0.07, 0.08, 0.09, 0.1),
                   ks = c(40, 60, 70, 80, 100), mc.cores = 1){

  # results =
    start = proc.time()
    parallel::mclapply(mbetas, function(mbeta) {
      lapply(ks, function(curk) {pstest(nsim = nsim, mbeta = mbeta, kperc = curk)})},
      mc.cores = mc.cores)
  proc.time() - start

  # for(mbeta in mbetas){
  #   for(curk in ks){
  #
  #     pstest(nsim = nsim, mbeta = mbeta, kperc = curk)
  #
  #   }
  # }


}
