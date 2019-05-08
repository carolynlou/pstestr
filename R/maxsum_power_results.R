# maxkat_power_results. Just combines the power results
# into one file. Will probably also create some plots.
# Assumes you ran maxkat_altsim_setup.R
# in the same directory because it uses setupfile variable
# and other variables defined from running that script.
#setwd('~/Documents/work/maxkat')
#source('maxkat_functions.R')
#
# powdir = './power_results/'
# setupfile = './power_files/setup_n100_p1000_modelnormal_sigma1_nsim1000_alpha0.05_rho0.9_seed2016.rdata'
# load(setupfile)
# #ks= c(10, 20, 40, 80)
#
# # get first power results file to get format of power output table
# powout = load(sub('setup_', '', sub('.rdata', paste('_k', ks[1], '_mbeta', mbetas[1], '_betasp', betasp, '.rdata', sep=''), sub('power_files', 'power_results', setupfile) ) ) )
# fullpow = powresults[-1,]
#
# # combines them all into one data frame
# for( curk in ks){
#   for(mbeta in mbetas){
#     powout = sub('setup_', '', sub('.rdata', paste('_k', curk, '_mbeta', mbeta, '_betasp', betasp, '.rdata', sep=''), sub('power_files', 'power_results', setupfile) ) )
#     load(powout)
#     fullpow = rbind(fullpow, powresults)
#   }
# }
#
# pdffile = sub('rdata', 'pdf', sub('setup_', '', setupfile))
# #shapes = expand.grid(1:6, 1:10)
# shapes = cbind(c(rep(1:2, each=length(rs)), rep(3, ncol(fullpow) - length(rs)*2) ), c(rep(1:length(rs), 2), 1:(ncol(fullpow) - length(rs)*2)) )
# pdf(pdffile, height=6, width=8)
# for( k in ks){
# 	subpow = fullpow[ fullpow$k==k,]
# 	subpow = subpow[ ! duplicated(subpow$mbeta) ,]
# 	subpow = subpow[ order(subpow$mbeta),]
# 	plot(y=subpow$aSPU, x=subpow$mbeta, xlab='Maximum Beta', ylab='Power',
# 		type='n', ylim=c(0,1), main=paste('n =',n, 'p =', p, '2k =', k))
# 	sapply(names(fullpow)[1:(ncol(fullpow)-2)], function(x) points(subpow$mbeta, subpow[,x],
# 		col=shapes[which(names(fullpow)==x),2], pch=shapes[which(names(fullpow)==x),1], type='b') )
#   legend('bottomright', legend=names(fullpow)[1:(ncol(fullpow)-2)], col=shapes[1:(ncol(fullpow)-2),2], pch=shapes[1:(ncol(fullpow)-2), 1],
#          y.intersp=0.9, bg='white')
# }
# dev.off()
