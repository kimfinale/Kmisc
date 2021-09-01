#
#
# # plot(res$time, res$I, type="s")
# # plot( sim_res[[1]]$time, sim_res[[1]]$I, col="#A9A9A9", type="l", ylim=c(0,20), xlim=c(0,300))
# #
# # for( i in 2:nsim ){
# #   lines( sim_res[[i]]$time, sim_res[[i]]$I, col="#A9A9A9" )
# # }
#
# plot( sim_res[[1]]$time, sim_res[[1]]$I, col=rgb(0.5,0.5,0.5,alpha=0.2), type="l", ylim=c(0,200), xlim=c(0,200),
#       xlab="day since introduction", ylab="cumulative number of infecteds (confirmed=red)", main="R0=0.6, I0=10, Pop=1e4")
#
# lines( sim_res[[1]]$time, sim_res[[1]]$H, col=rgb(1,0,0,alpha=0.2) )
#
# for( i in 2:nsim ){
#   lines( sim_res[[i]]$time, sim_res[[i]]$I, col=rgb(0.5,0.5,0.5,alpha=0.2) )
#   lines( sim_res[[i]]$time, sim_res[[i]]$H, col=rgb(1,0,0,alpha=0.2) )
# }
#
#
# # d <- matrix( NA, nrow=500, ncol=200 )
# # d <- unlist( sim_res )
#
# df <- data.frame( t=sim_res[[1]]$time, I=sim_res[[1]]$I, H=sim_res[[1]]$H )
# for( i in 2:length(sim_res) ){
#   df <- rbind( df, data.frame(t=sim_res[[i]]$time, I=sim_res[[i]]$I, H=sim_res[[i]]$H) )
# }
#
# ids <- which( df$t==0 ) # starting point of each simulation
# ext_times <- rep( NA, length(ids) )
#
# for( i in 2:(length(ids)+1) ){
#   nrow <- 1
#   if( i > length(ids) ){
#     nrow <- length(df$t) - ids[i-1] + 1
#   }
#   else{
#     # nrow <- ids[i] - ids[i-1] - 1
#     nrow <- ids[i] - ids[i-1]
#   }
#   for( j in 1:nrow ){
#     if( df$I[ ids[i-1]+ j-1 ] == 0 ){
#       ext_times[ i-1 ] <- df$t[ ids[i-1]+j-1 ]
#       break
#     }
#   }
# }
# ext_times <- ext_times[ !is.na(ext_times) ] # remove NAs
# time_cutoffs <- seq(0, 199, by=10)
# xx <- sapply( 1:length(time_cutoffs), function(x) sum(ext_times < time_cutoffs[x]) )
# probs <- xx/nsim
# plot( seq(0,199,10), probs, xlab="day since introduction", ylab="extinction probability", main="R0=0.6, I0=10, Pop=1e4")
# lines( seq(0,199,10), probs )

