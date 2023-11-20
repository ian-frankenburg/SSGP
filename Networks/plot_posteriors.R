for(i in 1:ncol(theta_design)){
  matplot((fullySpatial$mc_beta[i,]),type="l",
          main=paste("GP Parameter", i),ylab="y",xlab="MCMC draw")
}
matplot(fullySpatial$calibrate,type="l")
# hist((spgasp$mc_sigma2), xlim=c(0,2),
#      col=rgb(0/255,0/255,.502,.7), xlab="", breaks=50,
#      probability = 1, main=TeX('$\\sigma_z $ posterior'), add=0)
# hist(rgamma(nrow(spgasp$calibrate),1,1), xlim=c(0,2), breaks=50, probability = 1, xlab=TeX(
#   paste0('$\\eta_', i, "$")), main="Calibration Posterior",
#   col=rgb(.678, .847, .902, alpha=.5),add=T)
# abline(v=.2, col="red", lwd=3)
# apply(spgasp$calibrate, 2, acf, lag=200)
plot(sqrt(t(fullySpatial$mc_sigma2)),type="l")
# test
# apply(spgasp$calibrate,2,quantile)
test
save = fullySpatial$calibrate
round(colMeans(fullySpatial$calibrate),3)
{
  library(latex2exp)
  par(mfrow=c(2,2))
  for(i in 1:length(test)){
    drw = matrix(rbeta(length(fullySpatial$calibrate),1,1), ncol=ncol(fullySpatial$calibrate))
    if(i==1){
      hist(fullySpatial$calibrate[,i], breaks=10, main="", ylab="",
           col=rgb(0/255,0/255,.502,.7), xlim=c(0,1), cex.lab=2,
           freq = 0, xlab=TeX(
             paste0('$\\eta_', i, "$")), add=0)
      hist(ind_spatial$calibrate[,i], breaks=10, main="", ylab="",
           col=rgb(255/255,0/255,0,.5), xlim=c(0,1), cex.lab=2,
           freq = 0, xlab=TeX(
             paste0('$\\eta_', i, "$")), add=1)
      hist(drw[,i], breaks=10, freq = 0, xlab=TeX(
        paste0('$\\eta_', i, "$")),
        col=rgb(.678, .847, .902, alpha=.2),add=T)
    }else{
      hist(fullySpatial$calibrate[,i], breaks=10, main="",ylab="",
           col=rgb(0/255, 0/255, .502, .7), xlim=c(0,1), cex.lab=2,
           freq = 0, xlab=TeX(
             paste0('$\\eta_', i, "$")), add=0)
      hist(ind_spatial$calibrate[,i], breaks=10, main="", ylab="",
           col=rgb(255/255,0/255,0,.5), xlim=c(0,1), cex.lab=2,
           freq = 0, xlab=TeX(
             paste0('$\\eta_', i, "$")), add=1)
      hist(drw[,i], xlim=c(0,1), breaks=10, freq = 0, xlab=TeX(
        paste0('$\\eta_', i, "$")), main="",
        col=rgb(.678, .847, .902, alpha=.2),add=T)
    }
    points(y=0,x=test[i], col="red", lwd=5,pch=19)
  }
  hist(sqrt(fullySpatial$mc_sigma2), xlim=c(0,.25), cex.lab=2,ylab="",
       col=rgb(0/255,0/255,.502,.7), breaks=10, main="",
       freq = 0, xlab=TeX('$\\sigma_z $'), add=0)
  hist(sqrt(ind_spatial$mc_sigma2), breaks=10, main="", ylab="",
       col=rgb(255/255,0/255,0,.5), xlim=c(0,1), cex.lab=2,
       freq = 0, xlab=TeX(
         paste0('$\\eta_', i, "$")), add=1)
  hist(sqrt(rgamma(nrow(fullySpatial$calibrate),1,1)), xlim=c(0,2), breaks=10, freq = 0, xlab=TeX(
    paste0('$\\eta_', i, "$")),
    col=rgb(.678, .847, .902, alpha=.5),add=T)
  points(y=0,x=sigma_z, col="red", lwd=4,pch=19)
  plot(c(0,0), xlim=c(-10,10),ylim=c(-10,10),
       col="white",xlab="",ylab="",xaxt="n",yaxt="n",frame.plot = FALSE)
  legend(-10, 8, c("Prior", "Spatial Posterior", "Non-spatial Posterior"),
         fill=c(rgb(.678, .847, .902),
                rgb(0/255,0/255,.502, .7),
                rgb(1, 0, 0, alpha=.7)),
         bty = 'n',text.font=2,cex=1.5)
}
# points(rbind(test[1:2]),pch=19,cex=2,col="red")
# theta_loc = 7
# (test=theta_design[theta_loc,] + runif(4,0,0.05))
# goback = c()
# for(i in 1:ncol(theta_design)){
#   goback[i] = test[i]*(max(unscaled[,i])-min(unscaled[,i])) + min(unscaled[,i])
# }
# goback
# # 
# RcppParallel::setThreadOptions(numThreads = 4)
# a = goback[1]
# b = goback[2]
# c = goback[3]
# library(ggplot2)
# result <- spreadr(pnet, start_run, time = t_max,
#                   decay = a,
#                   retention = b,
#                   suppress = c,
#                   include_t0=F)
# ggplot(result, aes(x=time, y=activation, color=node)) +
#   geom_point() +
#   geom_line()