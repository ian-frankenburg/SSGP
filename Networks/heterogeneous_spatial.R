st_time = Sys.time();
p2=1
ind_spatial = parallel_mvcalibrator_discrepancy(z = z[-1,], 
                           y = ysp2[,-1,], 
                           niter = 5000, burnin = 5000/2,
                           m0 = matrix(rep(0,p2)), 
                           C0=diag(p2), 
                           FF = FFsp, G=1*diag(p2),
                           params = theta_design,
                           a0 = 1, b0=1, c0=1,
                           tune = rep(.01, ncol(theta_design)), 
                           tune2 = rep(.01, ncol(theta_design)), 
                           tune3 = .5,
                           nugget = 1e-6, S = S,
                           alpha0 = 1, beta0 = 1, 
                           yinit = ysp2[1,1,],
                           delta_w = .9,
                           delta_vw = .1, delta_ww = .1,
                           delta_v = .9,
                           a0z = 1, b0z = 1, sptune=.5,
                           discrep = 0, sigma_z = .1);
# zrep = ind_spatial$mc_yeta
# for(i in 1:ncol(ind_spatial$mc_sigma2)){
#   zrep[,i,] = ind_spatial$mc_yeta[,i,] + ind_spatial$mc_sigma2[i]
# }
# mu=(apply(zrep,c(1,3),mean))
# sd = (apply(zrep,c(1,3),sd))
# # matplot(mu,type="l")
# # matplot(z,type="l")
# MLmetrics::RMSE(z[-1,], mu)
# (GRS = -sum((z[-1,]-mu)^2)-2*sum(log(sd)))
# # # 
# zrep = fullySpatial$yeta
# for(i in 1:ncol(fullySpatial$mc_sigma2)){
#   zrep[,i,] = fullySpatial$yeta[,i,] + fullySpatial$mc_sigma2[i]
# }
# mu=(apply(zrep,c(1,3),mean))
# sd = (apply(zrep,c(1,3),sd))
# # matplot(mu,type="l")
# # matplot(z,type="l")
# RMSE(z[-1,], mu)
# (GRS = -sum((z[-1,]-mu)^2)-2*sum(log(sd)))
# # # 
# zrep = nonspatial$yeta
# for(i in 1:ncol(nonspatial$mc_sigma2)){
#   zrep[,i,] = nonspatial$yeta[,i,] + nonspatial$mc_sigma2[i]
# }
# mu=(apply(zrep,c(1,3),mean))
# sd = (apply(zrep,c(1,3),sd))
# matplot(mu,type="l")
# matplot(z,type="l")
(GRS = -sum((z[-1,]-mu)^2)-2*sum(log(sd)))
# # 

(Sys.time()-st_time)
for(i in 1:ncol(theta_design)){
  matplot((ind_spatial$mc_beta[i,]),type="l",
          main=paste("GP Parameter", i),ylab="y",xlab="MCMC draw")
}
matplot(ind_spatial$calibrate,type="l",ylim=c(0,1))
# hist((ind_spatial$mc_sigma2), xlim=c(0,2),
#      col=rgb(0/255,0/255,.502,.7), xlab="", breaks=50,
#      probability = 1, main=TeX('$\\sigma_z $ posterior'), add=0)
# hist(rgamma(nrow(ind_spatial$calibrate),1,1), xlim=c(0,2), breaks=50, probability = 1, xlab=TeX(
#   paste0('$\\eta_', i, "$")), main="Calibration Posterior",
#   col=rgb(.678, .847, .902, alpha=.5),add=T)
# abline(v=.2, col="red", lwd=3)
# apply(ind_spatial$calibrate, 2, acf, lag=200)
# plot(t(ind_spatial$mc_sigma2))
# test
# apply(ind_spatial$calibrate,2,quantile)
test
round(colMeans(ind_spatial$calibrate),3)
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
      hist(ind_spatial$calibrate[,i],, breaks=10, main="", ylab="",
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
      hist(ind_spatial$calibrate[,i],, breaks=10, main="", ylab="",
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
  hist(sqrt(ind_spatial$mc_sigma2),, breaks=10, main="", ylab="",
       col=rgb(255/255,0/255,0,.5), xlim=c(0,1), cex.lab=2,
       freq = 0, xlab=TeX(
         paste0('$\\eta_', i, "$")), add=1)
  hist(sqrt(rgamma(nrow(fullySpatial$calibrate),1,1)), xlim=c(0,2), breaks=10, freq = 0, xlab=TeX(
    paste0('$\\eta_', i, "$")),
    col=rgb(.678, .847, .902, alpha=.5),add=T)
  points(y=0,x=sigma_z, col="red", lwd=4,pch=19)
  plot(c(0,0),col="white",xlab="",ylab="",xaxt="n",yaxt="n",frame.plot = FALSE)
  legend(1, 1, c("Prior", "Spatial", "Non-spatial"),
         fill=c(rgb(.678, .847, .902),
                rgb(0/255,0/255,.502, .7),
                rgb(1, 0, 0, alpha=.7)),
         bty = 'n',text.font=2,cex=1.5)
}
# points(rbind(test[1:2]),pch=19,cex=2,col="red")
# theta_loc = 7
# (test=theta_design[theta_loc,] + runif(4,0,0.05))
goback = c()
for(i in 1:ncol(theta_design)){
  goback[i] = test[i]*(max(unscaled[,i])-min(unscaled[,i])) + min(unscaled[,i])
}
goback
# 
RcppParallel::setThreadOptions(numThreads = 4)
a = goback[1]
b = goback[2]
c = goback[3]
library(ggplot2)
result <- spreadr(pnet, start_run, time = t_max,
                  decay = a,
                  retention = b,
                  suppress = c,
                  include_t0=F)
# ggplot(result, aes(x=time, y=activation, color=node)) +
#   geom_point() +
#   geom_line()