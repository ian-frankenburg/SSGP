S2=1
FF = array(0, dim=c(S*sims, S2, length(times)-1))
for(i in 1:dim(FF)[3]){
  l <- split(FF1[,,i], rep(1:ncol(FF1[,,i]), each=nrow(FF1[,,i])))
  FF[,,i] = as.matrix(unlist(l))
}
yin = matrix(0,nrow=dim(FF)[1], ncol=dim(ysp)[3])
for(i in 1:dim(ysp)[3]){
  yin[, i] = c(ysp[,,i])
}
p2 = 1
st_time = Sys.time();
nonSpatial = DLMGP(z=z[,-1], y=yin[,-1],
                          yy = array(ysp2[,-1,],
                          dim=c(nrow(ysp2),ncol(ysp2)-1,dim(ysp2)[3])),
                          FF=FF, FFc=FFsp,
                          m0=matrix(rep(0,1)), 
                          C0=matrix(1), W=matrix(1),
                          yinit = ysp[1,,1], alpha0=.1, beta0=.1,
                          delta_v=.95, delta_w=.95,
                          params=theta_design,
                          pred_params = test,
                          tune=rep(.001,3), tune2=rep(.001,3), 
                          tune3=.25, tune4=rep(.01,2),
                          niter=5000, burnin=5000/2,
                          nugget=0, sigma_z=rep(1,1),
                          a0z = 1, b0z=1, S=1, S2=S)
(Sys.time()-st_time);
ztest = z
mu=apply(nonSpatial$yeta,c(1,3),mean)
matplot(mu,type="l")
sd=apply(nonSpatial$yeta,c(1,3),sd)
(GRS = -sum(((ztest[-1,]-mu)/sd)^2)-2*sum(log(sd)))


# dim(apply(nonSpatial$yeta,c(1,3),mean))
MLmetrics::RMSE(apply(nonSpatial$yeta,c(1,3),mean), (z[-1,]))
matplot((apply(nonSpatial$yeta,c(1,3),mean)), type="l")


heatmap(apply(nonSpatial$mc_C,1:2,mean), Rowv = NA,Colv = NA)
(Sys.time()-st_time)
for(i in 1:ncol(theta_design)){
  matplot((nonSpatial$mc_beta[i,]),type="l",
          main=paste("GP Parameter", i),ylab="y",xlab="MCMC draw")
}
matplot(nonSpatial$calibrate,type="l")
library(latex2exp)
# hist((nonSpatial$mc_sigma2), xlim=c(0,2),
#      col=rgb(0/255,0/255,.502,.7), xlab="", breaks=50,
#      probability = 1, main=TeX('$\\sigma_z $ posterior'), add=0)
# hist(rgamma(nrow(nonSpatial$calibrate),1,1), xlim=c(0,2), breaks=50, probability = 1, xlab=TeX(
#   paste0('$\\eta_', i, "$")), main="Calibration Posterior",
#   col=rgb(.678, .847, .902, alpha=.5),add=T)
# abline(v=sigma_z, col="red", lwd=3)
# plot(t(nonSpatial$mc_sigma2))
# test
# apply(nonSpatial$calibrate,2,quantile)
{
library(latex2exp)
par(mfrow=c(2,2))
colMeans(nonSpatial$calibrate)
for(i in 1:length(test[])){
  drw = matrix(rbeta(length(nonSpatial$calibrate),1,1), ncol=ncol(nonSpatial$calibrate))
  if(i==1){
    hist(nonSpatial$calibrate[,i], xlim=c(0,1), breaks=15, 
         col=rgb(0/255,0/255,.502,.7),
         freq = 0, xlab=TeX(
           paste0('$\\eta_', i, "$")), main="Calibration Posterior", add=0)
    hist(drw[,i], xlim=c(0,1), breaks=10, freq = 0, xlab=TeX(
         paste0('$\\eta_', i, "$")), main="Calibration Posterior",
       col=rgb(.678, .847, .902, alpha=.2),add=T)
    # legend(.1, 6, c("Prior", "Posterior", "Truth"), 
    #        fill=c(rgb(.678, .847, .902),
    #               rgb(0/255,0/255,.502, .7),
    #               "red"))
  }else{
    hist(nonSpatial$calibrate[,i], xlim=c(0,1), breaks=15, 
         col=rgb(0/255, 0/255, .502, .7),
         freq = 0, xlab=TeX(
           paste0('$\\eta_', i, "$")), main="", add=0)
    hist(drw[,i], xlim=c(0,1), breaks=10, freq = 0, xlab=TeX(
      paste0('$\\eta_', i, "$")), main="",
      col=rgb(.678, .847, .902, alpha=.2),add=T)
  }
  abline(v=test[i], col="red", lwd=3)
}
hist((nonSpatial$mc_sigma2),
     col=rgb(0/255,0/255,.502,.7), breaks=30, main="",
     freq = 0, xlab=TeX('$\\sigma_z $'), add=0)
hist(rgamma(nrow(nonSpatial$calibrate),1,1), xlim=c(0,2), breaks=50, freq = 0, xlab=TeX(
  paste0('$\\eta_', i, "$")),
  col=rgb(.678, .847, .902, alpha=.5),add=T)
abline(v=sigma_z, col="red", lwd=3)
}



{
  library(latex2exp)
  par(mfrow=c(2,2))
  colMeans(nonSpatial$calibrate)
  for(i in 1:length(test[])){
    if(i==1){
      hist(nonSpatial$mc_beta[i,], breaks=15, 
           col=rgb(0/255,0/255,.502,.7),
           freq = 0, xlab=TeX(
             paste0('$\\beta_', i, "$")), main="GP Correlation Parameters", add=0)
      legend(.1, 6, c("Prior", "Posterior", "Truth"), 
             fill=c(rgb(.678, .847, .902),
                    rgb(0/255,0/255,.502, .7),
                    "red"))
    }else{
      hist(nonSpatial$mc_beta[i,], breaks=15, 
           col=rgb(0/255, 0/255, .502, .7),
           freq = 0, xlab=TeX(
             paste0('$\\beta_', i, "$")), main="", add=0)
    }
  }
}

# matplot(apply(nonSpatial$mc_w,1:2,median), type="l")
# matplot(z[-1,],type="l")
# matplot(apply(nonSpatial$mc_yeta,3,rowMeans)+apply(nonSpatial$mc_w,1:2,mean)[-1,],type="l",add=T)
# matplot(apply(nonSpatial$mc_yeta,3,rowMeans),type="l",add=T)

# 
# matplot(nonSpatial$calibrate,type="l")
# par(oma = c(0,0,1,0))
# {
#   nrep=1000
#   emloc = round(runif(1,1,dim(ysp)[3]))
#   thetaloc = round(runif(1,1,nrow(theta_design)))
#   emulated = matrix(0,nrow=ncol(ysp)-1,ncol=nrep)
#   (test = c(runif(3,.3,.7)))
#   plot(theta_design[,1:2])
#   points(rbind(test[1:2]),pch=19,cex=2,col="red")
#   #(test=theta_design[thetaloc,] + runif(3,0,0.1))
#   goback = c()
#   for(i in 1:3){
#     goback[i] = test[i]*(max(unscaled[,i])-min(unscaled[,i])) + min(unscaled[,i])
#   }
#   beta_draws = nonSpatial$mc_beta
#   theta_draws = nonSpatial$mc_theta
#   K_draws = nonSpatial$mc_K
#   v_draws = nonSpatial$v
#   for(i in 1:nrep){
#     temp = spemulate(ysp[,-1,emloc], pred_params = test,
#                      params = theta_design, 
#                      beta = beta_draws[,sample(1:ncol(beta_draws),1)],
#                      theta = matrix(theta_draws[[sample(1:ncol(beta_draws),1)]][,,emloc],nrow=1), 
#                      F = FFsp[[emloc]],nugget=0,
#                      K = K_draws[,,sample(1:dim(K_draws)[3],1)]+0*diag(nrow(theta_design)), 
#                      yinit = ysp[thetaloc,1,emloc], 
#                      v_draws[,emloc,sample(1:ncol(beta_draws),1)])
#     emulated[,i] = temp$ypred
#   }
#   # matplot((emulated[,colSums(emulated<0)<=0])^1,type="l",add=0,col=rgb(0,0,0,alpha=1),
#   # main=paste("Spatial Location", emloc))
#   emulated = emulated[,colSums(emulated<0)<=0]
#   q05 = apply(emulated,1,quantile,.05,na.rm=T)
#   q5 = apply(emulated,1,quantile,.5,na.rm=T)
#   q95 = apply(emulated,1,quantile,.95,na.rm=T)
#   # s = simulation(beta=goback[1],gamma = goback[2],x = 1:time, N = N)[,]
#   beta = goback[1]
#   gamma = goback[2]
#   alpha = goback[3]
#   out <- ode.2D(y = yini, parms = NULL, func = SIR,
#                 nspec = 2, dimens = c(Nx, Ny), times = times,
#                 lrw = 2000000, names=c("X1", "X2"))
#   out = out[newtimes,]
#   matplot((out[-1,S+1+emloc])^(1/2),col="lightblue",type="l",lwd=3,
#           main=paste("Spatial Location"),ylim=c(0,max((out[-1,S+1+emloc])^(1/2))+10))
#   matplot((q05),type="l",col="darkorange",lwd=3,add=1, 
#           main=paste("Spatial Location"))
#   matplot((q5),type="l",col="darkorange",lwd=3,add=1)
#   matplot((q95),type="l",col="darkorange",lwd=3,add=T)
#   matplot((out[-1,S+1+emloc])^(1/2),col="darkblue",type="l",lwd=1,add=T,
#           main=paste("Spatial Location"),ylim=c(0,max((out[-1,S+1+emloc])^(1/2))+10))
#   # matplot((s[-1,2])^(1/4),add=1,col="darkorange",lwd=3,type="l")
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# library(cmdstanr)
# model = cmdstan_model(stan_file="marginalizedEmulator_DiscountFactors.stan")
# fit=model$sample(data = list("T"=ncol(yin[,]),
#                              num_series=nrow(yin[,]),
#                              'y' = yin[,],
#                              p=2, inits=as.array(1),
#                              theta=theta_design,
#                              theta_pred=test,
#                              lags=1, yinit=1,
#                              F=FF), chains = 1, iter_warmup = 500, 
#                  iter_sampling = 500, refresh = 1)
# fit$summary("beta")
# matplot(fit$summary("ypred")$mean,col="darkblue",type="l",lwd=6)
# matplot(fit$summary("ypred")$q5,col="darkblue",type="l",lwd=6,add=0 )
# matplot(fit$summary("ypred")$q95,col="darkblue",type="l",lwd=6,add=T )
# 
# 
# matplot(t(yin)[,1],type="l",col="pink",lwd=3,add=T)
# matplot(sqrt(simulation(goback[1,1],goback[1,2],1:60,N=N)[,2]),type="l",col="darkred",lwd=2,add=1)
# 
# colMeans(t(gasp$mc_beta))
# fit$summary("beta")$mean
# 
# fit$summary("tau")$mean
# plot(fit$summary("phi")$median)
# matplot(matrix(fit$summary("errors")$median,nrow=ncol(yin[,-1])),type="l")
# 
# matplot(t(matrix(fit$summary("meanfunc")$mean,nrow=nrow(theta_design))),type="l")
# matplot(t(yin),type="l",col="black",add=1)
# 
# 
# 
# library(cmdstanr)
# model = cmdstan_model(stan_file="stantest.stan")
# fit=model$sample(data = list("T"=ncol(yin[,]),
#                              num_series=nrow(yin[,]),
#                              'y' = yin[,],
#                              p=4, inits=as.array(1),
#                              theta=theta_design,
#                              theta_pred=test,
#                              lags=1, yinit=1,
#                              F=FF), chains = 1, iter_warmup = 200, 
#                  iter_sampling = 200, refresh = 1)