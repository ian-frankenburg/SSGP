FF = array(0, dim=c(S*sims, S, length(times)-1))
for(i in 1:dim(FF)[3]){
  l <- split(FF1[,,i], rep(1:ncol(FF1[,,i]), each=nrow(FF1[,,i])))
  FF[,,i] = as.matrix(Matrix::bdiag(l))
}
yin = matrix(0,nrow=dim(FF)[1], ncol=dim(ysp)[3])
for(i in 1:dim(ysp)[3]){
  yin[, i] = c(ysp[,,i])
}
z = log(yfull+1)
for(s in 1:S){
  matplot(t(ysp[,s,]),type="l",add=0)
  points((z[,s]), pch=19)
}
st_time = Sys.time();
p2=1
fullySpatial = spDLMGP(z=t(z)[,-1], y=yin[,-1], 
                       y_cube = array(ysp2[,-1,], dim=c(nrow(ysp2),ncol(ysp2)-1,dim(ysp2)[3])),
                       yinit = ysp2[1,1,],
                       niter=5000, burnin=5000/2,
                       FF=FF, FF_cube=FFsp,
                       C0=diag(2), W=matrix(c(1,.9,.9,1),2,2,byrow=TRUE),
                       params=theta_design,
                       delta_v=.9, delta_w=.95,
                       tune=rep(.001,ncol(theta_design)), tune2=rep(.1,ncol(theta_design)), 
                       nugget=1e-5, z_a=1, z_b = 1)

round(apply(fullySpatial$mc_C,1:2,mean),10)

##### FIT STAN #####
N <- length(lynx_hare_df$Year) - 1
ts <- 1:N
y_init <- c(lynx_hare_df$Hare[1], lynx_hare_df$Lynx[1])
y <- as.matrix(lynx_hare_df[2:(N + 1), 2:3])
y <- cbind(y[, 2], y[, 1])
# hare, lynx order
lynx_hare_data <- list(
  N = N,
  ts = ts,
  y0 = y_init,
  t0 = 0,
  y = y
)
library(cmdstanr)
# model = cmdstan_model(stan_file="LV_case_study.stan")
# fit=model$sample(data = lynx_hare_data,
#                  chains = 1,
#                  iter_warmup = 2500,
#                  iter_sampling = 2500, refresh = 10, parallel_chains = 4)
###### END FIT STAN #####

fit$summary("theta")
stan_draws = matrix(fit$draws("theta")[,1,],ncol=4)
library(latex2exp)
tempTest = matrix(0, nrow = length(stan_draws[,2]), ncol=4)
foo = unscaled
for(i in 1:4){
  tempTest[,i] = (stan_draws[,i]-min(foo[,i]))/(max(foo[,i])-min(foo[,i]))
}
par(mfrow=c(2,2))
for(i in 1:length(test)){
  drw = (rbeta(length(fullySpatial$calibrate[,i]),1,1))
  if(i==1){
    hist(tempTest[,i], xlim=c(0,1), breaks=10, probability = 1, xlab=TeX(
      paste0('$\\eta_', i, "$")), main="",
      col=rgb(1, 0, 0, alpha=.5),add=0)
    hist(fullySpatial$calibrate[,i], xlim=c(0,1), breaks=10, 
         col=rgb(0/255,0/255,.502,.7),
         probability = 1, xlab=TeX(
           paste0('$\\eta_', i, "$")), main="Calibration Posterior", add=1)
    hist(drw, xlim=c(0,1), breaks=10, probability = 1, xlab=TeX(
      paste0('$\\eta_', i, "$")), main="Calibration Posterior",
      col=rgb(.678, .847, .902, alpha=.4),add=T)
    # legend("topright", c("Prior", "Emulator/Calibrator", "Stan"), 
    #        col=c(rgb(.678, .847, .902),
    #               rgb(0/255,0/255,.502, .7),
    #               "red"), bty='n', text.font=2 ,cex=1, lty=1, lwd=10)
  }else{
    hist(tempTest[,i], xlim=c(0,1), breaks=10, probability = 1, xlab=TeX(
      paste0('$\\eta_', i, "$")), main="",
      col=rgb(1, 0, 0, alpha=.5),add=0)
    hist(fullySpatial$calibrate[,i], xlim=c(0,1), breaks=10, 
         col=rgb(0/255, 0/255, .502, .7),
         probability = 1, xlab=TeX(
           paste0('$\\eta_', i, "$")), main="", add=1)
    hist(drw, xlim=c(0,1), breaks=10, probability = 1, xlab=TeX(
      paste0('$\\eta_', i, "$")), main="",
      col=rgb(.678, .847, .902, alpha=.4),add=T)
  }
}

