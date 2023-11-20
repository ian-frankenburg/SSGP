library(lhs)
library(mbr)
library(deSolve)
library(ReacTran)
library(latex2exp)
# https://arxiv.org/pdf/1310.0505.pdf
SIR <- function(t, y, parms) {
  S <- matrix(nrow = Nx, ncol = Ny, data = y[1:(Nx*Ny)])
  I <- matrix(nrow = Nx, ncol = Ny, data = y[(Nx*Ny+1) : (2*Nx*Ny)])
  dS <- -beta*S*I/N 
  dI <- beta*S*I/N - gamma*I +
    tran.2D(C = I, D.x = alpha1, D.y = alpha1,
            dx = Gridx, dy = Gridy,
            C.x.up = 0,
            C.y.up = 0,
            C.x.down = 0,
            C.y.down = 0)$dC
  list(c(dS, dI))
}
Nx <-90
Ny <-90
Gridx <- setup.grid.1D(x.up = 0, x.down = Nx, N = Nx)
Gridy <- setup.grid.1D(x.up = 0, x.down = Ny, N = Ny)
beta=.6; gamma=.2; alpha1=1; alpha2 = 0

N = 100000
X1ini <- matrix(nrow = Nx, ncol = Ny, data = N)
X2ini <- matrix(nrow = Nx, ncol = Ny, data = 0) 
X2ini[50,50] = 10


yini <- c(X1ini, X2ini)
times <- seq(0,100,by=10)
print(system.time(
  out <- ode.2D(y = yini, parms = NULL, func = SIR,
                nspec = 2, dimens = c(Nx, Ny), times = times,
                lrw = 2000000, names=c("X1", "X2"))
))
par(oma = c(0,0,1,0))
image(out, which = "X2", xlab = "x", ylab = "y",
      mfrow = c(3, 3), ask = FALSE,
      main = paste("t = ", times),
      grid = list(x = Gridx$x.mid, y = Gridy$x.mid))
mtext(side = 3, outer = TRUE, cex = 1.25, line = -1,
      "2-D SIR")
temp = expand.grid(x = 1:Nx, y = 1:Nx)
spPts = 15
locs = temp[round(qunif(maximinLHS(spPts,k = 1),1,Nx^2)),]
plot(locs, pch=19, col="red")
locs = as.numeric(rownames(locs))
sims=25
cube = maximinLHS(sims,k = 4)
set.seed(1)
b = qunif(cube[,1],.5,1)
set.seed(1)
g = qunif(cube[,2],0.01,.5)
set.seed(1)
a1 = qunif(cube[,3],0,.01)

S = length(locs)
newtimes = seq(1,nrow(out), by=1)
# ysp = array(0, dim=c(sims,length(newtimes),length(locs)))
ysp_old = array(0, dim=c(sims,length(newtimes),length(locs)))
ysp = array(0, dim=c(sims,S,length(times)))
ysp2 = array(0, dim=c(sims,length(times),S))
FF1 = array(0, dim=c(sims, S, length(times)-1))
FFc = array(0, dim=c(sims, 1, length(times)-1, S))
FF = array(0, dim=c(sims, 1, length(newtimes)-1, length(locs)))
theta_design = cbind(b, g, a1)
plot(theta_design)
for(i in 1:sims){
  beta = b[i]
  gamma = g[i]
  alpha1 = a1[i]
  out <- ode.2D(y = yini, parms = NULL, func = SIR,
                nspec = 2, dimens = c(Nx, Ny), times = times,
                lrw = 2000000, names=c("X1", "X2"))[,-1]
  # par(oma = c(0,0,1,0))
  # image(out, which = "X1", xlab = "x", ylab = "y",
  #       mfrow = c(3, 3), ask = FALSE,
  #       main = paste("t = ", times),
  #       grid = list(x = Gridx$x.mid, y = Gridy$x.mid))
  ysp_ = matrix((out[,(Nx*Ny+1) : (2*Nx*Ny)]), ncol=(Nx*Ny))
  # ysp_ = s[,2]
  ysp_ = log(ysp_+1)
  for(sp in 1:length(locs)){
    ysp_old[i,,sp] = (ysp_[,locs[sp]])
    FF[i,1,,sp] = (ysp_old[i,-length(newtimes),sp])
    
    ysp[i,sp,] = (ysp_[,locs[sp]])
    FF1[i,sp,] = (ysp[i,sp,-length(times)])
    
    ysp2[i,,sp] = (ysp_[,locs[sp]])
    FFc[i,1,,sp] = (ysp2[i,-length(times),sp])
  }
}
tempArr = array(0,dim=c(dim(FF)[1], dim(FF)[2], dim(FF)[3]))
FFsp = list()
for(i in 1:dim(ysp_old)[3]){
  for(j in 1:dim(tempArr)[3]){
    tempArr[,,j] = FF[,1,j,i]
  }
  FFsp[[i]] = tempArr
}
unscaled = theta_design
x = theta_design
for(i in 1:ncol(theta_design)){
  x[,i] = (theta_design[,i])
  temp = (x[,i]-min(x[,i]))/(max(x[,i])-min(x[,i]))
  att.temp = attributes(temp)
  theta_design[,i] = temp
  # center[i] = att.temp$`scaled:center`
  # scl[i] = att.temp$`scaled:scale`
}
matplot(t(ysp_old[,,1]),add=0,type="l")
for(i in 2:S){
  matplot(t(ysp_old[,,i]),add=T, type="l")
}
set.seed(1)
(test = c(runif(ncol(theta_design),.3,.7)))
pairs(rbind(theta_design,test),
      cex=c(rep(1, nrow(theta_design)), 2),
      pch=c(rep(1, nrow(theta_design)), 19), 
      col=c(rep("black",nrow(theta_design)),"red"),
      labels = c(TeX('$\\beta$'),TeX('$\\gamma$'), TeX('$\\alpha$')))
# (test=theta_design[4,] + runif(3,0,0.05))
goback = c()
for(i in 1:ncol(theta_design)){
  goback[i] = test[i]*(max(unscaled[,i])-min(unscaled[,i])) + min(unscaled[,i])
}
goback
# 
RcppParallel::setThreadOptions(numThreads = "auto", stackSize = "auto")
RcppParallel::defaultNumThreads()
beta = goback[1]
gamma= goback[2]
alpha1= goback[3]
alpha2= goback[4]
out <- ode.2D(y = yini, parms = list(alpha1=alpha1,
                                     alpha2=alpha2,
                                     beta=beta,
                                     gamma=gamma), func = SIR,
              nspec = 2, dimens = c(Nx, Ny), times = times,
              lrw = 2000000, names=c("X1", "X2"))[,-1]
zz = log(out[,(Nx*Ny+1) : (2*Nx*Ny)]+1)
sigma_z = 1
set.seed(1)
z = zz[,locs] + matrix(rnorm(length(dim(ysp)[2]*dim(ysp)[3]), 0, sigma_z), nrow(zz), length(locs))
matplot(z,type="l")
st_time = Sys.time();
p2 = 1
spgasp = parallel_mvcalibrator_discrepancy(z = z[-1,], 
                                           y = ysp_old[,-1,], 
                                           niter = 5000, burnin = 5000/2,
                                           m0=matrix(rep(0,p2)), C0=1*diag(p2), 
                                           FF = FFsp, G=1*diag(p2),
                                           params = theta_design,
                                           a0 = 1,b0=1,c0=1,
                                           tune=rep(.001,ncol(theta_design)), 
                                           tune2=rep(.001,ncol(theta_design)), 
                                           tune3 = .25,
                                           nugget = 0, S = S,
                                           alpha0 = 1, beta0 = 1, 
                                           yinit = ysp_old[1,1,],
                                           delta_w = .9,
                                           delta_vw = .1, delta_ww = .1,
                                           delta_v = .9, sigma_z=1,
                                           a0z = 1, b0z = 1, sptune=.1,
                                           discrep = 0);
matplot(spgasp$calibrate,type="l")

ztest=(z)
zrep = spgasp$mc_yeta
for(i in 1:length(spgasp$mc_sigma2)){
  zrep[,i,] = spgasp$mc_yeta[,i,] + 
    matrix(rnorm(length(z[-1,]),0,spgasp$mc_sigma2[i]), nrow(z)-1, ncol(z)) +
    spgasp$mc_sigma2[i]
}
mu=(apply(zrep,c(1,3),mean))
sd = (apply(zrep,c(1,3),sd))
matplot(mu,type="l",ylim=c(-10,20))
matplot(z[-1,],type="l",add=T,col="black",lty=1)
(GRS = -sum((z[-1,]-mu)^2)-2*sum(log(sd)))
# # 
# mu=apply(fullySpatial$yeta,c(1,3),mean)
# matplot(mu,type="l")
# matplot(z[-1,],type="l",col="black",add=T)
# sd=apply(fullySpatial$yeta,c(1,3),sd)
# (GRS = -sum(((ztest[-1,]-mu)/sd)^2)-2*sum(log(sd)))
# #
# mu=apply(nonSpatial$yeta,c(1,3),mean)
# matplot(mu,type="l")
# sd=apply(nonSpatial$yeta,c(1,3),sd)
# (GRS = -sum(((ztest[-1,]-mu)/sd)^2)-2*sum(log(sd)))


(Sys.time() - st_time)
# l = sample(1:dim(spgasp$mc_theta[[1]])[3],1)
# print(l)
# matplot(spgasp$mc_theta[[1]][,,l],type="l",add=0,ylim=c(-2,2))
# for(i in 2:2500){
#   matplot(spgasp$mc_theta[[i]][,,l],type="l",add=1,ylim=c(-2,2))
# }
# (Sys.time()-st_time)
# for(i in 1:ncol(theta_design)){
#   matplot((spgasp$mc_beta[i,]),type="l",
#           main=paste("GP Parameter", i),ylab="y",xlab="MCMC draw")
# }
matplot(spgasp$calibrate,type="l")
# hist((spgasp$mc_sigma2), xlim=c(0,2),
#      col=rgb(0/255,0/255,.502,.7), xlab="", breaks=50,
#      probability = 1, main=TeX('$\\sigma_z $ posterior'), add=0)
# hist(rgamma(nrow(spgasp$calibrate),1,1), xlim=c(0,2), breaks=50, probability = 1, xlab=TeX(
#   paste0('$\\eta_', i, "$")), main="Calibration Posterior",
#   col=rgb(.678, .847, .902, alpha=.5),add=T)
# abline(v=.2, col="red", lwd=3)
# apply(spgasp$calibrate, 2, acf, lag=200)
plot((t(spgasp$mc_sigma2)))
# test
# apply(spgasp$calibrate,2,quantile)
save = spgasp$calibrate
{
  library(latex2exp)
  par(mfrow=c(2,2))
  colMeans(spgasp$calibrate)
  for(i in 1:length(test[])){
    drw = matrix(rbeta(length(spgasp$calibrate),1,1), ncol=ncol(spgasp$calibrate))
    if(i==1){
      hist(spgasp$calibrate[,i], xlim=c(0,1), breaks=20, cex.main=2,
           col=rgb(0/255,0/255,.502,.7),
           freq = 0, xlab=TeX(
             paste0('$\\eta_', i, "$")),
           main="Calibration Parameters", add=0)
      hist(drw[,i], xlim=c(0,1), breaks=10, freq = 0, xlab=TeX(
        paste0('$\\eta_', i, "$")),
        col=rgb(.678, .847, .902, alpha=.2),add=T)
    }else{
      hist(spgasp$calibrate[,i], xlim=c(0,1), breaks=20, 
           col=rgb(0/255, 0/255, .502, .7),
           freq = 0, xlab=TeX(
             paste0('$\\eta_', i, "$")), main="", add=0)
      hist(drw[,i], xlim=c(0,1), breaks=20, freq = 0, xlab=TeX(
        paste0('$\\eta_', i, "$")), main="",
        col=rgb(.678, .847, .902, alpha=.2),add=T)
    }
    abline(v=test[i], col="red", lwd=3)
  }
  hist(sqrt(spgasp$mc_sigma2), xlim=c(0,2),
       col=rgb(0/255,0/255,.502,.7), breaks=20, main="",
       freq = 0, xlab=TeX('$\\sigma_z $'), add=0)
  hist(sqrt(rgamma(nrow(spgasp$calibrate),1,1)), xlim=c(0,2),
       breaks=10, freq = 0, xlab=TeX(
         paste0('$\\eta_', i, "$")),
       col=rgb(.678, .847, .902, alpha=.5),add=T)
  points(y=0,x=sigma_z, col="red", lwd=4)
}



# {
#   library(latex2exp)
#   par(mfrow=c(2,2))
#   colMeans(spgasp$calibrate)
#   for(i in 1:length(test[])){
#     if(i==1){
#       hist(spgasp$mc_beta[i,], breaks=15, 
#            col=rgb(0/255,0/255,.502,.7),
#            freq = 0, xlab=TeX(
#              paste0('$\\beta_', i, "$")), main="GP Correlation Parameters", add=0)
#       legend(.1, 6, c("Prior", "Posterior", "Truth"), 
#              fill=c(rgb(.678, .847, .902),
#                     rgb(0/255,0/255,.502, .7),
#                     "red"))
#     }else{
#       hist(spgasp$mc_beta[i,], breaks=15, 
#            col=rgb(0/255, 0/255, .502, .7),
#            freq = 0, xlab=TeX(
#              paste0('$\\beta_', i, "$")), main="", add=0)
#     }
#   }
# }

# matplot(apply(spgasp$mc_w,1:2,median), type="l")
# matplot(z[-1,],type="l")
# matplot(apply(spgasp$mc_yeta,3,rowMeans)+apply(spgasp$mc_w,1:2,mean)[-1,],type="l",add=T)
# matplot(apply(spgasp$mc_yeta,3,rowMeans),type="l",add=T)

# 
# matplot(spgasp$calibrate,type="l")
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
#   beta_draws = spgasp$mc_beta
#   theta_draws = spgasp$mc_theta
#   K_draws = spgasp$mc_K
#   v_draws = spgasp$v
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
# matplot(log(simulation(goback[1,1],goback[1,2],1:60,N=N)[,2]),type="l",col="darkred",lwd=2,add=1)
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