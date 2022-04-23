set.seed(1)
library(MLmetrics)
LotVmod <- function (Time, State, Pars){
  with(as.list(c(State, Pars)), {
    dx = x*(alpha - beta*y)
    dy = -y*(gamma - delta*x)
    return(list(c(dx, dy)))
  })
}
names = c("alpha","beta","gamma","delta")
Pars = list()
test <- c(alpha = .55, beta = .025, gamma = .8, delta = .024)
# Pars2 = list()
# for(i in 1:4){
#   Pars[names[i]] = foo[i]*(max(unscaled[,i])-min(unscaled[,i])) + min(unscaled[,i])
#   Pars2[names[i]] = test[i]#fit$summary("theta")$median[i]
#   print(foo2[,i]*(max(unscaled[,i])-min(unscaled[,i])) + min(unscaled[,i]))
# }
State <- c(x = 30, y = 4)
Time <- seq(0, 20, by = 1)
# # initial state 
# out <- as.data.frame(ode(func = LotVmod, y = State, parms = Pars, times = Time))[,-1]
# matplot(log(out),type="l",ylim=c(0,10),col=c("orange","blue"),lty=1,lwd=3)
# matplot(log(rbind(y_init,y)),type="p",pch=19,add=T)
# MAE(yfull[,1], out[,1]);MSE(yfull[,2], out[,2])
# 
# out <- as.data.frame(ode(func = LotVmod, y = State, parms = Pars2, times = Time))[,-1]
# matplot(log(out),type="l",ylim=c(0,10),add=0,col=c("gold","darkblue"),lty=1,lwd=3)
# matplot(log(rbind(y_init,y)),type="p",pch=19,add=T)
# MAE(yfull[,1], out[,1]);MSE(yfull[,2], out[,2])


# x_design = sample(0:10,size = nsims,replace = T) # field observations
# x_design = sample(0:10,size = nsims,replace = T) # field observations
library(lhs)
library(truncnorm)
sims=50
A <- maximinLHS(sims, 4) 
theta_design = matrix(0,nrow=sims,ncol=4)
# design = fit$summary("theta")
test=c(.55,.028,.84,.024)
d1 = .6
d2 = .015
theta_design[,1] <- qtruncnorm(A[,1], a=0,b=Inf,mean = 1,sd=.5)
theta_design[,2] <- qtruncnorm(A[,2], a=0,b=Inf,mean = .05,sd=.05)
theta_design[,3] <- qtruncnorm(A[,3], a=0,b=Inf,mean = 1,sd=.5)
theta_design[,4] <- qtruncnorm(A[,4], a=0,b=Inf,mean = .05,sd=.05)

# theta_design[,1] <- qunif(A[,1], 0.01, 1.5)
# theta_design[,2] <- qunif(A[,2], 0.01, .1)
# theta_design[,3] <- qunif(A[,3], 0.01, 1.5)
# theta_design[,4] <- qunif(A[,4], 0.01, .1)
# apply(theta_design, 2, hist, xlim=c(0,1), ylim=c(0,1))

colMeans(theta_design)
# apply(theta_design,2,hist, breaks=100)

library(mbr)
library(deSolve)
lynx_hare_df <- read.csv("lynx_hare_df.csv")
N <- length(lynx_hare_df$Year) - 1
ts <- 1:N
y_init <- c(lynx_hare_df$Hare[1], lynx_hare_df$Lynx[1])
y <- as.matrix(lynx_hare_df[2:(N + 1), 2:3])
y <- cbind(y[ , 2], y[ , 1]); # hare, lynx order
yfull = rbind(y_init,y)
lynx_hare_data <- list(N = N, ts = ts, y_init = y_init, y = y)

N = nrow(yfull)
time=nrow(lynx_hare_df)
times = seq(1, time, by = 1)
times = Time
p1=matrix(0,nrow(theta_design), length(times))
p2=matrix(0,nrow(theta_design), length(times))
Nx = 2; Ny = 1
S = Nx*Ny
ysp = array(0, dim=c(sims,length(times),Nx*Ny))
FF = array(0, dim=c(sims, 1, length(times)-1, Nx*Ny))
for(i in 1:sims){
  Pars <- c(alpha = theta_design[i,1], 
            beta = theta_design[i,2], 
            gamma = theta_design[i,3],
            delta = theta_design[i,4])
  out <- as.data.frame(ode(func = LotVmod, y = State, parms = Pars, times = Time))[,-1]+1
  ysp_ = log(out+1)
  for(sp in 1:(Nx*Ny)){
    ysp[i,,sp] = (ysp_[,sp])
    FF[i,1,,sp] = (ysp[i,-length(times),sp])
  }
}
l=seq(1,length(Time),length.out=max(Time)+1)
for(s in 1:dim(ysp)[3]){
  matplot(t(ysp[,,s]), type="l", lwd=2, lty=1, main="Experimental Design", 
          ylab="log-population size", xlab="Year", xaxt="n")
  axis(1, at = seq(1,length(Time),length.out=max(Time)), 
       labels = seq(1,20),
       las=2)
  points(round(l),log(yfull[,s]+1), pch=19, cex=2)
}
tempArr = array(0,dim=c(dim(FF)[1], dim(FF)[2], dim(FF)[3]))
FFsp = list()
for(i in 1:dim(ysp)[3]){
  for(j in 1:dim(tempArr)[3]){
    tempArr[,,j] = FF[,1,j,i]
  }
  FFsp[[i]] = tempArr
}
library(LICORS)
unscaled = theta_design
x = theta_design
tempTest = test
for(i in 1:ncol(theta_design)){
  x[,i] = (theta_design[,i])
  temp = (x[,i]-min(x[,i]))/(max(x[,i])-min(x[,i]))
  tempTest[i] = (test[i]-min(x[,i]))/(max(x[,i])-min(x[,i]))
  theta_design[,i] = temp
}
tempTest

p=nrow(c); p2 = ncol(FF)
start_time = Sys.time()
swtch = 250
RcppParallel::setThreadOptions(numThreads = 4)
# B <- c(maximinLHS(1, 4))
# theta_test = c()
# theta_test[1] <- qtruncnorm(B[1], a=0,b=Inf, mean = test[1], sd = .05)
# theta_test[2] <- qtruncnorm(B[2], a=0,b=Inf, mean = test[2], sd = .05)
# theta_test[3] <- qtruncnorm(B[3], a=0,b=Inf, mean = test[3], sd = .05)
# theta_test[4] <- qtruncnorm(B[4], a=0,b=Inf, mean = test[4], sd = .05)
# Pars <- c(alpha = theta_test[1],
#           beta = theta_test[2],
#           gamma = theta_test[3],
#           delta = theta_test[4])
# out <- as.data.frame(ode(func = LotVmod, y = State, parms = Pars, times = Time))[,-1]
# z = log(out+1)
z = log(yfull+1)
for(s in 1:dim(ysp)[3]){
  matplot(t(ysp[,,s]),type="l",add=0)
  points((z[,s]), pch=19)
}
st_time = Sys.time();

# spgasp = parallel_mvcalibrator(z=matrix(z[-1,],ncol=dim(ysp)[3]), 
#                                   y = array(ysp[,-1,],dim=c(nrow(ysp),ncol(ysp)-1,dim(ysp)[3])), 
#                                   niter = 5000, burnin=5000/2,
#                                   m0=matrix(rep(0,p2)), C0=diag(p2), 
#                                   FF = FFsp, G=diag(p2),
#                                   params = theta_design,
#                                   a0 = 1,b0=1,c0=1,
#                                   tune=rep(.001,4), tune2=rep(.001,4), tune3 = .1,
#                                   nugget = 0, S = dim(ysp)[3],
#                                   alpha0 = 1, beta0 = 1,yinit = y_init,
#                                   delta_w = .99,
#                                   delta_v = .99,
#                                   a0z = 1, b0z = 1)
p2=1
spgasp = parallel_mvcalibrator_discrepancy(z=as.matrix(z[-1,],ncol=2), 
             y = array(ysp[,-1,],dim=c(nrow(ysp),ncol(ysp)-1,dim(ysp)[3])), 
             niter = 10000, burnin = 10000/2,
             m0=matrix(rep(0,p2)), C0=diag(p2), 
             FF = FFsp, G=1*diag(p2),
             params = theta_design,
             a0 = 1,b0=1,c0=1,
             tune=rep(.001,4), tune2=rep(.001,4), tune3 = .1,
             nugget = 0, S = 2,
             alpha0 = 1, beta0 = 1, yinit = ysp[1,1,],
             delta_w = .9, delta_vw = .8, delta_ww = .8,
             delta_v = .9, sigma_z=1,
             a0z = 1, b0z = 1, sptune=.1, discrep = 0);


v = apply(spgasp$v,1:2,median)
matplot(v)

matplot((apply(spgasp$mc_w,1:2,median)),type="l",ylim=c(0,1))
matplot((apply(spgasp$mc_w,1:2,quantile,.95)),type="l",add=T)
matplot((apply(spgasp$mc_w,1:2,quantile,.5)),type="l",add=T)

dim(spgasp$mc_theta[[1]])
matplot(spgasp$mc_theta[[1]][,,2],type="l",add=0,ylim=c(-2,2))
for(i in 2:2500){
  matplot(spgasp$mc_theta[[i]][,,2],type="l",add=1,ylim=c(-2,2))
}

library(latex2exp)
tempTest = matrix(0, nrow = length(stan_draws[,i]), ncol=4)
foo = unscaled
for(i in 1:4){
  tempTest[,i] = (stan_draws[,i]-min(foo[,i]))/(max(foo[,i])-min(foo[,i]))
}
par(mfrow=c(2,2))
for(i in 1:length(test)){
  drw = (rbeta(length(spgasp$calibrate[,i]),1,1))
  if(i==1){
    hist(tempTest[,i], xlim=c(0,1), breaks=10, probability = 1, xlab=TeX(
      paste0('$\\eta_', i, "$")), main="",
      col=rgb(1, 0, 0, alpha=.5),add=0)
    hist(spgasp$calibrate[,i], xlim=c(0,1), breaks=10, 
         col=rgb(0/255,0/255,.502,.7),
         probability = 1, xlab=TeX(
           paste0('$\\eta_', i, "$")), main="Calibration Posterior", add=1)
    hist(drw, xlim=c(0,1), breaks=10, probability = 1, xlab=TeX(
      paste0('$\\eta_', i, "$")), main="Calibration Posterior",
      col=rgb(.678, .847, .902, alpha=.4),add=T)
    legend(.35, 10, c("Prior", "spDLMGP Posterior", "Stan Posterior"), 
           fill=c(rgb(.678, .847, .902),
                  rgb(0/255,0/255,.502, .7), 
                  "red"), cex=.6)
  }else{
    hist(tempTest[,i], xlim=c(0,1), breaks=10, probability = 1, xlab=TeX(
      paste0('$\\eta_', i, "$")), main="",
      col=rgb(1, 0, 0, alpha=.5),add=0)
    hist(spgasp$calibrate[,i], xlim=c(0,1), breaks=10, 
         col=rgb(0/255, 0/255, .502, .7),
         probability = 1, xlab=TeX(
           paste0('$\\eta_', i, "$")), main="", add=1)
    hist(drw, xlim=c(0,1), breaks=10, probability = 1, xlab=TeX(
      paste0('$\\eta_', i, "$")), main="",
      col=rgb(.678, .847, .902, alpha=.4),add=T)
  }
}


# matplot((spgasp$thetaw[,,1]),type="l")
# matplot(rowMeans(spgasp$mc_Vw[,1,]),type="l")
# matplot(rowMeans(spgasp$mc_Vw[,2,]),type="l")
# 
# matplot(t(spgasp$thetaw[,,1])%*%matrix(c(1,spgasp$calibrate[sample(1:250,1),])),type="l",add=1,ylim=c(-20,20))
# matplot(spgasp$mc_w[,2,sample(1:2500,1)],type="l",ylim=c(-5,5),add=T)
# for(i in 1:250){
#   matplot(spgasp$mc_w[,2,sample(1:2500,1)],type="l",add=T)
# }


# matplot(t(spgasp$mc_beta),type='l')
# matplot(t(spgasp$v),type='l')


# matplot(t(apply(spgasp$v,1:2,quantile,c(.975,.5,.025))[,,1]),type="l")
# matplot(t(apply(spgasp$v,1:2,quantile,c(.975,.5,.025))[,,2]),type="l",add=T)
# matplot(sqrt(t(spgasp$mc_sigma2)), type="l")
# 
# 
# # 
# matplot(t(apply(spgasp$mc_w,1:2,median)),type="l",ylim=c(-5,5))
# matplot(t(apply(spgasp$mc_w,1:2,quantile,.025)),type="l",add=1)
# matplot(t(apply(spgasp$mc_w,1:2,quantile,.975)),type="l",add=1)
# points(z-z2, pch=19, col="darkred")
# 
# (Sys.time() - st_time)
# apply(spgasp$mc_Vw,1:2,median)
# apply(spgasp$mc_Vw,1:2,mean)
# apply(spgasp$mc_Vw,1:2,median)

# apply(t(spgasp$mc_beta), 2, matplot, type="l")
{
  thinCal = spgasp$calibrate[seq(1,nrow(spgasp$calibrate), by=1),]
  apply(thinCal,2,median)
  tempTest
  matplot(thinCal,type="l")
  goback = c()
  foo = rbind(test,unscaled)
  for(i in 1:4){
    tempTest[i] = (test[i]-min(foo[,i]))/(max(foo[,i])-min(foo[,i]))
    hist(thinCal[,i],breaks=10,xlim=c(0,1), add=0, probability = 0)
    hist(rbeta(nrow(thinCal), 1, 1), breaks=10,add=1,
         col="darkblue",xlim=c(0,1), ylim=c(0,1))
    abline(v=tempTest[1],lwd=3,col="darkred")
  }
}


# names = c("alpha","beta","gamma","delta")
# (foo = apply(thinCal,2,mean))
# (foo2 = apply(thinCal,2,quantile,c(.1,.5,.9)))
# Pars = list()
# test <- c(alpha = .55, beta = .025, gamma = .8, delta = .024)
# Pars2 = list()
# for(i in 1:4){
#   Pars[names[i]] = foo[i]*(max(unscaled[,i])-min(unscaled[,i])) + min(unscaled[,i])
#   Pars2[names[i]] = test[i]#fit$summary("theta")$median[i]
#   print(foo2[,i]*(max(unscaled[,i])-min(unscaled[,i])) + min(unscaled[,i]))
# }
# State <- c(x = 30, y = 4)
# Time <- seq(0, 20, by = 1)
# out <- as.data.frame(ode(func = LotVmod, y = State, parms = Pars2, times = Time))[,-1]
# matplot(apply(spgasp$mc_yeta,c(1,3),median)+(apply(spgasp$mc_w,1:2,quantile,.5))[-1,],type="l",ylim=c(0,10),add=0,col=c("darkblue","gold"),lty=1,lwd=3)
# matplot(apply(spgasp$mc_yeta,c(1,3),median)+(apply(spgasp$mc_w,1:2,quantile,.9))[-1,],type="l",add=1,col=c("darkblue","gold"),lty=1,lwd=3)
# matplot(apply(spgasp$mc_yeta,c(1,3),median)+(apply(spgasp$mc_w,1:2,quantile,.1))[-1,],type="l",add=1,col=c("darkblue","gold"),lty=1,lwd=3)
# matplot(log(rbind(y)),type="p",pch=19,add=T)
# matplot(apply(spgasp$mc_yeta,c(1,3),mean)[-1,],type="l",add=0,col=c("darkblue","gold"),lty=1,lwd=3,ylim=c(0,10))
# matplot(log(rbind(y)),type="p",pch=19,add=T)

median = c(apply(spgasp$mc_yeta,c(1,3),median)+(apply(spgasp$mc_w,1:2,quantile,.5))[-1,])
median2 = c(apply(spgasp$mc_yeta,c(1,3),median))

MSE(median, c(z[-1,]))
q5 = c(apply(spgasp$mc_yeta,c(1,3),quantile,.05)+(apply(spgasp$mc_w,1:2,quantile,.5))[-1,])
q95 = c(apply(spgasp$mc_yeta,c(1,3),quantile,.95)+(apply(spgasp$mc_w,1:2,quantile,.5))[-1,])

df2 = data.frame(cbind(rep(1:20,2),c(z[-1,]),c(median),c(q5),c(q95),rep(c(2,1),each=20)))
colnames(df2) = c("Measurement","Cases","Median","q5","q95","group")

library(ggplot2)
par(mfrow=c(1,2))
{
  h <- ggplot(data=df2, aes(y=Median, x=Measurement,color=factor(group),fill=factor(group)))+ggtitle("spDLMGP Fit")+
    geom_ribbon(aes(ymin=q5, ymax=q95,color=factor(group),alpha=factor(group)), linetype=0) +
    theme_classic()+coord_cartesian(ylim = c(0, 2*max(median)))+
    geom_line(aes(y=median, color=factor(group)), size=1,alpha=1)+
    scale_color_manual(name="Predictive Median",labels = c("Predator","Prey"),
                       values=c("#e6550d","#3182bd"))+
    scale_fill_manual(name="95 % Credible Interval",labels = c("Predator","Prey"),
                      values=c("#e6550d","#3182bd"))+
    theme(plot.title = element_text(size=30, hjust=.5))+
    guides(colour = guide_legend(override.aes = list(size = 6,alpha=1)))+
    scale_alpha_ordinal(name="95 % Credible Interval",range=c(.25,.4), labels = c("Predator","Prey"))
  
  h1=h+geom_point(data = df2, aes(y=Cases, color=factor(group)),size=3,alpha=1,show.legend = F)+ylab("Log-Species Count")
}

matplot(apply(spgasp$mc_yeta,c(1,3),median)+(apply(spgasp$mc_w,1:2,quantile,.5))[-1,], lty=1,lwd=4, type="l",
        ylim = c(0,5), col=c("darkblue","lightblue"), ylab="Log Species Count", xlab="Measurement")
matplot(apply(spgasp$mc_yeta,c(1,3),median), add=T, type="l",col=c("darkred","darkorange"),lty=2,lwd=4)
legend(5,1,legend=c("Prey (w/ bias term, MSE: .06)", "Pred (w/ bias term, MSE: .06)", "Prey (no bias, MSE: 0.21)", "Pred (no bias, MSE: 0.21)"),
       fill=c("darkblue","lightblue", "darkred","darkorange"), cex=.8)

MSE(median, c(z[-1,]))
MSE(median2, c(z[-1,]))
MSE(median_stan, c(z[-1,]))


median_stan = log(c(fit$summary("y_rep")$median))
q5 = log(c(fit$summary("y_rep")$q5))
q95 = log(c(fit$summary("y_rep")$q95))
df2 = data.frame(cbind(rep(1:20,2),c(z[-1,]),median,q5,q95,rep(c(2,1),each=20)))
colnames(df2) = c("Measurement","Cases","Median","q5","q95","group")
{
  h <- ggplot(data=df2, aes(y=Median, x=Measurement,color=factor(group),fill=factor(group)))+ggtitle("Stan Fit")+
    geom_ribbon(aes(ymin=q5, ymax=q95,color=factor(group),alpha=factor(group)), linetype=0) +
    theme_classic()+coord_cartesian(ylim = c(0, 2*max(median)))+
    geom_line(aes(y=median), size=1,alpha=1)+
    scale_color_manual(name="Predictive Median",labels = c("Predator","Prey"),
                       values=c("#e6550d","#3182bd"))+
    scale_fill_manual(name="95 % Credible Interval",labels = c("Predator","Prey"),
                      values=c("#e6550d","#3182bd"))+
    theme(plot.title = element_text(size=30, hjust=.5))+
    guides(colour = guide_legend(override.aes = list(size = 6,alpha=1)))+
    scale_alpha_ordinal(name="95 % Credible Interval",range=c(.25,.4), labels = c("Predator","Prey"))
  h2=h+geom_point(data = df2, aes(y=Cases, color=factor(group)),size=3,alpha=1,show.legend = F)+ylab("Log-Species Count")
}
library(gridExtra)
grid.arrange(h1,h2,nrow=1)
library(MLmetrics)



# matplot(matrix(log(cbind(out$x,out$y))), pch=1,add=T)#,matrix(log(rbind(y_init,y)),ncol=2))
# matplot(matrix(log(cbind(out$x,out$y)),ncol=2), type="l",add=T)#,matrix(log(rbind(y_init,y)),ncol=2))
# 
# matplot(log(rbind(y_init,y)),type="p",pch=19,add=T)
# matplot(-apply(spgasp$mc_w,1:2,quantile,.5),type="l")
# matplot(-apply(spgasp$mc_w,1:2,quantile,.9),type="l",add=1)
# matplot(-apply(spgasp$mc_w,1:2,quantile,.1),type="l",add=1)
# 

# {
#   nrep=5000
#   emloc = 1
#   thetaloc = 2
#   emulated = matrix(0,nrow=ncol(ysp)-1,ncol=nrep)
#   (test = c(runif(4,.3,.7)))
#   #(test=theta_design[thetaloc,] + runif(2,0,0))
#   plot(theta_design[,1:2])
#   points(rbind(test[1:2]),pch=19,cex=2,col="red")
#   goback = c()
#   for(i in 1:4){
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
#                      F = FFsp[[emloc]],
#                      K = K_draws[,,sample(1:dim(K_draws)[3],1)],
#                      yinit = ysp[thetaloc,1,emloc],
#                      v_draws[,emloc,sample(1:ncol(beta_draws),1)], nugget=0)
#     emulated[,i] = temp$ypred
#   }
#   matplot((emulated[,colSums(emulated<0)<=1])^1,type="l",add=0,col=rgb(0,0,0,alpha=.1))
#   q05 = apply(emulated[,colSums(emulated<0)<=1]^1,1,quantile,.05,na.rm=T)
#   q5 = apply(emulated[,colSums(emulated<0)<=1]^1,1,quantile,.5,na.rm=T)
#   q95 = apply(emulated[,colSums(emulated<0)<=1]^1,1,quantile,.95,na.rm=T)
#   matplot((q05),type="l",col="gold",lwd=2,add=T)
#   matplot((q5),type="l",col="gold",lwd=2,add=1)
#   matplot((q95),type="l",col="gold",lwd=2,add=T)
#   Pars <- c(alpha = goback[1],
#             beta = goback[2],
#             gamma = goback[3],
#             delta = goback[4])
#   out <- as.data.frame(ode(func = LotVmod, y = State, parms = Pars, times = Time))[,-1]
#   matplot(log(out[-1,emloc]),type="l",add=T,col="darkorange",lwd=5)
# }
# 
# 
# 
# 
N <- length(lynx_hare_df$Year) - 1
ts <- 1:N
y_init <- c(lynx_hare_df$Hare[1], lynx_hare_df$Lynx[1])
y <- as.matrix(lynx_hare_df[2:(N + 1), 2:3])
y <- cbind(y[ , 2], y[ , 1]); # hare, lynx order
lynx_hare_data <- list(N = N, ts = ts, y_init = y_init, y = y)

library(cmdstanr)
model = cmdstan_model(stan_file="LV_case_study.stan")
# int<lower = 0> N;           // num measurements
# real ts[N];                 // measurement times > 0
# real y_init[2];             // initial measured population
# real<lower = 0> y[N, 2];    // measured population at measurement times
fit=model$sample(data = lynx_hare_data,
                 chains = 1,
                 iter_warmup = 2500,
                 iter_sampling = 2500, refresh = 10, parallel_chains = 4)
# 
fit$summary("theta")
stan_draws = matrix(fit$draws("theta")[,1,],ncol=4)
apply(stan_draws,2,hist)

for(i in 1:250){
  a_draw = fit$draws("theta")[sample(1:length(test),1,1:4)]
  out = ode(y_init, 1:N, LotVmod, list(alpha = .55, beta = .028, gamma = .8, delta = .025))[,-1]
  matplot((out[,]),type="l",add=0)
}
out = ode(y_init, 1:N, LotVmod, list(alpha = .55, beta = .028, gamma = .8, delta = .025))[,-1]
matplot((out[,]),type="l",add=0)
matplot(rbind(y_init,y[,]),type="l",add=1,lty=2)


ztest = (z)
mu=apply(spgasp$mc_yeta,c(1,3),median)
sd=apply(spgasp$mc_yeta,c(1,3),sd)
(GRS = -sum(((ztest[-1,]-mu)/sd)^2)-2*sum(log(sd)))

stanModel = log(c(fit$summary("y_rep")$mean))
mu=matrix(log(c(fit$summary("y_rep")$mean)),ncol=2)
sd=log(matrix(c(fit$summary("y_rep"))$sd,ncol=2))
(GRS = -sum(((ztest[-1,]-mu)/sd)^2)-2*sum(log(sd)))


# 
