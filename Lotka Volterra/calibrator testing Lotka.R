set.seed(1)
library(MLmetrics)
library(Matrix)
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
sims=25
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
times = Time
p1=matrix(0,nrow(theta_design), length(times))
p2=matrix(0,nrow(theta_design), length(times))
Nx = 2; Ny = 1
S = Nx*Ny
ysp = array(0, dim=c(sims,Nx*Ny,length(times)))
ysp2 = array(0, dim=c(sims,length(times),Nx*Ny))
FF1 = array(0, dim=c(sims, Nx*Ny, length(times)-1))
FFc = array(0, dim=c(sims, 1, length(times)-1, Nx*Ny))
for(i in 1:sims){
  Pars <- c(alpha = theta_design[i,1], 
            beta = theta_design[i,2], 
            gamma = theta_design[i,3],
            delta = theta_design[i,4])
  out <- as.data.frame(ode(func = LotVmod, y = State, parms = Pars, times = Time))[,-1]+1
  ysp_ = log(out)
  for(sp in 1:(Nx*Ny)){
    ysp[i,sp,] = (ysp_[,sp])
    FF1[i,sp,] = (ysp[i,sp,-length(times)])
    
    ysp2[i,,sp] = (ysp_[,sp])
    FFc[i,1,,sp] = (ysp2[i,-length(times),sp])
  }
}
FF = array(0, dim=c(Nx*Ny*sims, Nx*Ny, length(times)-1))
for(i in 1:dim(FF)[3]){
  l <- split(FF1[,,i], rep(1:ncol(FF1[,,i]), each = nrow(FF1[,,i])))
  FF[,,i] = as.matrix(Matrix::bdiag(l))
}
tempArr = array(0,dim=c(dim(FFc)[1], dim(FFc)[2], dim(FFc)[3]))
FFsp = list()
for(i in 1:dim(ysp2)[3]){
  for(j in 1:dim(tempArr)[3]){
    tempArr[,,j] = FFc[,1,j,i]
  }
  FFsp[[i]] = tempArr
}

l=seq(1,length(Time),length.out=max(Time)+1)
for(s in 1:dim(ysp)[2]){
  matplot(t(ysp[,s,]), type="l", lwd=2, lty=1, main="Experimental Design", 
          ylab="log-population size", xlab="Year", xaxt="n")
  axis(1, at = seq(1,length(Time),length.out=max(Time)), 
       labels = seq(1,20),
       las=2)
  points(round(l),log(yfull[,s]), pch=19, cex=2)
}

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
# yinit = unique(ysp[,,1])
prtb =runif(4,0,0.001)
test = theta_design[8,]
Pars <- c(alpha = unscaled[8,1]+prtb[1],
          beta = unscaled[8,2]+prtb[2],
          gamma = unscaled[8,3]+prtb[3],
          delta = unscaled[8,4]+prtb[4])
out <- as.data.frame(ode(func = LotVmod, y = State, parms = Pars, times = Time))[,-1]
z = t(log(out))
matplot(t(z),type="l")
# z = log(t(yfull))
for(s in 1:dim(ysp2)[3]){
  matplot(t(ysp[,s,]),type="l",add=0)
  points((z[s,]), pch=19)
}
st_time = Sys.time();
p2=1
yin = matrix(0,nrow=dim(FF)[1], ncol=dim(ysp)[3])
for(i in 1:dim(ysp)[3]){
  yin[, i] = c(ysp[,,i])
}
p2=2

a=FF[,,1]%*%matrix(c(.9,.1,.1,.9),2,2,byrow = 1)%*%t(FF[,,1])

# const mat& z, const mat& y,  const cube& y_cube, const vec yinit,
# const uword& niter, const uword& burnin,
# const cube& FF, const field<cube>& FF_cube, 
# const mat& C0, const mat& W,
# const mat& params,
# const double delta_v, const double delta_w,
# vec& tune, const vec& tune2, 
# const double nugget
gasp = spSSGP(z=z[,-1], y=yin[,-1], y_cube = array(ysp2[,-1,],dim=c(nrow(ysp2),ncol(ysp2)-1,dim(ysp2)[3])), 
              yinit = unique(yin[,1]), niter=5000, burnin=5000/2,
              FF=FF, FF_cube = FFsp,
              C0 = 1*diag(p2), W=diag(p2),
              params=theta_design,
              delta_v=.95, delta_w=.95,
              tune=rep(.01,4), tune2=rep(.01,4),
              nugget=0)

# em = emulate(y=yin[,-1],pred_params=theta_test,params=theta_design,
#         K = apply(gasp$mc_K,1:2,median), F = FF, yinit = yinit,
#         beta = rowMeans(gasp$mc_beta), theta = apply(gasp$mc_theta,1:2,mean),
#         v = rowMeans(gasp$v),nugget = .1,S = 2)
plot(rowMeans(gasp$yeta[,,1]), type="l")
points(ysp[8,1,-1])

plot(apply(gasp$v,1,mean))
matplot(gasp$v,type="l")

# View(apply(gasp$mc_K,1:2,median))
matplot(t(gasp$mc_beta),type = "l")
plot(t(gasp$mc_sigma2))
matplot(gasp$calibrate,type="l")
colMeans(gasp$calibrate)
test
matplot(t(apply(gasp$mc_theta,1:2,median)), type="l",ylim=c(-5,5))
for(i in 1:250){
  matplot(t(gasp$mc_theta[,,i]), type="l", add=T)
}
# matplot((apply(spgasp$mc_w,1:2,median)),type="l",ylim=c(0,1))
# matplot((apply(spgasp$mc_w,1:2,quantile,.95)),type="l",add=T)
# matplot((apply(spgasp$mc_w,1:2,quantile,.5)),type="l",add=T)


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
# fit=model$sample(data = lynx_hare_data,
#                  chains = 1,
#                  iter_warmup = 500,
#                  iter_sampling = 500, refresh = 10, parallel_chains = 4)
# # 
stan_draws = matrix(fit$draws("theta")[,1,],ncol=4)
# apply(stan_draws,2,hist)

# for(i in 1:250){
#   a_draw = fit$draws("theta")[sample(1:length(test),1,1:4)]
#   out = ode(y_init, 1:N, LotVmod, list(alpha = .55, beta = .028, gamma = .8, delta = .025))[,-1]
#   matplot((out[,]),type="l",add=0)
# }
# out = ode(y_init, 1:N, LotVmod, list(alpha = .55, beta = .028, gamma = .8, delta = .025))[,-1]
# matplot((out[,]),type="l",add=0)
# matplot(rbind(y_init,y[,]),type="l",add=1,lty=2)
# 


library(latex2exp)
tempTest = matrix(0, nrow = length(stan_draws[,1]), ncol=4)
foo = unscaled
for(i in 1:4){
  tempTest[,i] = (stan_draws[,i]-min(foo[,i]))/(max(foo[,i])-min(foo[,i]))
}
par(mfrow=c(2,2))
for(i in 1:length(test)){
  drw = (rbeta(length(gasp$calibrate),1,1))
  if(i==1){
    drw = (rbeta(length(gasp$calibrate[,i])*3,1,1))
    hist(tempTest[,i], xlim=c(0,1), breaks=10, probability = 1, xlab=TeX(
      paste0('$\\eta_', i, "$")), main="",
      col=rgb(1, 0, 0, alpha=.5),add=0)
    hist(gasp$calibrate[,i], xlim=c(0,1), breaks=10, 
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
    drw = (rbeta(length(gasp$calibrate[,i]),1,1))
    hist(tempTest[,i], xlim=c(0,1), breaks=10, probability = 1, xlab=TeX(
      paste0('$\\eta_', i, "$")), main="",
      col=rgb(1, 0, 0, alpha=.5),add=0)
    hist(gasp$calibrate[,i], xlim=c(0,1), breaks=10, 
         col=rgb(0/255, 0/255, .502, .7),
         probability = 1, xlab=TeX(
           paste0('$\\eta_', i, "$")), main="", add=1)
    hist(drw, xlim=c(0,1), breaks=10, probability = 1, xlab=TeX(
      paste0('$\\eta_', i, "$")), main="",
      col=rgb(.678, .847, .902, alpha=.4),add=T)
  }
}

