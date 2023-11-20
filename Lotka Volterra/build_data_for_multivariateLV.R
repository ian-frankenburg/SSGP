library(spreadr)
library(ggraph)
library(ggplot2)
library(igraph)
library(gifski)
library(png)
library(grid)
library(gridExtra)

library(lhs)
sims=3
S = 2
A <- maximinLHS(sims, 4) 
theta_design = matrix(0,nrow=sims,ncol=4)
test <- c(alpha = .55, beta = .025, gamma = .8, delta = .024)
theta_design[,1] <- qtruncnorm(A[,1], a=0, b=Inf, mean = 1, sd=.5)
theta_design[,2] <- qtruncnorm(A[,2], a=0, b=Inf, mean = .05, sd=.05)
theta_design[,3] <- qtruncnorm(A[,3], a=0, b=Inf, mean = 1, sd=.5)
theta_design[,4] <- qtruncnorm(A[,4], a=0, b=Inf, mean = .05, sd=.05)

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
Nx = 2; Ny = 1
S = Nx*Ny
t_max = max(times)
FF = array(0, dim=c(sims, 1, length(times)-1, S))
ysp = array(0, dim=c(sims,S,length(times)))
ysp2 = array(0, dim=c(sims,length(times),S))
FF1 = array(0, dim=c(sims, S, length(times)-1))
FFc = array(0, dim=c(sims, 1, length(times)-1, S))
pairs(theta_design)
for(i in 1:sims){
  Pars <- c(alpha = theta_design[i,1], 
            beta = theta_design[i,2], 
            gamma = theta_design[i,3],
            delta = theta_design[i,4])
  out <- as.data.frame(ode(func = LotVmod, y = State, parms = Pars, times = times))[,-1]
  ysp_ = log(out)
  for(sp in 1:S){
    ysp[i,sp,] = as.numeric(ysp_[,sp])
    FF1[i,sp,] = (ysp[i,sp,-length(times)])
    
    
    ysp2[i,,sp] = as.numeric(ysp_[,sp])
    FFc[i,1,,sp] = (ysp2[i,-length(times),sp])
  }
}
for(i in 1:S){
  matplot((t(ysp[,i,])), type="l")
}
tempArr = array(0,dim=c(dim(FF)[1], dim(FF)[2], dim(FF)[3]))
FFsp = list()
for(i in 1:dim(ysp2)[3]){
  for(j in 1:dim(tempArr)[3]){
    tempArr[,,j] = FFc[,1,j,i]
  }
  FFsp[[i]] = tempArr
}
# library(LICORS)
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
library(latex2exp)
set.seed(1)
(test = c(runif(ncol(theta_design),0.1,.9)))
pairs(rbind(theta_design,test),
      cex=c(rep(1, nrow(theta_design)), 2),
      pch=c(rep(1, nrow(theta_design)), 19), 
      col=c(rep("black",nrow(theta_design)),"red"),
      labels = c(TeX('$\\beta$'),TeX('$\\gamma$'), TeX('$\\alpha$')))

goback = c()
for(i in 1:ncol(theta_design)){
  goback[i] = test[i]*(max(unscaled[,i])-min(unscaled[,i])) + min(unscaled[,i])
}
goback
