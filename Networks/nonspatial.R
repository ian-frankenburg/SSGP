FF = array(0, dim=c(S*sims, 1, length(times)-1))
for(i in 1:dim(FF)[3]){
  l <- split(FF1[,,i], rep(1:ncol(FF1[,,i]), each=nrow(FF1[,,i])))
  FF[,,i] = as.matrix(unlist(l))
}
yin = matrix(0,nrow=dim(FF)[1], ncol=dim(ysp)[3])
for(i in 1:dim(ysp)[3]){
  yin[, i] = c(ysp[,,i])
}
st_time = Sys.time();
nonspatial = DLMGP(z=z[-1,], y=yin[,-1],
               yy = array(ysp2[,-1,],
                          dim=c(nrow(ysp2),ncol(ysp2)-1,dim(ysp2)[3])),
               FF=FF, FFc=FFsp,
               m0=matrix(rep(0,1)), 
               C0=matrix(1), W=matrix(1),
               yinit = ysp[1,,1], alpha0=.1, beta0=.1,
               delta_v=.95, delta_w=.95,
               params=theta_design,
               pred_params = test,
               tune=rep(.001,2), tune2=rep(.001,2), 
               tune3=.25, tune4=rep(.01,2),
               niter=500, burnin=500/2,
               nugget=0, sigma_z=rep(1,1),
               a0z = 1, b0z=1, S=1, S2=ncol(z))

(Sys.time()-st_time)

