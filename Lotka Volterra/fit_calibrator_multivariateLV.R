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
set.seed(1)
fullySpatial = spDLMGP(z=t(z)[,-1], y=yin[,-1], 
                       y_cube = array(ysp2[,-1,], dim=c(nrow(ysp2),ncol(ysp2)-1,dim(ysp2)[3])),
                       yinit = ysp2[1,1,],
                       niter=500, burnin=500/2,
                       FF=FF, FF_cube=FFsp,
                       C0=diag(2), W=diag(2),
                       params=theta_design,
                       delta_v=.9, delta_w=.95,
                       tune=rep(.001,2), tune2=rep(.1,2), 
                       nugget=1e-5)
