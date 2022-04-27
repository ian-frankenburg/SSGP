set.seed(1)
source("helper.R")
set.seed(1)
spSSGP = spDLMGP(z=zsp[,-1], y=yin[,-1],
                       yy = array(ysp2[,-1,],
                            dim=c(nrow(ysp2),ncol(ysp2)-1,dim(ysp2)[3])),
                       FF=FF, FFc=FFsp,
                       yinit = ysp[1,,1],
                       delta_v=.95, delta_w=.95,
                       params=theta_design,
                       tune=rep(.001,3), tune2=rep(.01,3),
                       niter=1000, burnin=1000/2,
                       nugget=0, S=S)
{
  library(latex2exp)
  par(mfrow=c(2,2))
  for(i in 1:length(test)){
    drw = matrix(rbeta(length(spSSGP$calibrate),1,1), ncol=ncol(spSSGP$calibrate))
    if(i==1){
      hist(spSSGP$calibrate[,i], breaks=10, main="", ylab="",
           col=rgb(0/255,0/255,.502,.7), xlim=c(0,1), cex.lab=2,
           freq = 0, xlab=TeX(
             paste0('$\\eta_', i, "$")), add=0)
      hist(drw[,i], breaks=10, freq = 0, xlab=TeX(
        paste0('$\\eta_', i, "$")),
        col=rgb(.678, .847, .902, alpha=.2),add=T)
      legend(.5, 5, c("Prior", "Posterior"),
             fill=c(rgb(.678, .847, .902),
                    rgb(0/255,0/255,.502, .7),
                    rgb(1, 0, 0, alpha=.7)),
             bty = 'n',text.font=2,cex=.7)
    }else{
      hist(spSSGP$calibrate[,i], breaks=10, main="",ylab="",
           col=rgb(0/255, 0/255, .502, .7), xlim=c(0,1), cex.lab=2,
           freq = 0, xlab=TeX(
             paste0('$\\eta_', i, "$")), add=0)
      hist(drw[,i], xlim=c(0,1), breaks=10, freq = 0, xlab=TeX(
        paste0('$\\eta_', i, "$")), main="",
        col=rgb(.678, .847, .902, alpha=.2),add=T)
    }
    points(y=0,x=test[i], col="red", lwd=5,pch=19)
  }
  hist(sqrt(spSSGP$mc_sigma2), cex.lab=2,ylab="", xlim=c(0,5),
       col=rgb(0/255,0/255,.502,.7), breaks=10, main="",
       freq = 0, xlab=TeX('$\\sigma_z $'), add=0)
  hist((rgamma(nrow(spSSGP$calibrate),1,1)), breaks=10, freq = 0, xlab=TeX(
    paste0('$\\eta_', i, "$")),
    col=rgb(.678, .847, .902, alpha=.5),add=T)
  points(y=0,x=sigma_z, col="red", lwd=4,pch=19)
  # plot(c(0,0),col="white",xlab="",ylab="",xaxt="n",yaxt="n",frame.plot = FALSE)
}


