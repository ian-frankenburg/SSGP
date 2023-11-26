library(spreadr)
library(ggraph)
library(ggplot2)
library(igraph)
library(gifski)
cut.at.n.tile <- function(X , n = 4){ 
  cut( X , breaks = quantile( X , 
                              probs = (0:n)/n , na.rm = TRUE ) ,
       labels = 1:n,
       include.lowest = TRUE )}
data("pnet")  # load and inspect the igraph object
plot(pnet)
start_run <- data.frame(
  node=c("beach", "speck"),
  activation=rep(1,34))

# # make an adjacency matrix and randomly fill some cells with 1s
N = 11
network_orig=network <- make_tree(N, children = 2,mode = "undirected")
plot(network)
l = layout_nicely(network)
mat = as_adjacency_matrix(network)
act = sample(c(1:10),size = N,replace = T)
act = rep(0,N)
act[1] = 10
initial_df <- data.frame(node=c(1,2:(N)), activation=act, stringsAsFactors=FALSE)
a = .01; b = .05; c = 0
results <- spreadr(mat, initial_df, include_t0 = T,
                   decay = a, retention = b, suppress = c)
results$activation=(results$activation+0)
dat = (reshape(results, idvar = "node", timevar = "time", direction = "wide")[,-1])
ggplot(results, aes(x=time, y=activation, color=node)) +
  geom_point() +
  geom_line()

dat = reshape(results, idvar = "node", timevar = "time", direction = "wide")[,-1]
(foo=cut.at.n.tile(as.numeric(unlist(dat)), n=2))
levels(foo) = 1:length(levels(foo))
foo = unlist(dat)
foo2 = matrix(foo, nrow(dat), ncol(dat), byrow = F)
pal = (heat.colors(length(levels(foo))))
# display.brewer.pal(n = length(levels(foo)), name = 'YlGnBu')
# pal = brewer.pal(n = length(levels(foo)), name = 'YlOrBr')
deg = deg2 = degree(network, mode="all")
mat=as_adjacency_matrix(network)
baz = rank(unlist(dat))
bazMat = matrix(baz, nrow(dat), ncol(dat), byrow = F)
# a = c(seq(1,0.1,length.out=45),rep(.1,5))
# for(i in seq(1,ncol(dat),length.out=9)){
#   (temp = exp(dat[, i]))
#   temp = bazMat[,i]
#   network = network_orig
#   V(network)$label = NA
#   V(network)$color <- as.numeric(foo2[,i])
#   V(network)$size <- 3*degree(network)
# 
#   deg2[deg>10] = 10*15
#   deg2[!deg>1]=10*8
#   if(i<10){
#     png(paste0("/Users/ianfrankenburg/Desktop/DyTSA/Paper 2/Networks/figs/0figs0",i,".png"),
#         width = 1200,height = 1200)
#     #print(paste0("0figs0",i,".png"))
#   }else{
#     #print(paste0("0figs",i,".png"))
#     png(paste0("/Users/ianfrankenburg/Desktop/DyTSA/Paper 2/Networks/figs/0figs",i,".png"),
#         width = 1200,height = 1200)
#   }
#   # plot.igraph(x=network, edge.arrow.size=0.2,
#   #             vertex.color = adjustcolor(pal[as.numeric(foo2[,i])], alpha.f = .7),
#   #             layout = l, add = 0, vertex.size=1/deg2*20, main=paste0("Time ", i))
#   #
#   ramp = colorRampPalette(c("lightblue","darkred"))
#   foo = factor(foo2[,i],ordered = 1,nmax = length(unique(foo2[,i])))
#   levels(foo) = (ramp(length(levels(foo))))
#   set.seed(1)
#   print(ggraph(network, layout = "lgl") +
#           geom_edge_link0(edge_colour = "black",
#                           edge_width = 0.7, edge_alpha = .8) +
#           geom_node_point(aes(fill = foo, size=(scale(foo2[,i]))),
#                           shape = 21, alpha=a[i]) +
#           scale_fill_manual(values=ramp(length(levels(foo))))+
#           ggtitle(paste("time", round(i)))+
#           theme(legend.position = "none",
#                 panel.background = element_rect(fill = "transparent"), # bg of the panel
#                 plot.background = element_rect(fill = "transparent"), # bg of the plot
#                 panel.grid.major = element_blank(),
#                 panel.grid.minor = element_blank())+
#           labs(axis.title.x=element_blank(),
#                axis.text.x=element_blank(),
#                axis.ticks.x=element_blank())+
#   scale_size(range = c(2, 30)))
#   deg2=deg
#   dev.off()
# }
# png_files <- list.files("/Users/ianfrankenburg/Desktop/DyTSA/Paper 2/Networks/figs", pattern = "*.png", full.names = TRUE)
# gifski(png_files, gif_file = "/Users/ianfrankenburg/Desktop/DyTSA/Paper 2/Networks/network.gif",
#        width = 8000, height = 8000, delay = 1)
# unlink("/Users/ianfrankenburg/Desktop/DyTSA/Paper 2/Networks/figs/*")

library(png)
library(grid)
library(gridExtra)

library(lhs)
sims=50
cube = maximinLHS(sims,k = 2)
dec = qunif(cube[,1],0,.2)
ret = qunif(cube[,2],0,.5)
S = nrow(dat)
t_max = 11
times = (0:t_max)
FF = array(0, dim=c(sims, 1, t_max, S))
ysp = array(0, dim=c(sims,S,length(times)))
ysp2 = array(0, dim=c(sims,length(times),S))
FF1 = array(0, dim=c(sims, S, length(times)-1))
FFc = array(0, dim=c(sims, 1, length(times)-1, S))
theta_design = cbind(dec,ret)
pairs(theta_design)
for(i in 1:sims){
  a = dec[i]
  b = ret[i]
  c = 0
  results <- spreadr(mat, initial_df, include_t0 = T,
                     decay = a, retention = b, suppress = c, time = t_max)
  results$activation=(results$activation+0)
  dat = reshape(results, idvar = "node", timevar = "time", direction = "wide")
  ysp_ = log(dat[,-1]+1)
  for(sp in 1:S){
    ysp[i,sp,] = as.numeric(ysp_[sp,])
    FF1[i,sp,] = (ysp[i,sp,-length(times)])
    
    
    ysp2[i,,sp] = as.numeric(ysp_[sp,])
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
# 
RcppParallel::setThreadOptions(numThreads = 4)
a = goback[1]
b = goback[2]
c = 0
library(ggplot2)
results <- spreadr(mat, initial_df,include_t0 = T,
                   decay = a, retention = b, suppress = c, time = t_max)
results$activation=(results$activation+0)
ggplot(results, aes(x=time, y=activation, color=node)) +
  geom_point() +
  geom_line()

z = (t(reshape(results, idvar = "node", timevar = "time", direction = "wide")[,-1]))
sigma_z = .25
z = log(z+1)
z = matrix(rnorm(nrow(z)*ncol(z),0,sigma_z),nrow(z),ncol(z)) + z
matplot(z,type="l")

