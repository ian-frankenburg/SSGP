D = diag(degree(network))
D_sqrt_inv = diag(1/sqrt(diag(D))) 
A = as.matrix(mat)
rho = 0.99
W = diag(nrow(A)) + rho * D_sqrt_inv %*% A %*% D_sqrt_inv
heatmap(W,Rowv = NA,Colv = NA)
matrixcalc::is.positive.definite(W)