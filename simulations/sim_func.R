#####################################################################################
# Codes for simulation studies: 
# The functions for generating simulated data & evaluating performances,
# which are used to support numerical simulation studies.
#####################################################################################

#####################################################################################
# Functions for generation of simulated data
#####################################################################################
generate.data = function(n, K.true=2, r=5, rsn=0.2, Mu.true, unbalan=F, ka=0.5){
  if(!unbalan){
    Z.true = NULL
    B.true = NULL
    for (k in 1:K.true) {
      mu.k = Mu.true[k,]
      sigma = rsn^2
      Z.true = rbind(Z.true, mvrnorm(n/K.true, mu.k, diag(sigma,r)))
      B.true = rbind(B.true, t(matrix(rep(mu.k,n/K.true),r,n/K.true)))
    }
    cluster.matrix.true = Matrix(1,n/K.true,n/K.true,sparse = T)
    for (k in 2:K.true) {
      cluster.matrix.true = bdiag(cluster.matrix.true,matrix(1,n/K.true,n/K.true))
    }
  } else {
    if(K.true!=2){print("Warning: the true K is not 2")}
    nn = n
    n.big = floor(nn/(1+ka))
    n.small1 = nn - n.big
    n.vec = c(n.small1, n.big)
    Z.true = NULL
    B.true = NULL
    for (k in 1:length(n.vec)) {
      mu.k = Mu.true[k,]
      sigma = rsn^2
      Z.true = rbind(Z.true, mvrnorm(n.vec[k], mu.k, diag(sigma,r)))
      B.true = rbind(B.true, t(matrix(rep(mu.k,n.vec[k]),r,n.vec[k])))
    }
    cluster.matrix.true = Matrix(1,n.vec[1],n.vec[1],sparse = T)
    for (k in 2:K.true) {
      cluster.matrix.true = bdiag(cluster.matrix.true,matrix(1,n.vec[k],n.vec[k]))
    }
  }
  
  Theta.true = Z.true %*% t(Z.true)
  P.true = exp(Theta.true)/(1+exp(Theta.true))
  A = Matrix(rbinom(length(P.true),1,P.true),n,n,sparse = T)
  A[lower.tri(A)] <- t(A)[lower.tri(A)]
  diag(A) = 0
  
  
  return(list(A=A, Z.true=Z.true, B.true=B.true, Mu.true=Mu.true,
         cluster.matrix.true=cluster.matrix.true,
         Theta.true=Theta.true, P.true=P.true))
}


################################################################################
# Functions for evaluating performances of proposed methods
################################################################################
evaluation = function(Z.hat, Z.true, cluster.matrix.hat, cluster.matrix.true, 
                      P.true, Theta.true, K.hat=4, K.true=4, SBM=F, P.hat.SBM){
  if(length(Z.hat) > 0){
    n = dim(Z.hat)[1]
    p = dim(Z.hat)[2]
  } else {
    n = dim(cluster.matrix.true)[1]
  }
  
  diff.clu.mat = as.matrix(cluster.matrix.hat - cluster.matrix.true)
  group_error_rate = sum(abs(diff.clu.mat[upper.tri(diff.clu.mat)]))/(n*(n-1)/2)
  
  if(!SBM){
    Theta.hat = Z.hat %*% t(Z.hat)
    P.hat = exp(Theta.hat)/(1+exp(Theta.hat))
    prob.error = sum((P.hat-P.true)^2) / sum((P.true)^2)
    Theta.error = sum((Theta.hat-Theta.true)^2) / sum((Theta.true)^2)
  } else {
    P.hat = P.hat.SBM
    prob.error = sum((P.hat-P.true)^2) / sum((P.true)^2)
    Theta.error = 1
  }
  res = c(as.numeric(K.hat==K.true), K.hat, group_error_rate, prob.error, Theta.error)
  res = as.data.frame(t(res))
  names(res) = c("prop", "K.hat", "cluster.error", "prob.error", "Theta.error")
  
  return(res)
}

evaluation.alter = function(res.alter){
  ECV.LSM =  evaluation(res.alter$res.ECV$res.LSM$Z.hat, Z.true, 
             res.alter$res.ECV$res.LSM$cluster.matrix.hat, cluster.matrix.true, 
             P.true, Theta.true, res.alter$res.ECV$res.LSM$K.hat, K.true)
  
  NCV.LSM =  evaluation(res.alter$res.NCV$res.LSM$Z.hat, Z.true, 
             res.alter$res.NCV$res.LSM$cluster.matrix.hat, cluster.matrix.true, 
             P.true, Theta.true, res.alter$res.NCV$res.LSM$K.hat, K.true)
  
  ECV.SSP =  evaluation(res.alter$res.ECV$res.SSP$Z.hat, Z.true, 
             res.alter$res.ECV$res.SSP$cluster.matrix.hat, cluster.matrix.true, 
             P.true, Theta.true, res.alter$res.ECV$res.SSP$K.hat, K.true, 
             SBM=T, res.alter$res.ECV$res.SSP$P.hat.DCSBM)
  
  NCV.SSP =  evaluation(res.alter$res.NCV$res.SSP$Z.hat, Z.true, 
             res.alter$res.NCV$res.SSP$cluster.matrix.hat, cluster.matrix.true, 
             P.true, Theta.true, res.alter$res.NCV$res.SSP$K.hat, K.true, 
             SBM=T, res.alter$res.NCV$res.SSP$P.hat.DCSBM)
  

  index.alter = list(ECV.LSM=ECV.LSM, ECV.SSP=ECV.SSP, NCV.LSM=NCV.LSM,
                     NCV.SSP=NCV.SSP)

  return(index.alter)
}

evaluation.pro = function(res, Z.true, cluster.matrix.true, 
                          P.true, Theta.true, K.true){
  n_lam = res$Opt_num
  Z.hat = res$Z.list[[n_lam]]
  cluster.matrix.hat = res$cluster.matrix.list[[n_lam]]
  K.hat = res$K.vec[n_lam]
  index = evaluation(Z.hat, Z.true, cluster.matrix.hat, cluster.matrix.true, P.true, Theta.true, K.hat, K.true)
  return(index)
}




