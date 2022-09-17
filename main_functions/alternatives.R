estimate.K = function(A, K.max=8, cv.rep=1){
  K.rep = rep(0,cv.rep)
  K.ecv = K.rep
  K.ncv = K.rep
  
  for (j in 1:cv.rep) {
    ######## select K ########
    # edge cv
    ecv = ECV.block(A,K.max)
    K.ecv[j] = max(which.min(ecv$l2),which.min(ecv$dev))
    
    # node cv
    ncv = NCV.select(A,K.max)
    K.ncv[j] = max(which.min(ncv$l2),which.min(ncv$dev))
  
  }
  
  K.ecv0 = min(K.ecv)
  K.ncv0 = min(K.ncv)
  
  K.all =as.data.frame(t(c(K.ecv0, K.ncv0)))
  names(K.all) = c("K.ecv", "K.ncv")

  return(K.all)
}

reg.SSP.SBM = function(A, R=8, K.hat.alter, lap=TRUE){
  # detects communities by regularized spherical spectral clustering
  ssc = reg.SSP(A, K=K.hat.alter, lap=lap)
  memb = ssc$cluster
  cluster.matrix.hat = memb.cluster(memb)
  
  # estimates DCSBM parameters given community labels
  dcsbm = DCSBM.estimate(A, memb)
  P.hat.DCSBM = dcsbm$Phat
  # estimates SBM parameters given community labels
  sbm = SBM.estimate(A, memb)
  P.hat.SBM = sbm$Phat

  res.SSP = list(Z.hat=matrix(0,dim(A)[1],R), K.hat=K.hat.alter, memb=memb,
                 P.hat.DCSBM=P.hat.DCSBM, P.hat.SBM=P.hat.SBM,
                 cluster.matrix.hat=cluster.matrix.hat)
  return(res.SSP)
}

LSM = function(A, R=8, K.hat.alter, niter=50){
  # estimates inner product latent space model by projected gradient descent
  # clustering using K-means
  fit <- LSM.PGD(A, k=R, niter=niter)
  Z.lsm = as.matrix(fit$Z)
  kmeans.clust = kmeans(Z.lsm, K.hat.alter)
  memb = kmeans.clust$cluster
  cluster.matrix.hat = memb.cluster(memb)
  res.LSM = list(Z.hat=Z.lsm, K.hat=K.hat.alter, memb=memb,
                 cluster.matrix.hat=cluster.matrix.hat)
  return(res.LSM)
}

LSM.PGD = function(A,k,step.size=0.3,niter=500,trace=0){
  N <- nrow(A)
  ones = rep(1,N)
  M = matrix(1, N, N)
  Jmat <- diag(rep(1,N)) - M/N

  P.tilde <- USVT(A)
  P.tilde[P.tilde>(1-1e-5)] <- (1-1e-5)
  P.tilde[P.tilde< 1e-5] <- 1e-5

  Theta.tilde <- logit(P.tilde)

  alpha_0 <- solve(N*diag(rep(1,N))+M,rowSums(Theta.tilde))

  G <- Jmat%*%(Theta.tilde - outer(alpha_0,alpha_0,"+"))%*%Jmat

  eig <- eigs_sym(A=G,k = k)
  eig$values[eig$values<=0] <- 0

  Z_0 <- t(t(eig$vectors[,1:k])*sqrt(eig$values[1:k]))
  obj <- NULL
  step.size.z <- step.size/norm(Z_0,"2")^2
  step.size.alpha <- step.size/(2*N)
  for(i in 1:niter){
    Theta.hat <- alpha_0 %*% t(rep(1,N)) + rep(1, N) %*% t(alpha_0) + Z_0 %*% t(Z_0)
    Phat <- sigmoid(Theta.hat)
    tmp.obj <- (sum(A*log(Phat)) + sum((1-A)*log(1-Phat)) - sum(diag(log(1-Phat))))/2
    if(trace>0){
      print(tmp.obj)
    }
    obj <- c(obj,tmp.obj)
    Z <- Z_0 + 2*step.size.z*(A-Phat)%*%Z_0
    alpha <- alpha_0 + 2*step.size.alpha*(A-Phat)%*%matrix(rep(1,N))
    Z <- Jmat%*%Z

    Z_0 <- Z
    alpha_0 <- alpha
  }

  Theta.hat <- alpha_0 %*% t(rep(1,N)) + rep(1, N) %*% t(alpha_0) + Z_0 %*% t(Z_0)
  Phat <- sigmoid(Theta.hat)
  tmp.obj <- (sum(A*log(Phat)) + sum((1-A)*log(1-Phat)) - sum(diag(log(1-Phat))))/2
  obj <- c(obj,tmp.obj)
  return(list(Z=Z,alpha=alpha,Phat=Phat,obj=obj))

}

sigmoid = function(x, a = 1, b = 0){
  if (length(x) == 0)
    return(c())
  stopifnot(is.numeric(a), is.numeric(b))
  a <- a[1]
  b <- b[1]
  return(1/(1 + exp(-a * (x - b))))
}

alternatives = function(A, R=8, K.max=8, cv.rep=1, niter=50, lap=TRUE){
  K.all = estimate.K(A, K.max=K.max, cv.rep=cv.rep)
  # ECV 
  K.hat.alter = K.all$K.ecv
  res.LSM = LSM(A, R=R, K.hat.alter=K.hat.alter, niter=niter)
  res.SSP = reg.SSP.SBM(A, R=R, K.hat.alter=K.hat.alter, lap=lap)
  res.ECV = list(res.LSM=res.LSM, res.SSP=res.SSP)
  # NCV 
  K.hat.alter = K.all$K.ncv
  res.LSM = LSM(A, R=R, K.hat.alter=K.hat.alter, niter=niter)
  res.SSP = reg.SSP.SBM(A, R=R, K.hat.alter=K.hat.alter, lap=lap)
  res.NCV = list(res.LSM=res.LSM, res.SSP=res.SSP)
  
  res = list(res.ECV=res.ECV, res.NCV=res.NCV)
  return(res)
}


