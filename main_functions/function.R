##################################################################################
# This document includes main functions of proposed methods, 
# which are used to support numerical simulation studies and real data analysis
# in the paper
# "Consistent Estimation of the Number of Communities via Regularized Network Embedding"
##################################################################################

############################# Functions for main algorithms ############################
fusion.network = function(A, Omega.A, sample.index.A, gamma.index, lambda1, 
                          lambda2, lambda3, Z.int, B.int, a=3, kappa=1, alpha=1,
                          eps = 5e-2, niter = 20, niter.Z=5, update.B="ADMM", 
                          local.oppro=F, ad.BIC=F, merge.all=T){
  n = as.integer(dim(A)[1])
  p = dim(Z.int)[2]
  nc = n*(n-1)
  K.int = length(unique(apply(B.int*B.int,1,sum)))
    
  gamma.int = Matrix(0,nrow=dim(Omega.A)[2],ncol=p,sparse = T)
  V.int = t(t(B.int) %*% Omega.A)
  Z.new = Z.int
  B.new = B.int
  gamma.new = gamma.int
  V.new = V.int
  
  iter = 0
  K.hat = K.int
  niter.K = 3
  iter.K = 0
  while( K.hat == K.int && iter.K<niter.K ){
    iter.K = iter.K+1
    eps.diff.Z = 10
    eps.diff.V = 10

    while( eps.diff.Z+eps.diff.V >= eps && iter<niter )
    { 

      iter = iter+1
      Z.old = Z.new
      B.old = B.new
      gamma.old = gamma.new
      V.old = V.new
      
      Theta.old = Z.old %*% t(Z.old)
      P.old = exp(Theta.old)/(1+exp(Theta.old))
      P.old[Theta.old > 700] = 1
      diag(P.old) = 0
      
      # update Z
      if(!local.oppro){
        alphaZ = 2*alpha
        obj.value.old = 1
        obj.value.new = 10
        iter.Z=0
        while( obj.value.new > obj.value.old  && iter.Z<niter.Z ){
          iter.Z = iter.Z+1
          alphaZ = alphaZ/2
          g.L = - 2/n * t(A - P.old) %*% Z.old
          g.J1 = 2*lambda1 * (Z.old - B.old)
          Z.m = sqrt(apply(Z.old*Z.old, 2, sum))
          g.J3 = t( t(Z.old) * c(mcp_d(Z.m,lambda3) / (Z.m+10^(-10))) )
          Z.new = Z.old - alphaZ*(g.L+g.J1+g.J3)
          obj.value.old = obj.value(A, Z.old, B.old, lambda1, lambda2, lambda3, Omega.A)
          obj.value.new = obj.value(A, Z.new, B.new, lambda1, lambda2, lambda3, Omega.A)
        }
      } else {
        Z.m = apply(Z.old, 2, function(x) sqrt(sum(x^2)))
        c.Z1 = (lambda1 + mcp_d(Z.m,lambda3) / (Z.m+10^(-10)) / 2)^(-1)
        Z.new = t(t( lambda1 * B.old + 1/nc * t(A - P.old) %*% Z.old ) * c.Z1)
      }
      
      obj.value.old.B = obj.value(A, Z.new, B.old, lambda1, lambda2, lambda3, Omega.A)
      # update B
      if(update.B=="AMA"){
        gamma.up = Matrix(t(t(gamma.old) %*% gamma.index),sparse = T)
        B.new = Z.new + gamma.up
      } else if(update.B=="ADMM"){
        kappa.lam = kappa/2/(lambda1+10^(-10))
        gamma.up = Matrix(t(t(gamma.old + kappa.lam*V.old) %*% gamma.index),sparse = T)
        Z.new.mean = n*kappa.lam*apply(Z.new, 2, mean)/(1+n*kappa.lam)
        Y.new = 1/(1+n*kappa.lam)*(Z.new + gamma.up)
        B.new = t(t(Y.new) + Z.new.mean)
      } else {print("The method of updating B must be AMA or ADMM.")}
  
      # update V
      V.old = t(t(B.new) %*% Omega.A)
      delta = V.old - gamma.old/kappa
      V.new = Matrix(t(apply(delta, 1, MCP_soft, lambda=lambda2, kappa=kappa)),sparse = T)
      # update gamma
      gamma.new = gamma.old + kappa*(V.new - V.old)
      
      eps.diff.Z = sqrt(sum((Z.new - Z.old)^2)/length(Z.old))
      eps.diff.V = sqrt(sum((V.new - V.old)^2)/length(V.old))
      
      obj.value.new.B = obj.value(A, Z.new, B.new, lambda1, lambda2, lambda3, Omega.A)

  }
    
    Z.hat = Z.new
    B.hat = B.new
    V.hat = V.new
    
    V.norm = apply(V.hat*V.hat,1,sum)
    diff.Bij = rep(1,nc/2)
    diff.Bij[sample.index.A[3,V.norm == 0]] = 0
    cluster.results = get.cluster(diff.Bij, n, sample.index.n, merge.all=merge.all)
    K.hat = cluster.results$K.hat
    BIC.vec = BIC.value(A, Z.hat, K.hat, adjust=ad.BIC, B.hat)
  }
  
  Z.hat = Z.new
  B.hat = B.new
  V.hat = V.new
  
  V.norm = apply(V.hat*V.hat,1,sum)
  diff.Bij = rep(1,nc/2)
  diff.Bij[sample.index.A[3,V.norm == 0]] = 0
  cluster.results = get.cluster(diff.Bij, n, sample.index.n, merge.all=merge.all)
  K.hat = cluster.results$K.hat
  BIC.vec = BIC.value(A, Z.hat, K.hat, adjust=ad.BIC, B.hat)
  
  return(list(Z.hat=Z.hat, B.hat=B.hat, V.hat=V.hat, K.hat=K.hat,
              cluster.results=cluster.results, BIC.vec=BIC.vec))
 }


obj.value = function(A, ZZ, BB, lambda1, lambda2, lambda3, Omega.A){
  n = dim(ZZ)[1]
  p = dim(ZZ)[2]
  
  Theta.hat = ZZ %*% t(ZZ)
  likelihood.mat =  log(1 + exp(Theta.hat)) - Theta.hat * A
  likelihood.mat[Theta.hat > 700] = Theta.hat[Theta.hat > 700] - Theta.hat[Theta.hat > 700] * as.matrix(A)[as.matrix(Theta.hat) > 700]
  L = 2*sum(likelihood.mat[upper.tri(likelihood.mat)])/n/(n-1)
  
  J1 = lambda1*sum((ZZ - BB)^2)
  
  BB.diff = t(t(BB) %*% Omega.A)
  BB.diff.norm = Matrix(sqrt(apply(BB.diff*BB.diff, 1, sum)), sparse = T)
  J2 = sum(MCP_func(BB.diff.norm, lambda2))
  
  Zm.norm = Matrix(sqrt(apply(ZZ*ZZ, 1, sum)), sparse = T)
  J3 = sum(MCP_func(Zm.norm, lambda3))
  
  obj.val = L + J1 + J2 + J3
  return(obj.val)
}


############################## Functions for tuning lambdas ############################
network.comm.num = function(A, sample.index.n, lambda, Z.int, B.int, 
                                 a=3, kappa=1, alpha=1, eps=5e-2, niter=20, 
                                 niter.Z=5, update.B="ADMM",local.oppro=F, merge.all=T,
                                 ad.BIC=F, Fully.Connected=T, trace=F, line.search=T){
  n  = as.integer(dim(A)[1])
  p  = dim(Z.int)[2]
  
  Omega.list     = gen.weight(A, sample.index.n, Fully.Connected=Fully.Connected)
  Omega.A        = Omega.list$Omega.A
  sample.index.A = Omega.list$sample.index.A
  gamma.index    = Omega.list$gamma.index
  
  if(line.search){
    lam1     = lambda$lambda1
    lam2     = lambda$lambda2
    lam3     = lambda$lambda3
    L1       = length(lam1)
    L2       = length(lam2)
    L3       = length(lam3)
    L        = L1+L2+L3
    aBIC     = rep(0,L)
    K.list   = rep(0,L)
    Z.list   = list()
    B.list   = list()
    V.list   = list()
    BIC.list = list()
    lam.list = list()
    member.list = list()
    cluster.matrix.list = list()
    if(L == 3){
      l=1
      aBIC = rep(0,l)
      K.list  = rep(0,l)
      lambda1 = lam1; lambda2 = lam2; lambda3 = lam3
      lam.all = c(round(lambda1,2),round(lambda2,2),round(lambda3,2))
      res.list = fusion.network(A, Omega.A, sample.index.A, gamma.index, 
                                lambda1, lambda2, lambda3, Z.int, B.int, 
                                a=a, kappa=kappa, alpha=alpha, eps = eps, 
                                niter = niter, niter.Z=niter.Z, ad.BIC=ad.BIC,
                                update.B=update.B, local.oppro=local.oppro, 
                                merge.all=merge.all)
      Z.list[[l]]   = res.list$Z.hat
      B.list[[l]]   = res.list$B.hat
      V.list[[l]]   = res.list$V.hat
      BIC.list[[l]] = res.list$BIC.vec
      aBIC[l]       = res.list$BIC.vec$BIC
      K.list[l]     = res.list$K.hat
      lam.list[[l]] = lam.all
      member.list[[l]] = res.list$cluster.results$member
      cluster.matrix.list[[l]] = res.list$cluster.results$cluster.matrix
      if(trace){
        print(c(cat(paste(l,"lam1 lam2 lam3 ="), lam.all,":"),
                paste("K =",as.numeric(res.list$K.hat))))
      }
    } else {
      # search lambda2
      lambda1 = 1;lambda3 = median(lam3)
      for (l in 1:L2) {
        lambda2 = lam2[l]
        
        lam.all = c(round(lambda1,2),round(lambda2,2),round(lambda3,2))
        res.list = fusion.network(A, Omega.A, sample.index.A, gamma.index, 
                                  lambda1, lambda2, lambda3, Z.int, B.int, 
                                  a=a, kappa=kappa, alpha=alpha, eps = eps, 
                                  niter = niter, niter.Z=niter.Z, ad.BIC=ad.BIC,
                                  update.B=update.B, local.oppro=local.oppro, 
                                  merge.all=merge.all)
        Z.list[[l]]   = res.list$Z.hat
        B.list[[l]]   = res.list$B.hat
        V.list[[l]]   = res.list$V.hat
        BIC.list[[l]] = res.list$BIC.vec
        aBIC[l]       = res.list$BIC.vec$BIC
        K.list[l]     = res.list$K.hat
        lam.list[[l]] = lam.all
        member.list[[l]] = res.list$cluster.results$member
        cluster.matrix.list[[l]] = res.list$cluster.results$cluster.matrix
        if(trace){
          print(c(cat(paste(l,"lam1 lam2 lam3 ="), lam.all,":"),
                  paste("K =",as.numeric(res.list$K.hat))))
        }
      }
      aBIC[is.na(aBIC)] = 10^10
      n_lam2 = which(aBIC[1:L2] == min(aBIC[1:L2]))[1];lambda2 = lam2[n_lam2]
      
      # search lam1
      for (l1 in 1:L1) {
        lambda1 = lam1[l1];l = L2+l1
        
        lam.all = c(round(lambda1,2),round(lambda2,2),round(lambda3,2))
        res.list = fusion.network(A, Omega.A, sample.index.A, gamma.index, 
                                  lambda1, lambda2, lambda3, Z.int, B.int, 
                                  a=a, kappa=kappa, alpha=alpha, eps = eps, 
                                  niter = niter, niter.Z=niter.Z, ad.BIC=ad.BIC,
                                  update.B=update.B, local.oppro=local.oppro, 
                                  merge.all=merge.all)
        Z.list[[l]]   = res.list$Z.hat
        B.list[[l]]   = res.list$B.hat
        V.list[[l]]   = res.list$V.hat
        BIC.list[[l]] = res.list$BIC.vec
        aBIC[l]       = res.list$BIC.vec$BIC
        K.list[l]     = res.list$K.hat
        lam.list[[l]] = lam.all
        member.list[[l]] = res.list$cluster.results$member
        cluster.matrix.list[[l]] = res.list$cluster.results$cluster.matrix
        if(trace){
          print(c(cat(paste(l,"lam1 lam2 lam3 ="), lam.all,":"),
                  paste("K =",as.numeric(res.list$K.hat))))
        }
      }
      aBIC[is.na(aBIC)] = 10^10
      n_lam1 = which(aBIC[(L2+1):(L2+L1)] == min(aBIC[(L2+1):(L2+L1)]))[1];lambda1 = lam1[n_lam1]
      
      # search lam3
      for (l3 in 1:L3) {
        lambda3 = lam3[l3];l = L2+L1+l3
        
        lam.all = c(round(lambda1,2),round(lambda2,2),round(lambda3,2))
        res.list = fusion.network(A, Omega.A, sample.index.A, gamma.index, 
                                  lambda1, lambda2, lambda3, Z.int, B.int, 
                                  a=a, kappa=kappa, alpha=alpha, eps = eps, 
                                  niter = niter, niter.Z=niter.Z, ad.BIC=ad.BIC,
                                  update.B=update.B, local.oppro=local.oppro, 
                                  merge.all=merge.all)
        Z.list[[l]]   = res.list$Z.hat
        B.list[[l]]   = res.list$B.hat
        V.list[[l]]   = res.list$V.hat
        BIC.list[[l]] = res.list$BIC.vec
        aBIC[l]       = res.list$BIC.vec$BIC
        K.list[l]     = res.list$K.hat
        lam.list[[l]] = lam.all
        member.list[[l]] = res.list$cluster.results$member
        cluster.matrix.list[[l]] = res.list$cluster.results$cluster.matrix
        if(trace){
          print(c(cat(paste(l,"lam1 lam2 lam3 ="), lam.all,":"),
                  paste("K =",as.numeric(res.list$K.hat))))
        }
      }
      aBIC[is.na(aBIC)] = 10^10
      n_lam3 = which(aBIC[(L2+L1+1):(L2+L1+L3)] == min(aBIC[(L2+L1+1):(L2+L1+L3)]))[1];lambda3 = lam3[n_lam3]
      
    }
  } else {
    lam1     = lambda$lambda1
    lam2     = lambda$lambda2
    lam3     = 0
    L1       = length(lam1)
    L2       = length(lam2)
    L3       = length(lam3)
    L        = L1*L2
    aBIC     = rep(0,L)
    K.list   = rep(0,L)
    Z.list   = list()
    B.list   = list()
    V.list   = list()
    BIC.list = list()
    lam.list = list()
    member.list = list()
    cluster.matrix.list = list()
    
    l = 1
    for (l2 in 1:L2) {
      for (l1 in 1:L1) {
        lambda1 = lam1[l1]
        lambda2 = lam2[l2]
        lambda3 = 0
        
        lam.all = c(round(lambda1,2),round(lambda2,2),round(lambda3,2))
        res.list = fusion.network(A, Omega.A, sample.index.A, gamma.index, 
                                  lambda1, lambda2, lambda3, Z.int, B.int, 
                                  a=a, kappa=kappa, alpha=alpha, eps = eps, 
                                  niter = niter, niter.Z=niter.Z, ad.BIC=ad.BIC,
                                  update.B=update.B, local.oppro=local.oppro, 
                                  merge.all=merge.all)
        Z.list[[l]]   = res.list$Z.hat
        B.list[[l]]   = res.list$B.hat
        V.list[[l]]   = res.list$V.hat
        BIC.list[[l]] = res.list$BIC.vec
        aBIC[l]       = res.list$BIC.vec$BIC
        K.list[l]     = res.list$K.hat
        lam.list[[l]] = lam.all
        member.list[[l]] = res.list$cluster.results$member
        cluster.matrix.list[[l]] = res.list$cluster.results$cluster.matrix
        if(trace){
          print(c(cat(paste(l,"lam1 lam2 lam3 ="), lam.all,":"),
                  paste("K =",as.numeric(res.list$K.hat))))
        }
        l = l+1
      }
    }
    aBIC[is.na(aBIC)] = 10^10
  }
  
  
  n_lam = which(aBIC == min(aBIC))[1]
  Opt_aBIC = min(aBIC)
  Opt_lambda = lam.list[[n_lam]]
  Opt_K = K.list[n_lam]
  Opt_Z = Z.list[[n_lam]]
  Opt_B = B.list[[n_lam]]
  Opt_member = member.list[[n_lam]]
  Opt_cluster.matrix = cluster.matrix.list[[n_lam]]
  result = list(BIC=aBIC, K.vec=K.list, Z.list=Z.list, B.list=B.list, 
                BIC.list=BIC.list, lam.list=lam.list,
                member.list=member.list, cluster.matrix.list=cluster.matrix.list,
                Opt_num=n_lam, Opt_Z=Opt_Z, Opt_B=Opt_B, 
                Opt_K=Opt_K, Opt_member=Opt_member,
                Opt_aBIC=Opt_aBIC, Opt_lambda=Opt_lambda, 
                Opt_cluster.matrix=Opt_cluster.matrix)
  return(result)
  
}

BIC.value = function(A, Z.hat, K.hat, adjust=F, B.hat){
  n = dim(Z.hat)[1]
  p = dim(Z.hat)[2]
  r.hat = 1
  
  Theta.hat = Z.hat %*% t(Z.hat)
  likelihood.mat =  log(1 + exp(Theta.hat)) - Theta.hat * A
  fitness = 2*sum(likelihood.mat[upper.tri(likelihood.mat)])/n/(n-1)
  
  if(adjust){Cn = log(n*p)} else {Cn = 1}
  degree = Cn*r.hat*K.hat*log(n)  / n
  BICvalue = fitness + degree
  return(list(BIC=BICvalue, fitness=fitness, degree=degree))
}

genelambda.obo = function(nlambda1=10,lambda1_max=1,lambda1_min=0.05,
                           nlambda2=10,lambda2_max=1,lambda2_min=0.01,
                           nlambda3=10,lambda3_max=5,lambda3_min=0.5){
  lambda1 = exp(seq(log(lambda1_max),log(lambda1_min),len= nlambda1))
  lambda2 = exp(seq(log(lambda2_max),log(lambda2_min),len= nlambda2))
  lambda3 = exp(seq(log(lambda3_max),log(lambda3_min),len= nlambda3))
  lambda = list(lambda1=lambda1,lambda2=lambda2,lambda3=lambda3)
  return(lambda)
}

genelambda.obo0 = function(nlambda1=10,lambda1_max=1,lambda1_min=0.05,
                           nlambda2=10,lambda2_max=1,lambda2_min=0.01,
                           nlambda3=10,lambda3_max=5,lambda3_min=0.5){
  lambda1 = seq(lambda1_max,lambda1_min,len= nlambda1)
  lambda2 = seq(lambda2_max,lambda2_min,len= nlambda2)
  lambda3 = seq(lambda3_max,lambda3_min,len= nlambda3)
  lambda = list(lambda1=lambda1,lambda2=lambda2,lambda3=lambda3)
  return(lambda)
}

############################# Some fundamental supporting functions ############################
gen.weight = function(A, sample.index.n, Fully.Connected=F){
  n = dim(A)[1]
  A.con = A
  if(Fully.Connected){A.con[A==0] = 1}
  omega.n = NULL
  for (i in 1:(n-1)) {
    index.i = as.matrix(sample.index.n[,sample.index.n[1,] == i])
    omega.n = c(omega.n,index.i[3,match(c(i+which(A.con[i,(i+1):n] == 1)),index.i[2,])])
  }
  Omega = Matrix(apply(sample.index.n[1:2,],2,dMatrixFun),sparse = T)
  Omega.A = Omega[,omega.n]
  sample.index.A = rbind(sample.index.n[,omega.n],1:length(omega.n))
  
  gamma.index = Matrix(0,length(omega.n),n, sparse = T)
  for (i in 1:n) {
    nn.ga.after = which(sample.index.A[1,] == i)
    nn.ga.before = which(sample.index.A[2,] == i)
    if(length(nn.ga.before)>0 & length(nn.ga.after)>0){
      gamma.index[nn.ga.after,i] = 1
      gamma.index[nn.ga.before,i] = -1
    }
    if(i==1){gamma.index[nn.ga.after,i] = 1}
    if(i==n){gamma.index[nn.ga.before,i] = -1}
  }
  
  return(list(Omega.A=Omega.A, sample.index.A=sample.index.A, gamma.index=gamma.index))
}

mcp_d = function(x,lambda,a=3){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: mcp_d
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Calculating the derivative of the MCP
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ x: a float value or a vector, the independent variable in the MCP.
  ## @ lambda: a float value, the tuning parameter in the MCP.
  ## @ a: a float value, regularization parameter in the MCP, the default setting is 3.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ rho: the derivative of the MCP.   
  ## -----------------------------------------------------------------------------------------------------------------
  
  if(lambda!=0){
    rho = lambda*( 1 > abs(x)/( lambda*a ) )*( 1 - abs(x)/( lambda*a ))
  } else{
    rho=0
  }
  return(rho)
}

S_soft = function(z,lambda){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: S_soft
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Define the soft threshold operator.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ z: a float value or a vector, the independent variable.
  ## @ lambda: a float value, the tuning parameter.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ The result of the soft threshold operator.   
  ## -----------------------------------------------------------------------------------------------------------------
  norm.z = sqrt(sum(z^2))
  if(norm.z!=0){
    n.x = 1 - lambda/norm.z
    rho = n.x*(n.x > 0)*z
  } else{
    rho = z
  }
  return(rho)
}

MCP_soft = function(z,lambda,a=3,kappa=1){
  
  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: MCP_soft
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Define the MCP threshold operator.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages: No
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ z: a float value or a vector, the independent variable in the MCP.
  ## @ lambda: a float value, the tuning parameter in the MCP.
  ## @ a: a float value, regularization parameter in the MCP, the default setting is 3.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ The result of the MCP threshold operator.   
  ## -----------------------------------------------------------------------------------------------------------------
  norm.z = sqrt(sum(z^2))
  return( S_soft(z,lambda/kappa)/(1-1/a/kappa) * (norm.z - a*lambda <= 0) + z * (norm.z - a*lambda > 0) )
}

MCP_func = function(z,lambda,a=3){
  MCP.value = (lambda*abs(z) - z^2/2/a) * (abs(z) <= lambda*a) +
    (a*lambda^2/2) * (abs(z) > lambda*a)
  return(MCP.value)
}

dMatrixFun = function(indx){
  e.vec=matrix(0,n,1)
  e.vec[indx[1],]=1
  e.vec[indx[2],]=(-1)
  return(e.vec)
}

get.cluster = function(diff.Bij, n, sample.index.n, merge.all=T){
  size = n
  diff.gfl=diff.Bij
  capgfl.matrix=Matrix(0,nrow=size,ncol=size,sparse = T)
  
  if(length(which(diff.gfl==0))==0){
    capgfl.matrix2 = capgfl.matrix+t(capgfl.matrix); diag(capgfl.matrix2)=1
    group.num.gf = size
    memb = 1:n
  }else if(length(which(diff.gfl==0))==size*(size-1)/2){
    capgfl.matrix[upper.tri(capgfl.matrix)]=1
    capgfl.matrix2 = capgfl.matrix+t(capgfl.matrix); diag(capgfl.matrix2)=1
    group.num.gf = 1
    memb = rep(1,n)
  }else{
    sample.index.gfl=sample.index.n[1:2,which(diff.gfl==0)]
    if(length(which(diff.gfl==0))==1){
      capgfl.matrix[sample.index.gfl[1],sample.index.gfl[2]]=1
    }else{
      for(i in 1:length(which(diff.gfl==0))){
        capgfl.matrix[sample.index.gfl[1,i],sample.index.gfl[2,i]]=1
      }
    }
    capgfl.matrix2 = capgfl.matrix+t(capgfl.matrix); diag(capgfl.matrix2)=1
    group.num.gf = nrow(unique(as.matrix(capgfl.matrix2)))
    
    
    cap = capgfl.matrix2
    cap.cap = apply(cap, 1, function(a){which(a == 1)})
    if(is.list(cap.cap)){num_subgroup = unique(cap.cap)}
    if(is.matrix(cap.cap)){
      uni.cap = unique(t(cap.cap))
      num_subgroup = list()
      for (i in 1:dim(uni.cap)[1]) {
        num_subgroup[[i]] = as.numeric(uni.cap[i,])
      }
      }
    non_inter_list = list()
    vv = 1
    non_inter = c(1:length(num_subgroup))
    repeat{
        a = num_subgroup[[non_inter[1]]]
        KK_k = setdiff(non_inter,non_inter[1])
        non_inter = c()
        i=1
        for (k2 in KK_k) {
          if(length(intersect(a,num_subgroup[[k2]])) > 0){
            a = union(a,num_subgroup[[k2]])
          } else {
            non_inter[i] = k2
            i=i+1
          }
        }
        non_inter_list[[vv]] = a
        vv = vv+1
        if(length(non_inter) == 0){break}
      }
    
    if(merge.all){
      for (i in 1:dim(cap)[1]) {
        for (k in 1:length(non_inter_list)) {
          if(length(intersect(which(cap[i,]==1),non_inter_list[[k]])) > 0){
            cap[i,non_inter_list[[k]]] = 1
          }
        }
      }
      capgfl.matrix2 = cap
      group.num.gf = nrow(unique(as.matrix(capgfl.matrix2)))
    }
    
    memb = rep(0,n)
    for (k in 1:length(non_inter_list)) {
      memb[non_inter_list[[k]]] = k
    }
    
  }
  
  return(list(K.hat=group.num.gf, cluster.matrix=capgfl.matrix2, member=memb))
}

memb.cluster = function(memb){
  n = length(memb)
  capgfl.matrix = matrix(0,n,n)
  for(i in 1:(n-1)){
    for (j in (i+1):(n)) {
      capgfl.matrix[i,j] <- as.numeric(memb[i] == memb[j])
    }
  }
  capgfl.matrix2 = capgfl.matrix+t(capgfl.matrix); diag(capgfl.matrix2)=1
  return(capgfl.matrix2)
}

cluster.memb = function(cluster.matrix){
  n = dim(cluster.matrix)[1]
  memb = rep(0,n)
  num_subgroup = unique(apply(cluster.matrix, 1, function(a){which(a == 1)}))
  if(is.list(num_subgroup)){
    K = length(num_subgroup)
    for (l in 1:K) {
      memb[num_subgroup[[l]]] = l
    }
  } else {
    num_subgroup = unique(t(apply(cluster.matrix, 1, function(a){which(a == 1)})))
    K = dim(num_subgroup)[1]
    for (l in 1:K) {
      memb[num_subgroup[l,]] = l
    }
  }
  return(memb)
}

############################# Functions for generating the initial values ############################
network.Z = function(A, Z.int, lambda3=0, a=3, kappa=1, alpha=1,
                     eps = 1e-2, niter = 20, niter.Z=5){
  n = as.integer(dim(A)[1])
  p = dim(Z.int)[2]
  nc = n*(n-1)
  
  Z.new = Z.int
  
  iter = 0
  eps.diff.Z = 10
  while( eps.diff.Z >= eps  && iter<niter )
  { 
    iter = iter+1
    Z.old = Z.new
    Theta.old = Z.old %*% t(Z.old)
    P.old = exp(Theta.old)/(1+exp(Theta.old))
    P.old[Theta.old > 700] = 1
    diag(P.old) = 0
    
    # update Z
    alphaZ = 2*alpha
    obj.value.old = 1
    obj.value.new = 10
    iter.Z=0
    while( obj.value.new > obj.value.old  && iter.Z<niter.Z ){
      iter.Z = iter.Z+1
      alphaZ = alphaZ/2
      g.L = - 2/n * t(A - P.old) %*% Z.old
      Z.m = sqrt(apply(Z.old*Z.old, 2, sum))
      g.J3 = t( t(Z.old) * c(mcp_d(Z.m,lambda3) / (Z.m+10^(-10))) )
      Z.new = Z.old - alphaZ*(g.L+g.J3)
      obj.value.old = obj.value.L(A, Z.old)
      obj.value.new = obj.value.L(A, Z.new)
    }
    
    eps.diff.Z = sqrt(sum((Z.new - Z.old)^2)) / max(sqrt(sum((Z.new)^2)), sqrt(sum((Z.old)^2)) )
  }

  Z.hat = Z.new
  return(list(Z.hat=Z.hat))
}

obj.value.L = function(A, ZZ){
  n = dim(ZZ)[1]
  p = dim(ZZ)[2]
  
  Theta.hat = ZZ %*% t(ZZ)
  likelihood.mat =  log(1 + exp(Theta.hat)) - Theta.hat * A
  likelihood.mat[Theta.hat > 700] = Theta.hat[Theta.hat > 700] - Theta.hat[Theta.hat > 700] * as.matrix(A)[as.matrix(Theta.hat) > 700]
  L = 2*sum(likelihood.mat[upper.tri(likelihood.mat)])/n/(n-1)
  
  obj.val = L
  return(obj.val)
}

gen.int = function(A, R=8, K.max0=8,rand.seed=123,
                   lambda3=0, a=3, kappa=1, alpha=1,
                   eps = 1e-2, niter = 20, niter.Z=5){
  n = as.integer(dim(A)[1])
  set.seed(rand.seed)
  Z.int = matrix(rnorm(n*R,0,1),n,R)
  Z.int.list = network.Z(A, Z.int, lambda3=lambda3, a=a, kappa=kappa, 
                         alpha=alpha, eps = eps, niter = niter, niter.Z=niter.Z)
  Z.int = as.matrix(Z.int.list$Z.hat)
  kmeans.clust = kmeans(Z.int, K.max0)
  memb = kmeans.clust$cluster
  B.int = Z.int
  K.memb = as.data.frame(table(memb))
  for (i in 1:dim(K.memb)[1]) {
    nemb.k = which(memb == as.numeric(K.memb[i,1]))
    if(length(nemb.k) == 1){
      B.int[nemb.k,] = t(matrix(rep(apply(t(Z.int[nemb.k,]),2,mean),length(nemb.k)),ncol=length(nemb.k)))
    } else {
      B.int[nemb.k,] = t(matrix(rep(apply(Z.int[nemb.k,],2,mean),length(nemb.k)),ncol=length(nemb.k)))
    }
    
  }
  return(list(Z.int=Z.int,B.int=B.int,memb=memb))
}

  