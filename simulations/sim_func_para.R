lam.max=3
lam.min=0.5
lam1.s=2/log(n)
lam2.s=sqrt(8*log(n)/n)
lam3.s=1/8/log(n)/sqrt(n)
lambda = genelambda.obo(nlambda1=3,lambda1_max=lam.max*lam1.s,lambda1_min=lam.min*lam1.s,
                        nlambda2=10,lambda2_max=lam.max*lam2.s,lambda2_min=lam.min*lam2.s,
                        nlambda3=1,lambda3_max=lam.max*lam3.s,lambda3_min=lam.min*lam3.s)

set.seed(2)
Mu.true = matrix(rnorm(K.true*r,0,1),K.true,r)

if(sim.type=="balance"){
  unbalance = F; ka=1
  filename = paste0(c(sim.type,paste0("K",K.true),paste0("sigma0",10*rsn)), collapse = "-")
}
if(sim.type=="unbalance"){
  unbalance = T
  filename = paste0(c(sim.type,paste0("K",K.true),paste0("ka",ka),paste0("sigma0",10*rsn)), collapse = "-")
}




