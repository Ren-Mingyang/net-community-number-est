precision.ADHD.com = precision$precision.ADHD.com
precision.ADHD.ina = precision$precision.ADHD.ina
n = dim(precision.ADHD.com)[1]
sample.index.n = rbind(combn(n,2),1:(n*(n-1)/2))

## the function saving heatmap figure
save_pheatmap = function(x, filename, width=8, height=8) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

#####################################
########## ADHD.combined ############
K.max = 8
R = 10
# ADHD.com
precision.mat = precision.ADHD.com
int.list = gen.int(precision.mat, R=R, K.max0=K.max)
Z.int = int.list$Z.int
B.int = int.list$B.int
lam1.s=2/log(n)
lam2.s=sqrt(log(n)/n)
lam3.s=1e-10
lam.max1=5
lam.min1=0.5
lam.max2=3
lam.min2=2.5
nlam1 = 3
nlam2 = 8
lambda = genelambda.obo0(nlambda1=nlam1,lambda1_max=lam.max1*lam1.s,lambda1_min=lam.min1*lam1.s,
                         nlambda2=nlam2,lambda2_max=lam.max2*lam2.s,lambda2_min=lam.min2*lam2.s,
                         nlambda3=1,lambda3_max=lam3.s,lambda3_min=lam3.s)
t1 = proc.time()
res = network.comm.num(precision.mat, sample.index.n, lambda, Z.int, B.int, line.search = F)
t2 = proc.time()
t2 - t1

n.o = res$Opt_num
res$K.vec[n.o]
Opt_member = res$member.list[[n.o]]
table(Opt_member)

precision.plot = precision.mat
memb = Opt_member
t.memb = as.data.frame(table(memb))
t.memb = t.memb[order(t.memb$Freq,decreasing = T),]
memb.vec = NULL
for (k in 1:dim(t.memb)[1]) {
  memb.k = which(memb == t.memb[k,1])
  memb.vec = c(memb.vec, memb.k)
  precision.k = precision.plot[memb.k, memb.k]
  if(dim(t.memb)[1] %% 2 == 1){
    precision.k[precision.k!=0] = (-1)^k * ceiling(k/2) + as.numeric((-1)^k > 0) +(-1)^k*0.5
    precision.plot[memb.k, memb.k] = precision.k
  } else {
    precision.k[precision.k!=0] = (-1)^k * ceiling(k/2) + 2*as.numeric((-1)^k > 0)-1 +(-1)^k*0.5
    precision.plot[memb.k, memb.k] = precision.k
  }
  
}
A.memb = precision.plot[memb.vec,memb.vec]
table(A.memb)
library(pheatmap)
map = pheatmap(A.memb, cluster_cols =F, cluster_rows =F, legend = FALSE)
save_pheatmap(map, "ADHD.adjacency.com.eps")
res.ADHD.com = res






#####################################
######### ADHD-Inattentive ##########
K.max = 8
R = 10
# ADHD.ina
precision.mat = precision.ADHD.ina
int.list = gen.int(precision.mat, R=R, K.max0=K.max)
Z.int = int.list$Z.int
B.int = int.list$B.int
lam1.s=2/log(n)
lam2.s=sqrt(log(n)/n)
lam3.s=1e-10
lam.max1=5
lam.min1=0.5
lam.max2=3.75
lam.min2=2.5
nlam1 = 3
nlam2 = 8
lambda = genelambda.obo0(nlambda1=nlam1,lambda1_max=lam.max1*lam1.s,lambda1_min=lam.min1*lam1.s,
                         nlambda2=nlam2,lambda2_max=lam.max2*lam2.s,lambda2_min=lam.min2*lam2.s,
                         nlambda3=1,lambda3_max=lam3.s,lambda3_min=lam3.s)
t1 = proc.time()
res = network.comm.num(precision.mat, sample.index.n, lambda, Z.int, B.int, line.search = F)
t2 = proc.time()
t2 - t1

n.o = res$Opt_num
res$K.vec[n.o]
Opt_member = res$member.list[[n.o]]
table(Opt_member)

precision.plot = precision.mat
memb = Opt_member
t.memb = as.data.frame(table(memb))
t.memb = t.memb[order(t.memb$Freq,decreasing = T),]
memb.vec = NULL
for (k in 1:dim(t.memb)[1]) {
  memb.k = which(memb == t.memb[k,1])
  memb.vec = c(memb.vec, memb.k)
  precision.k = precision.plot[memb.k, memb.k]
  if(dim(t.memb)[1] %% 2 == 1){
    precision.k[precision.k!=0] = (-1)^k * ceiling(k/2) + as.numeric((-1)^k > 0) +(-1)^k*0.5
    precision.plot[memb.k, memb.k] = precision.k
  } else {
    precision.k[precision.k!=0] = (-1)^k * ceiling(k/2) + 2*as.numeric((-1)^k > 0)-1 +(-1)^k*0.5
    precision.plot[memb.k, memb.k] = precision.k
  }
  
}
A.memb = precision.plot[memb.vec,memb.vec]
table(A.memb)
library(pheatmap)
map = pheatmap(A.memb, cluster_cols =F, cluster_rows =F, legend = FALSE)
save_pheatmap(map, "ADHD.adjacency.ina.eps")
res.ADHD.ina = res


res.ADHD = list(res.ADHD.com=res.ADHD.com,
                res.ADHD.ina=res.ADHD.ina)
ADHD.analysis = list(precision.ADHD=precision, res.ADHD=res.ADHD) 


save(ADHD.analysis, file = paste0("./results/correlation-ADHD.results.analysis.RData"))



stat.memb = function(memb, edge){
  t.memb = as.data.frame(table(memb))
  t.memb = t.memb[order(t.memb$Freq,decreasing = T),]
  K = dim(t.memb)[1]
  # edge.num: Degrees per node
  # edge.num.range: Minimum and maximum node degrees in different communities
  # node.num.range: Minimum and maximum number of nodes in different communities
  # within.edge.average: Average degree of nodes in the community
  # betw.edge.average: Average degree of nodes between communities
  
  edge.num = apply(edge,1,sum) - 1
  edge.num.list = list()
  edge.num.range = as.data.frame(matrix(0, K, 2))
  for (k in 1:K) {
    memb.k = which(memb == t.memb[k,1])
    edge.num.list[[k]] = edge.num[memb.k]
    edge.num.range[k,] = c(min(edge.num[memb.k]),max(edge.num[memb.k]))
  }
  row.names(edge.num.range) = as.character(t.memb[,1])
  names(edge.num.range) = c("min", "max")
  node.num.range = c(min(t.memb[,2]),max(t.memb[,2]))
  stat.num = list(K=K, node.num.range=node.num.range, edge.num.range=edge.num.range,
                  edge.num.list=edge.num.list, t.memb=t.memb)
  return(stat.num)
}


