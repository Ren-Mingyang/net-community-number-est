################################################################################################
#Codes for conducting simulations
################################################################################################
rm(list = ls(all = TRUE))
ls()
setwd("./")
library(MASS)
library(Matrix)
library(randnet)
library(pracma)
library(RSpectra)
source("../main_functions/function.R")
source("../main_functions/alternatives.R")
source("sim_func.R")
n = 360
sample.index.n = rbind(combn(n,2),1:(n*(n-1)/2))
r = 5



############## balance design #################
sim.type = "balance"
K.true = 4 # 2, 3, 4
rsn = 0.4  # 0.0, 0.2, 0.4
source("sim_func_para.R")


############## unbalance design #################
sim.type = "unbalance"
K.true = 2
ka = 0.1   # 0.1, 0.5, 0.9
rsn = 0.4 # 0.0, 0.4, 0.6, 0.8
source("sim_func_para.R")


# generate simulated network
sim.data.list = generate.data(n, K.true=K.true, r=r, rsn=rsn, Mu.true=Mu.true, unbalan=unbalance, ka=ka)
A = sim.data.list$A
Z.true=sim.data.list$Z.true
B.true=sim.data.list$B.true
Mu.true=sim.data.list$Mu.true
cluster.matrix.true=sim.data.list$cluster.matrix.true
Theta.true=sim.data.list$Theta.true
P.true=sim.data.list$P.true

#### proposed method ####
int.list = gen.int(A)
Z.int = int.list$Z.int
B.int = int.list$B.int
res = network.comm.num(A, sample.index.n, lambda, Z.int, B.int)
i.pro = evaluation.pro(res, Z.true, cluster.matrix.true, P.true, Theta.true, K.true)


#### alternatives ####
res.alter = alternatives(A)
index.alter = evaluation.alter(res.alter)
i.ECV.LSM = index.alter$ECV.LSM
i.NCV.LSM = index.alter$NCV.LSM
i.ECV.SSP = index.alter$ECV.SSP
i.NCV.SSP = index.alter$NCV.SSP

# True
res.LSM = LSM(A, K.hat.alter=K.true)
res.SSP = reg.SSP.SBM(A, K.hat.alter=K.true)
res.trueK = list(res.LSM=res.LSM, res.SSP=res.SSP)
i.trueK.LSM = evaluation(res.trueK$res.LSM$Z.hat, Z.true, 
                       res.trueK$res.LSM$cluster.matrix.hat, cluster.matrix.true, 
                       P.true, Theta.true, res.trueK$res.LSM$K.hat, K.true)
i.trueK.SSP =  evaluation(res.trueK$res.SSP$Z.hat, Z.true, 
                        res.trueK$res.SSP$cluster.matrix.hat, cluster.matrix.true, 
                        P.true, Theta.true, res.trueK$res.SSP$K.hat, K.true, 
                        SBM=T, res.trueK$res.SSP$P.hat.DCSBM)

i.pro           # the proposed method
i.ECV.LSM       # ECV+LSM
i.NCV.LSM       # NCV+LSM
i.ECV.SSP       # ECV+DCBM
i.NCV.SSP       # NCV+DCBM
i.trueK.LSM     # the oracle method: True K + LSM
i.trueK.SSP     # the oracle method: True K + DCBM
