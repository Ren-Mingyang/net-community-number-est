################################################################################################
#Codes for conducting real data analysis
################################################################################################
rm(list = ls(all = TRUE))
ls()
setwd("./")
source("../main_functions/function.R")
source("../main_functions/alternatives.R")

###############################
library(MASS)
library(Matrix)
library(randnet)
library(pracma)
library(RSpectra)
library(pheatmap)
load("./correlation-AAL-NYU.RData")


#############################################################################
## Apply the proposed method for community detection, and Output the 
## figure of adjacency matrix of ADHD brain functional connectivity networks
## Running time: about 5 minutes
source("community_detection.R")
#############################################################################


#######################################
## some numerical features of networks 
#######################################
###### ADHD.combined
edge.com = ADHD.analysis$precision.ADHD$precision.ADHD.com
memb.com = ADHD.analysis$res.ADHD$res.ADHD.com$Opt_member
stat.num = stat.memb(memb.com, edge.com)
stat.num$K
stat.num$node.num.range # Minimum and maximum number of nodes in different communities
stat.num$edge.num.range # Minimum and maximum node degrees in different communities
table(memb.com)

###### ADHD-Inattentive
edge.ina = ADHD.analysis$precision.ADHD$precision.ADHD.ina
memb.ina0 = ADHD.analysis$res.ADHD$res.ADHD.ina$Opt_member
memb.ina = memb.ina0
memb.ina[memb.ina0==3] = 4
memb.ina[memb.ina0==4] = 3
stat.num = stat.memb(memb.ina, edge.ina)
stat.num$K
stat.num$node.num.range
stat.num$edge.num.range
table(memb.ina)

####### the clustering consistency measure between community structures of the two subtypes
cap.mat.com = memb.cluster(memb.com)
cap.mat.ina = memb.cluster(memb.ina)
n = length(memb.com)
diff.clu.mat = cap.mat.com - cap.mat.ina
1 - sum(abs(diff.clu.mat[upper.tri(diff.clu.mat)]))/(n*(n-1)/2)



###### Comparison of community structures of two subtypes and node information 
node.info = read.csv("./data/Node_AAL116_info.csv")
if(0 == sum(abs(match(ADHD.analysis$precision.ADHD$AAL.nodes, as.character(node.info$Tzourio_Mazoyer_code)) - c(1:n)))){
  part.name = unlist(lapply(strsplit(node.info$Tzourio_Mazoyer_name, "_"), function(x) x[1]))
  member = cbind(memb.com, memb.ina, part.name) ## community structures of two subtypes
  node.summary = cbind(node.info[,1:3], member, node.info[,c(7,10,6,8,9)]) ## node information 
  
  #### Output the .node and .edge files needed to draw the brain network figures
  node.cor = node.info[,c(1:5,9)]
  ###### ADHD.combined
  node.cor[,4] = memb.com
  node.cor[,6] = as.numeric(node.cor[,6])
  write.table(node.cor, "./results/com.Node_AAL.node", row.names = F, col.names = F)
  write.table(edge.com, "./results/com.Edge_AAL.edge", row.names = F, col.names = F)
  
  ###### ADHD-Inattentive
  node.cor[,4] = memb.ina
  node.cor[,6] = as.numeric(node.cor[,6])
  write.table(node.cor, "./results/ina.Node_AAL.node", row.names = F, col.names = F)
  write.table(edge.ina, "./results/ina.Edge_AAL.edge", row.names = F, col.names = F)
  
  edge.diff.com = edge.com - edge.ina
  edge.diff.ina = - edge.com + edge.ina
  edge.diff.com[edge.diff.com < 0] = 0
  edge.diff.ina[edge.diff.ina < 0] = 0
  write.table(edge.diff.com, "./results/edge.diff.com_AAL.edge", row.names = F, col.names = F)
  write.table(edge.diff.ina, "./results/edge.diff.ina_AAL.edge", row.names = F, col.names = F)
  
}






