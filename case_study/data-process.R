################################################################################################
#Codes for conducting real data analysis
################################################################################################
rm(list = ls(all = TRUE))
ls()
###############################
setwd("./")
parcellation_site = "AAL-NYU"
NYU_phenotypic = read.csv("./data/NYU_phenotypic.csv")
group = read.csv("./data/NYU-group.csv")
###############################
group = group[,c("Subject","Group")]
NYU_phenotypic = NYU_phenotypic[,c("ScanDir.ID","QC_Rest_1","QC_Rest_2")]




# the quality control
QC.ind = function(x){
  x = as.character(x)
  if(sum(x=="1") > 0){
    n.x = which(x==1)[1]
  } else {n.x = 0}
  return(n.x)
}
QC.mat = as.data.frame(cbind(NYU_phenotypic[,1], apply(NYU_phenotypic[,c(2,3)],1,QC.ind)))
names(QC.mat) = c("sub","QC")
for (i in 1:dim(QC.mat)[1]) {
  if(nchar(QC.mat[i,1]) == 5){
    QC.mat[i,1] = paste0("00", QC.mat[i,1])
  }
}
out.sub = "3619797"
QC.mat$QC[match(out.sub,QC.mat$sub)] = 0
sub.name.gr = group[,1]
group.na  =function(x){
  xx = x[2]
  if(nchar(xx) == 5){
    xxx = paste0("00", xx)
  } else {xxx = xx}
  return(xxx)
}
group[,1] = unlist(lapply(strsplit(sub.name.gr, "_"), group.na))

# generating correlation matrix
datapath = "./data/Time_Courses"
fileset = dir(datapath)
QC.mat = QC.mat[match(fileset, QC.mat[,1]),]
precision.list = list()
index = rep(1,length(fileset))
for (i in 1:length(fileset)) {
  filenamei = dir(paste0(datapath,"/",fileset[i]))
  QC.num = QC.mat$QC[match(fileset[i], QC.mat$sub)]
  
  if(length(filenamei) == 0){
    index[i] = 0
    print(paste0(i,"-warning"))
  } else if(QC.num > 0){
    
    if(length(filenamei) == 1){
      TC.pathfilenamei = paste0(datapath,"/",fileset[i],"/",filenamei)
    } else {
      file.list = strsplit(filenamei, "_")
      nf = length(file.list[[1]])
      sfi = unlist(file.list)
      snwmrda.i = (sub("^([[:alpha:]]*).*", "\\1", sfi[(c(1:(length(sfi)/nf))-1)*nf + 1]) == "snwmrda")
      rest.i = (sfi[(c(1:(length(sfi)/nf))-1)*nf + match("rest",file.list[[1]]) + 1] == as.character(QC.num))
      TC.filenamei = filenamei[snwmrda.i &  rest.i]
      TC.pathfilenamei = paste0(datapath,"/",fileset[i],"/",TC.filenamei)
    }
    
    TC.i = read.table(TC.pathfilenamei,head=TRUE)
    data.i = TC.i[,-c(1,2)]
    ii = match(fileset[i],group$Subject)
    
    if(dim(data.i)[1] > 30){
        precision.list[[i]] = cor(data.i)
    } 
    
  } else {
    index[i] = 0
    print(paste0(i,"-warning"))
  }
  
  
  print(i)
}

file.QC = as.data.frame(cbind(fileset,as.numeric(index)))
names(file.QC) = c("sub","QC")
group0 = group[match(file.QC$sub,group$Subject),]
group.QC = as.data.frame(cbind(group0,file.QC$QC))
names(group.QC)[3] = "QC"
match(fileset,group.QC[,1])
table(group.QC[,-1])

p = dim(precision.list[[1]])[1]
n.all = dim(group.QC)[1]
precision.array = array(0, dim = c(p, p, n.all))
for (i in 1:n.all) {
  if(length(precision.list[[i]]) > 0){
    precision.array[,,i] = precision.list[[i]]
  }
}

precision.array.ADHD.com = precision.array[,,group.QC$QC == 1 & group.QC$Group == "ADHD-Combined"]
precision.array.ADHD.hyp = precision.array[,,group.QC$QC == 1 & group.QC$Group == "ADHD-Hyperactive/Impulsive"]
precision.array.ADHD.ina = precision.array[,,group.QC$QC == 1 & group.QC$Group == "ADHD-Inattentive"]
precision.array.TDC = precision.array[,,group.QC$QC == 1 & group.QC$Group == "Typically Developing"]


precision.ADHD.com = apply(precision.array.ADHD.com, c(1,2), mean)
precision.ADHD.hyp = precision.array.ADHD.hyp
precision.ADHD.ina = apply(precision.array.ADHD.ina, c(1,2), mean)
precision.TDC = apply(precision.array.TDC, c(1,2), mean)



# binarizing by retaining connections with the average correlation matrix larger than 0.2.
precision.mat = precision.ADHD.com
(sum(abs(precision.mat) >= 0.2 & precision.mat!=0) - p)/2 / p
precision.mat[abs(precision.mat) < 0.2] = 0
precision.mat[precision.mat!=0] = 1
sort(apply(precision.mat, 1, sum), decreasing = T)
precision.ADHD.com = precision.mat
  
precision.mat = precision.ADHD.ina
(sum(abs(precision.mat) >= 0.2 & precision.mat!=0) - p)/2 / p
precision.mat[abs(precision.mat) < 0.2] = 0
precision.mat[precision.mat!=0] = 1
sort(apply(precision.mat, 1, sum), decreasing = T)
precision.ADHD.ina = precision.mat
  


AAL.nodes = unlist(lapply(strsplit(names(data.i), "_"),function(x) x[2]))
save.image(paste0("ALL.correlation-",parcellation_site,".RData"))
precision = list(precision.ADHD.com=precision.ADHD.com, precision.ADHD.hyp=precision.ADHD.hyp,
                   precision.ADHD.ina=precision.ADHD.ina, precision.TDC=precision.TDC,
                   AAL.nodes=AAL.nodes)
save(precision, file = paste0("correlation-",parcellation_site,".RData"))

