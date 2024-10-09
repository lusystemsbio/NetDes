

### load data and meta data
RAWresults<- read_csv("~/iQTL/raw_counts.csv")
meta.data<-read.delim("~/iQTL/cell_metadata_cols.tsv",stringsAsFactors= F)
names<-as.character(as.matrix(RAWresults[,1]))
RAWresults<-as.matrix(RAWresults[,2:36045])
rownames(RAWresults)=names
RAWresults=RAWresults[,rownames(meta.data)]
roundRAW=round(RAWresults)

pbmc <- CreateSeuratObject(counts = roundRAW)
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")# store mitochondrial percentage in object meta data
pbmc2 <- SCTransform(pbmc, method = "nb",ncells = 5000, variable.features.n = 3000,verbose = TRUE,return.only.var.genes = F,theta_estimation_fun="theta.mm")

rownames(roundRAW)=rownames(pbmc)
roundRAW=as.matrix(roundRAW[,rownames(meta.data)])
roundRAW.log=as.matrix(log2(roundRAW+1))
sce=SingleCellExperiment(list(counts=roundRAW,logcounts=roundRAW.log),colData=DataFrame(meta.data))
design <- model.matrix(~experiment, data = colData(sce))
# define highly variable genes
alt.fit <- modelGeneVar(sce, design = design)
top500=names(alt.fit@metadata[["var"]][order(alt.fit@listData[["bio"]],decreasing = TRUE)])[1:500]
top5000=names(alt.fit@metadata[["var"]][order(alt.fit@listData[["bio"]],decreasing = TRUE)])[1:5000]

logu=pbmc2@assays[["SCT"]]@scale.data  ### dataset
logu.500=logu[which(rownames(logu)%in%top500),]#get the top500 gene dataset
save(logu.500, file = "top500pbmc.RData")

sce<-SingleCellExperiment(list(logcounts=logu.500))
#rd <- irlba::prcomp_irlba(t(logcounts(sce)), n = 5, scale. = FALSE)$x[, 1:5]
pca_res <- prcomp(t(logcounts(sce)), center = TRUE,scale.=TRUE)
rd <- pca_res$x[, 1:5]

reducedDims(sce) <- SimpleList(PCA = rd)
colData(sce)$cl <- 1
fit_ori <- slingshot(sce, reducedDim = 'PCA', clusterLabels = "cl")
#fit_ori <- slingshot(rd, reducedDim = 'PCA')

pseudotime.u = rescale(colData(fit_ori)$slingPseudotime_1)
#ori_tbl <- tibble(cell = colnames(sce), pseudotime = rescale(colData(fit_ori)$slingPseudotime_1))
ori_tbl <- tibble(cell = colnames(sce), pseudotime = rescale(colData(fit_ori)$slingPseudotime_1))



####
#### delete the end points of pseudotime

cell.pick=ori_tbl[(ori_tbl[,2]<0.90)&(ori_tbl[,2]>0.05),1]

cell.pick=as.character(as.matrix(cell.pick))
meta.data2=meta.data[cell.pick,]
fit_ori2=colData(fit_ori)[rownames(colData(fit_ori))%in%cell.pick,]
ori_tbl2=tibble(cell = cell.pick, pseudotime = rescale(fit_ori2$slingPseudotime_1))
pseudotime.u = rescale(fit_ori2$slingPseudotime_1)


## plot pseudotime
p<-ggplot(data.frame(ori_tbl2), aes(x=pseudotime, y=0,colour= as.factor(meta.data2$day))) +
  geom_point(size = 0.1) +labs(x="psedotime",y="")+ theme(panel.background = element_blank(),
                                                          axis.text = element_blank(),
                                                          axis.ticks = element_blank(),)
pdf("testplot/slingshot_PC5_2.pdf", width = 8, height = 4)
p
dev.off()
###

roundRAW=roundRAW[,cell.pick]
roundRAW.log=as.matrix(log2(roundRAW+1))

RAW.sce=SingleCellExperiment(list(counts=roundRAW,logcounts=roundRAW.log))
system.time(res_fix.raw5000 <- PseudotimeDE::runPseudotimeDE(gene.vec = top5000,
                                                             ori.tbl = ori_tbl2,
                                                             sub.tbl = NULL, # Set as NULL to only get fix.pv
                                                             mat = RAW.sce,
                                                             model = "nb",mc.cores = 4,k=4,knots=c(0:3/3))) 


fit.nb5000=matrix(NA,ncol=ncol(RAW.sce),nrow =5000 )
colnames(fit.nb5000)=rownames(meta.data2)
rownames(fit.nb5000)=as.matrix(res_fix.raw5000[,1])
for( i in 1:5000){
  fit.nb5000[i,]=predict(res_fix.raw5000[[8]][[i]])
}
## get scale results 
save(fit.nb5000,file="fitnb30000cell2.RData")
df=fit.nb5000
for (i in 1:nrow(df)){
  df[i,]=scale(df[i,],center = T,scale = T)
}
#### get dif gene
ord.u=ori_tbl2[order(pseudotime.u),]
ord.u=ord.u[,1]
ord.u=as.character(as.matrix(ord.u))

sdvar.same=matrix(NA,ncol=2,nrow =nrow(fit.nb5000))
colnames(sdvar.same)=c("sd","var")
rownames(sdvar.same)=rownames(fit.nb5000)
for (i in 1:nrow(fit.nb5000)){
  test=fit.nb5000[i,]
  sdvar.same[i,1]=sd(test)
  sdvar.same[i,2]=var(test)
}
test=sdvar.same[order(sdvar.same[,1]),]
fitorder=rownames(fit.nb5000[order(sdvar.same[,1]),])
gene.dif=rownames(sdvar.same[which((sdvar.same[,1]>0.25)),])#### choose a cut off by sd 
gene.same=rownames(sdvar.same[which((sdvar.same[,1]<=0.25)),])

###
df.pse=df[,ord.u]
df.pse=df.pse[gene.dif,]

psu_ord=pseudotime.u[order(pseudotime.u)]
result=cluster_genes3(df.pse,mscut=0.2,pseudo=psu_ord)

gene_cluster=cbind(rownames(df.pse),result[[1]])

cluster.results=list()
cluster.results[[1]]=gene_cluster[gene_cluster[,2]==1,1]
cluster.results[[2]]=gene_cluster[gene_cluster[,2]==2,1]
cluster.results[[3]]=gene_cluster[gene_cluster[,2]==3,1]
cluster.results[[4]]=gene_cluster[gene_cluster[,2]==4,1]
cluster.results[[5]]=gene_cluster[gene_cluster[,2]==5,1]
cluster.results[[6]]=gene_cluster[gene_cluster[,2]==6,1]
cluster.results[[7]]=gene_cluster[gene_cluster[,2]==7,1]

saveRDS(cluster.results,file="2geneclusters_dis_02_025.rds")  ## 2 _ use8 as cut , 3_use 6 as cut
saveRDS(pseudotime.u,file="pseudotime.rds")
saveRDS(ori_tbl2,file="ori_tbl.rds")
#ori_tbl2=readRDS("ori_tbl.rds")
#pseudotime.u=readRDS("pseudotime.rds")
#load("fitnb30000cell2.RData")
#cls=readRDS(file="2geneclusters_dis_12_025.rds")


nrows <- c()
for (i in 0:200) {
  closest_to_i <- which.min(abs(pseudotime.u - i*0.005))
  nrows <- c(nrows, closest_to_i)
}

pseudotime.pick=pseudotime.u[nrows]
fit.nb5000_pick=fit.nb5000[,nrows]
df=fit.nb5000_pick
for (i in 1:nrow(df)){
  df[i,]=scale(df[i,],center = T,scale = T)
}



pdf("testplot/cluster_iPSC_sd002_contour.pdf", width = 12, height = 6)
par(mfrow=c(3,3), mar=c(5,5,2,2))
for (ii in 1:7){
  need=cluster.results[[ii]]
  fit=df[need,]
  fit <- data.frame(fit)
  fit_vector=as.vector(as.matrix(fit))
  x=rep(pseudotime.pick,each=length(cluster.results[[ii]]))
  y=fit_vector
  x <- as.numeric(as.vector(x))
  y <- as.numeric(as.vector(y))
  dat <- data.frame(x = x, y = y)
  den=kde2d(dat$x, dat$y,n=100)
  avg_traj <- apply(fit, 2, mean)
  plot(NA, xlim=c(0,1), ylim=c(-3,3), xlab="Pseudotime", ylab="log2u", cex.lab=1.5, cex.axis=1.5)
  # Check the class of 'den'
  class(den)
  # If 'den' is a list of vectors, extract the first vector
  myvector <- den[[1]]
  # Find the maximum value of the vector
  my_max <- max(myvector)
  # Modify the contour function to use the correct levels
  contour(den, xlab="Pseudotime", ylab="Y", drawlabels=FALSE, add=TRUE, lwd=0.5, levels=seq(0.1, my_max, length.out=10))
  lines(pseudotime.pick, avg_traj, col="red", lwd=2)
  title(paste("cluster", ii), cex.main=2)
  legend("topright", legend=c("Density plot", "Average trajectory"), lty=c(1, 1), col=c(1, "red"), lwd=c(1, 2), cex=0.5)                                                                                                                  
}
dev.off()




pdf("testplot/cluster_iPSC_sd25_dis02.pdf", width = 8, height = 9)
par(mfrow=c(3,3), mar=c(5,5,2,2))
for (ii in 1:7) {
  need = cluster.results[[ii]]
  fit = df[need, ]
  fit = data.frame(fit)
  
  x = rep(pseudotime.pick, each = length(cluster.results[[ii]]))
  y = as.vector(as.matrix(fit))
  
  x = as.numeric(as.vector(x))
  y = as.numeric(as.vector(y))
  
  plot(NA, xlim = c(0, 1), ylim = range(y), xlab = "Pseudotime", ylab = "log2u", cex.lab = 1.5, cex.axis = 1.5)
  
  # Add lines for each gene
  for (gene in 1:nrow(fit)) {
    lines(pseudotime.pick, fit[gene, ], col = "black", lwd = 0.5)
  }
  
  # Calculate and plot the average trajectory
  avg_traj = apply(fit, 2, mean)
  lines(pseudotime.pick, avg_traj, col = "red", lwd = 2)
  
  # Add the number of genes to the upper right of the
  
}
dev.off()
