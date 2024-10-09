library(readr)
library(SacoGraci)
Network_results<- read_csv("network_noadd_33.csv") 

Network_results

for( i in 1:nrow(Network_results)){
  if (Network_results[i,3]>1){Network_results[i,3]=1}
  else if (Network_results[i,3]<1){Network_results[i,3]=2}
  else{Network_results[i,3]=0}
}

Network_results=as.data.frame(Network_results[,1:3])
colnames(Network_results)=c("Source","Target","Type")
gene_list = list()
gene_list[[1]]=c("FOS","HES1","LEF1")
gene_list[[2]]=c("POU5F1","TWIST2","SOX2","MYCN","E2F4")
gene_list[[3]]=c("GATA3","RARA","ZEB1","KLF8")
data_exp = gen_RACIPE(Network_results, nModels=1000) ## nModels same to the modelsCGr in opt_SA
output_model_grouping = modClustHCA(logscData = data_exp,numbClust = 3)  # model clustering with HCA, we specify 3 clusters
data_reordered = output_model_grouping$dataRearr
clusterRef = output_model_grouping$clInd  # cluster indices of all models
output_gene_grouping = geneClustInd(logscData = data_reordered,numbGeneClust = 3)
output_processed = reordering(logscData = data_reordered, gene_list = gene_list) # reordered gene expression data by the clustering outcome
data_processed = output_processed$data
stat_clusters <-centMedVarCutDistPerc(data = data_processed, clusterRef = clusterRef, percThr = 0.01)
inTopsM <-gaInitial_gen(circuit_top = Network_results, gene_list = gene_list, numbNewTop = 90)

circuit4 = opt_SA(network_top = Network_results, data = data_processed, clusterRef = clusterRef, 
                  cenMedRef = stat_clusters$center, cutOffM = stat_clusters$radius, 
                  gene_list = gene_list, init_top = inTopsM[1,], 
                  output = "Results4_3", nRepeat= 1, modelsCGr = 1000, 
                  maxT=400, threshT=40, decayRate1=0.8, decayRate2=0.6, iter_per_temp=3)## modelsCGr same to the nModels in gen_RACIPE

write.csv(circuit4, "circuit42.csv", row.names = FALSE)