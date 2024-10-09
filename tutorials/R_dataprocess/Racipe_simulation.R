# Initialize an empty list to store the data
Networklist_noadd <- list()

# Define the common file path
file_path <- "~/Desktop/Singlecell6.25/Finally fitting/network/"

# Loop through the file indices and read the CSV files
for (i in 31:45) {
  file_name <- paste0(file_path, "network_noadd_", i, ".csv")
  Networklist_noadd[[46 - i]] <- read_csv(file_name)
}



NetDes_noadd_results=list()



for (i2 in 1:15){
  network=Networklist_noadd[[i2]]
  for( i in 1:nrow(network)){
    if (network[i,3]>1){network[i,3]=1}
    else if (network[i,3]<1){network[i,3]=2}
    else{network[i,3]=0}
  }
  network=as.data.frame(network[,1:3])
  colnames(network)=c("Source","Target","Type")
  plot_network(network)
  rSet <- sRACIPE::sracipeSimulate(circuit = network, numModels = 10000,
                                   plots = FALSE, integrateStepSize = 0.1, 
                                   simulationTime = 30)
  
  NetDes_noadd_results[[i2]]=rSet
}



numeric_vector <- 2^seq(from = -10, to = 5, by = 1)


dis_scale_NetDesnoadd=numeric(15)
percentage_NetDesnoadd=matrix(NA,ncol = 15,nrow=7)
percentage_NetDesnoadd_3D=matrix(NA,ncol = 15,nrow=7)
umap_vector_NetDesnoadd=matrix(NA,nrow=6,ncol=15)
umap_vector_random_NetDesnoadd=matrix(NA,nrow=6,ncol=15)
results_NetDesnoadd=matrix(NA,ncol = 15,nrow=7)

for (i in 1:15){
  rSet=NetDes_noadd_results[[i]]
  
  dat1=rSet@assays@data@listData[["deterministic"]]
  
  tf12RARA2=tf12RARA[rownames(dat1),]
  
  for (i2 in 1:nrow(dat1)){
    tf12RARA2[i2,]=scale(as.numeric(tf12RARA2[i2,]))
  }
  
  tf12RARA3 <- matrix(as.numeric(tf12RARA2),    # Convert to numeric matrix
                      ncol = ncol(tf12RARA2))
  rownames(tf12RARA3)=rownames(tf12RARA2)
  colnames(tf12RARA3)=colnames(tf12RARA2)
  
  opt_results=ob_bestscalevalue(RACIPE_data = dat1, real_data = tf12RARA3,numeric_vector=numeric_vector)
  
  for (ib in 1:nrow(dat1)){
    
    dat1[ib,]=log2(dat1[ib,])
  }
  stats_scale_2=tf12RARA3[,c(1,40,80,120,160,201)]
  dat_0_1=t(scale(t(dat1)))
  
  pca_sim=prcomp(t(dat_0_1))
  pca_percentage = pca_sim$sdev^2 / sum(pca_sim$sdev^2)
  pca_percentage = pca_percentage * 100
  
  pca_df=data.frame(PC1 = pca_sim$x[,1], PC2 = pca_sim$x[,2],PC3=pca_sim$x[,3])
  p <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
    geom_pointdensity(size = 3, shape = 17) +
    scale_color_viridis_c(guide = guide_colorbar(title = "# of neighboring cells")) + 
    xlab(paste("PC1 (", round(pca_percentage[1], 2), "%)", sep = "")) +
    ylab(paste("PC2 (", round(pca_percentage[2], 2), "%)", sep = "")) +
    theme_bw()
  
  
  p2 <- ggplot(pca_df, aes(x = PC1, y = PC3)) +
    geom_pointdensity(size = 3, shape = 17) +
    scale_color_viridis_c(guide = guide_colorbar(title = "# of neighboring cells")) + 
    xlab(paste("PC1 (", round(pca_percentage[1], 2), "%)", sep = "")) +
    ylab(paste("PC3 (", round(pca_percentage[3], 2), "%)", sep = "")) +
    theme_bw()
  
  set.seed(1)
  clust=cluster_nor(10000,stats=stats_scale_2,dat_RACIPE = dat_0_1,p=0.05)
  cutoffs=cluster_cutoffs(10000,stats=stats_scale_2,dat_RACIPE = dat_0_1,p=0.05)
  #clust_all2[[i]]=clust
  #clust=clust_all2[[i]]
  #clust=clust_all[[i]]
  percentage_NetDesnoadd[,i]=calculate_density_ratios(dat_0_1 =dat_0_1 ,k=15,clust=clust)
  results_NetDesnoadd[,i] <- c(length(which(clust == 0)) / 10000,length(which(clust == 1)) / 10000,length(which(clust == 2)) / 10000,
                                length(which(clust == 3)) / 10000,length(which(clust == 4)) / 10000,length(which(clust == 5)) / 10000,
                                length(which(clust == 6)) / 10000)
  percentage_NetDesnoadd_3D[,i]=calculate_density_ratios(dat_0_1 =t(pca_df) ,k=15,clust=clust)
  calculate_density_ratios(dat_0_1 =t(pca_df) ,k=15,clust=clust)
  colours_pick <- c(1, '#1f78b4', '#ff7f00', '#33a02c', '#6a3d9a', '#fb9a99','#a6a6a6')
  plot(pca_df[, 1], pca_df[, 2], col = colours_pick[clust + 1], xlab = "PC1", ylab = "PC2")
  
  colours_pick <- c(1, '#1f78b4', '#ff7f00', '#33a02c', '#6a3d9a', '#fb9a99', '#a6a6a6')
  
  # Create a ggplot object
  eigenvectors=pca_sim[["rotation"]]
  pca2= t(tf12RARA3) %*% eigenvectors
  PCAresults2=as.data.frame(pca2)
  
  
  p_line=p+   geom_path(data = PCAresults2, aes(x = PC1, y = PC2),col=2)
  
  p_line 
  p2_line=p2+   geom_path(data = PCAresults2, aes(x = PC1, y = PC3),col=2)
  p2_line
  plot(pca_df[, 1], pca_df[, 3], col = colours_pick[clust + 1], xlab = "PC1", ylab = "PC3")
  plot(pca_df[, 1], pca_df[, 2], col = colours_pick[clust + 1], xlab = "PC1", ylab = "PC2")
  ummap_results=non_map_scale_distance(clust=clust,stats=stats_scale_2,dat=dat_0_1,cutoffs=cutoffs)
  ummap_distance=ummap_results[[1]]
  umap_vector_NetDesnoadd[,i]=ummap_results[[2]]
  set.seed(1)
  data_random <- random_points(dat_0_1,10000)
  clust2=cluster_nor(10000,stats=stats_scale_2,dat_RACIPE = data_random,p=0.05)
  unmap_random_results=non_map_scale_distance(clust=clust2,stats=stats_scale_2,dat=data_random,cutoffs=cutoffs)
  unmap_random_dis=unmap_random_results[[1]]
  dis_scale_NetDesnoadd[i]=(ummap_distance-1)/(unmap_random_dis-1)
  umap_vector_random_NetDesnoadd[,i]=unmap_random_results[[2]]
}

entropy_network=results_NetDesnoadd[2:nrow(results_NetDesnoadd),]
NetDesnoadd_entropy=infor_matrix(entropy_network)



