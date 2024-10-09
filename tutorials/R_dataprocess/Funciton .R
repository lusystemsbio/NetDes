##  import Rcistarget funciton
library(RcisTarget)
library(readr)
library(Seurat)
library(ggplot2)
library(sctransform)
library(mgcv)
library(slingshot)
library(scales)
library(tibble)
library(stringr)
library(BiocParallel)
library(parallel)
library(PseudotimeDE)
library(multtest)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(visNetwork)
suppressWarnings(suppressPackageStartupMessages(library(ggthemes)))
suppressWarnings(suppressPackageStartupMessages(library(scran)))


BiocParallel::register(BiocParallel::SerialParam())


cluster<-function(x){
  as.matrix(table(sign(x))[1])
} ## test


genenames<-function(listgene){
  genes=list(length=length(listgene))
  for (i in 1:length(listgene)){
    gene <- gsub(" \\(.*\\). ", "; ", listgene[i], fixed=FALSE)
    k<- unique(unlist(strsplit(gene, "; ")))
    if(length(k)==0){
      genes[[i]]=c('NULL')}else{
        genes[[i]]=k
      }
    
  }
  genes}


targetnames<-function(listgene2){
  genes2=list(length=length(listgene2))
  for (i in 1:length(listgene2)){
    gene <- listgene2[i]
    k<- unique(unlist(strsplit(gene, ";")))
    if(length(k)==0){
      genes2[[i]]=c('NULL')}else{
        genes2[[i]]=k
      }
  } 
  genes2
}

updatalist<-function(gene,number,tfs){
  genef=gene
  for (i in number){
    genef[[i]]=intersect(genef[[i]],tfs)
  }
  genef
}




getnrow<-function(genes,tfs){
  
  number=c(NULL)     
  for (c in (1:length(genes))){
    if (!(length(intersect(genes[[c]],tfs))==0)){
      number=c(number,c)
    }
  }  
  number
}


tflinks<-function(number,genes,genes2){
  dtotal=matrix((NA),ncol=3)
  for(i in number){
    if(length(genes[[i]])==1){
      d=matrix(NA,ncol = 3,nrow = length(genes2[[i]]))
      d[,1]=genes[[i]]
      d[,2]=genes2[[i]]
      dtotal=rbind(dtotal,d)
    }else if (length(genes[[i]])>1){ for (t in 1:length(genes[[i]])){
      d=matrix(NA,ncol = 3,nrow = length(genes2[[i]]))
      d[,1]=genes[[i]][t]
      d[,2]=genes2[[i]]
      dtotal=rbind(dtotal,d)
    }
    } 
  }
  dtotal=dtotal[-1,]
  colnames(dtotal)=c('Source','Target','Interaction')
  dtotal
}



cluster_genes3 <- function(scdata,mscut,pseudo) {
  slope= matrix(NA,ncol=500,nrow=nrow(scdata))
  rownames(slope)=rownames(scdata)
  i = as.integer(ncol(scdata)/500)  
  for(k in 1:500){
    slope[,k]=rowMeans(scdata[,((k-1)*i+1):(k*i+1)])
  }
  slope2= matrix(NA,ncol=499,nrow=nrow(scdata))
  for(k in 1:499){
    slope2[,k]=(slope[,k+1]-slope[,k])/0.002
  }
  cluster_results <- numeric(nrow(scdata)) # Initialize the results vector
  original_change<- numeric(nrow(scdata)) 
  num_changes_all<- numeric(nrow(scdata))# Initialize the results vector
  for (r_n in 1:nrow(scdata)) {
    gene_score <- numeric(498) # Initialize the gene score vector for each gene
    for (k in 1:498) {
      gene_score[k] <- slope2[r_n, k + 1] * slope2[r_n, k]
    }
    # Identify where the slope changes sign (detecting change points)
    negative_positions <- which(gene_score < 0)
    original_change[r_n]=length(negative_positions)
    
    orichange<- length(negative_positions)
    # Reduce change points based on segment standard deviations
    num_changes <- orichange
    ms_list <- numeric(num_changes + 1)
    if (num_changes > 0) {
      split_points <- c(0, (negative_positions+1), ncol(slope2))
      cut_time_index <- c(0, (negative_positions+1)*as.integer(ncol(scdata)/500) , ncol(scdata))
      for (i in 1:(length(split_points) - 1)) {
        segment <- slope2[r_n, (split_points[i]):split_points[i + 1]]
        cut_time_number=max(pseudo[cut_time_index[i]:cut_time_index[i+1]])-min((pseudo[cut_time_index[i]:cut_time_index[i+1]]))
        
        ms_list[i] <- mean(abs(segment))*cut_time_number
        print(ms_list)
      }
      count_greater <- sum(ms_list > mscut)
      result_count <- count_greater - 1
      indices <- which(ms_list > mscut)
      if (length(indices) > 1) {
        for (j in 2:length(indices)) {
          gap_length <- indices[j] - indices[j-1] - 1  # Number of elements between the two greater elements
          if (gap_length %% 2 != 0 && gap_length > 0) {  # Check if the gap is odd and non-zero
            result_count <- result_count - 1
            break  # Adjust result once for the first odd gap found
          }
        }
      }
      num_changes=result_count
      
    }else if(num_changes==0){
      ms_list[1]=mean(abs(slope2[r_n,]))}
    for (i in 1:length(ms_list)) {
      if (ms_list[i] > mscut) {
        changeposition <- i
        break
      }
    }
    split_points2 <- c(1, (negative_positions+1), ncol(slope2))
    split_point <- split_points2[changeposition]
    split_point_next<-split_points2[changeposition+1]
    points1 <- mean(slope2[r_n,split_point:split_point_next])
    # Perform clustering based on the reduced number of changes
    if (num_changes == 0) {
      if (points1 > 0) {
        cluster_results[r_n] <- 1
      } else {
        cluster_results[r_n] <- 2
      }
    } else if (num_changes == 1) {
      if (points1>0) {
        cluster_results[r_n] <- 3}else{ cluster_results[r_n] <- 4 }
    }
    else if (num_changes ==2) {
      if (points1>0) { cluster_results[r_n] <- 5
      } else {cluster_results[r_n] <- 6}}
    else if (num_changes >2){
      cluster_results[r_n] <- 7
    }
    num_changes_all[r_n]=num_changes
  }
  return(list(cluster=cluster_results,new=num_changes_all,old=original_change))}





plot_network2 = function(tf_links = tf_links){
  # require(visNetwork)
  topology = data.frame(as.matrix(tf_links), stringsAsFactors = FALSE)
  
  node_list <- unique(c(topology[,1], topology[,2]))
  nodes <- data.frame(id = node_list, font.size = 30, value = c(rep(1, length(node_list))))  # Removed 'label' to hide node names
  
  #nodes <- data.frame(id = node_list, label = node_list, font.size =30,shape='circle',value=c(rep(1,length(node_list))))
  edge_col <- data.frame(c(1, 2), c("blue", "darkred"))
  colnames(edge_col) <- c("relation", "color")
  arrow_type <- data.frame(c(1, 2), c("arrow", "circle"))
  colnames(arrow_type) <- c("type", "color")
  
  edges <- data.frame(from = c(topology[,1]), to = c(topology[,2]),
                      arrows.to.type = arrow_type$color[c(as.numeric(topology[,3]))],
                      width = 3,
                      color = edge_col$color[c(as.numeric(topology[,3]))])
  
  visNetwork(nodes, edges, height = "1000px", width = "100%") %>%
    visEdges(arrows = "to") %>%
    visOptions(manipulation = TRUE) %>%
    visLayout(randomSeed = 123) %>%
    visPhysics(solver = "forceAtlas2Based", stabilization = FALSE)
  # file  <- paste("network_",file,".html",sep="")
  # visSave(network, file = file, selfcontained = F)
}

