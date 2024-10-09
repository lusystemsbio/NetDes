
back_rr=readRDS("background.rds")

data(hDB)
lens = sapply(hDB, length)
ids = which(lens >= 5)
hDB = hDB[ids]
hDB=c(hDB, addhDB)

cluster.results.rr=readRDS("ref_dis_06_025.rds") # use repro_dis_06_025 and ref_dis_06_025

b2=c(NULL)
for( i in 1:length(cluster.results.rr)){
  test=cluster.results.rr[[i]]
  #test.split2<-cbind(test,1)
  #name<-str_split_fixed(test.split2[,1],"-",2)
  #test.name<-name[,2]
  b2=c(b2,test)
}




TF=c(intersect(b2,names(hDB)))
TF.hDB = hDB[TF]
final.results = c(NULL)
for (k in 1:length(TF)) {
  test.name = TF.hDB[[k]]
  listA = cluster.results.rr
  
  # Dynamically create names for listA based on its length
  names(listA) = paste0("g", 1:length(listA))
  
  nameB = test.name
  background = back_rr
  number = vector(mode = "list", length = length(listA))
  
  for (i in 1:length(listA)) {
    a = length(background[background %in% listA[[i]] & background %in% nameB])
    b = length(background[background %in% listA[[i]] & !(background %in% nameB)])
    c = length(background[background %in% nameB & !(background %in% listA[[i]])])
    d = length(background[!(background %in% c(listA[[i]], nameB))])
    
    fish.mat <- matrix(c(a, c, b, d), nrow = 2, ncol = 2)
    fish <- fisher.test(fish.mat, alternative = "greater")
    number[[i]] = fish$p.value
  }
  
  number2 = as.matrix(number)
  k1 = cbind(names(TF.hDB)[[k]], names(listA), number2)
  final.results = rbind(final.results, k1)
}

q_scale_groups=list()
for (i in 1:length(cluster.results.rr)){
  name_group = paste0("g", i)
  gn=final.results[final.results[,2]==name_group,]
  q_value=p.adjust(as.numeric(gn[,3]), method = "BH")
  q_scale_groups[[i]]=cbind(gn,q_value)
}


final.results= do.call(rbind,q_scale_groups)

tfs=final.results[final.results[,4]<=0.10,]
tfs=unique(tfs[order(unlist(tfs[,4])),])
tfs=as.character(unique(tfs[tfs[,4]<=0.10,1]))



# Define your list of gene names
tf_results <- tfs




TF.hDB2=TF.hDB[tf_results]
final.results=c(NULL)
for (k in 1:length(tf_results)){
  test.name=TF.hDB2[tf_results[k]]
  
  listA=TF.hDB2[tf_results[-k]]
  nameB=test.name[[1]]
  background=back_gens
  number= vector(mode="list",length=length(listA))
  for (i in 1:length(listA)){
    a=length(background[background%in%listA[[i]]&background%in%nameB])
    b=length(background[background%in%listA[[i]]&!(background%in%nameB)])
    c=length(background[background%in%nameB&!(background%in%listA[[i]])])
    d=length(background[!(background%in%c(listA[[i]],nameB))])
    fish.mat <- matrix(c(a,c,b,d), nrow = 2, ncol=2)
    fish <- fisher.test(fish.mat, alternative = "greater")
    number[[i]]=fish$p.value
  }
  number2=as.matrix(number)
  k1=cbind(tf_results[k],names(listA),number2)
  #k1=k1[order(as.numeric(k1[,4])),
  final.results=rbind(final.results,k1)
}


q_scale_groups=list()
for (i in 1:length(tf_results)){
  name_group=tf_results[i]
  gn=final.results[final.results[,2]==name_group,]
  q_value=p.adjust(gn[,3], method = "BH")
  q_scale_groups[[i]]=cbind(gn,q_value)
}

final.results = do.call(rbind, q_scale_groups)

final.results[final.results[,4]<=10^-20,]

head(fit.nb5000)
# Filter final.results based on the q_value condition
filtered_results <- final.results[final.results[,4] <= 10^-20, ]
gene_pairs=filtered_results[,1:2]

# Initialize a vector to store the gene with the lower sd
lower_sd_gene <- vector("character", nrow(gene_pairs))

# Loop through each pair of genes
for (i in 1:nrow(gene_pairs)) {
  gene1 <- gene_pairs[i, 1][[1]]
  gene2 <- gene_pairs[i, 2][[1]]
  
  # Calculate the standard deviation for both genes
  sd_gene1 <- sd(fit.nb5000[gene1, ])
  sd_gene2 <- sd(fit.nb5000[gene2, ])
  
  # Determine which gene has the higher and lower sd
  if (sd_gene1 > sd_gene2) {
    higher_sd_gene[i] <- gene1
    lower_sd_gene[i] <- gene2
  } else {
    higher_sd_gene[i] <- gene2
    lower_sd_gene[i] <- gene1
  }
}



filtered_results_with_sd <- cbind(filtered_results, lower_sd_gene)

genes_to_remove=unlist(unique(filtered_results_with_sd[,5]))

tf_results2 <- tf_results[!tf_results %in% genes_to_remove]


saveRDS(tf_results2,"tfs.rds")

tf_fit=fit.nb5000[tf_results2,]
write.csv(tf_fit, file = "tf_fit.csv", row.names = TRUE) #save data for fitting


