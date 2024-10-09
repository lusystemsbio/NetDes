

tfs <- readRDS("~/Desktop/8.17/Singlecell6.25/7.27/fisher exact test/guohaotfs.rds")


data(motifAnnotations_hgnc)
motifRankings <- importRankings("~/Desktop/8.17/Rcistarget test/hg19-tss-centered-10kb-7species.mc9nr.feather")

motiftable <- list()
motifwgene <- list()

for(i in 1:length(cluster.results)) {
  tryCatch({
    motifEnrichmentTable_wGenes <- cisTarget(cluster.results[[i]], motifRankings,
                                             motifAnnot=motifAnnotations_hgnc)
    
    motifs_AUC <- calcAUC(cluster.results[[i]], motifRankings)
    
    # 2. Select significant motifs, add TF annotation & format as table
    motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, 
                                               motifAnnot=motifAnnotations_hgnc, nesThreshold = 1)
    
    # 3. Identify significant genes for each motif
    motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable, 
                                                       geneSets=cluster.results[[i]],
                                                       rankings=motifRankings, 
                                                       nCores=1,
                                                       method="aprox")
    
    motiftable[[i]] <- motifEnrichmentTable
    motifwgene[[i]] <- motifEnrichmentTable_wGenes
  }, error = function(e) {
    message(paste("Error in iteration", i, ":", e$message))
  })
}


motifwgenetotal=do.call(rbind,motifwgene)

motiftabletotal=do.call(rbind,motiftable)



library(dplyr)
library(NetAct)
tfs=tfs
lim=1
ae.wgene=motifwgenetotal[motifwgenetotol$NES>lim,]
ae=motiftabletotal[motiftabletotal$NES>lim,]
TFhighConf=split(ae.wgene$TF_highConf,
                 ae$geneSet)
listgene=TFhighConf$geneSet

TFhighConf2=split(ae.wgene$enrichedGenes,
                  ae$geneSet)
listgene2=TFhighConf2$geneSet

genestest=genenames(listgene = listgene)
genestest2=targetnames(listgene2 = listgene2)

numbertest=getnrow(genestest,tfs=tfs)
genestest=updatalist(gene = genestest,number = numbertest,tfs = tfs)
genestest2=updatalist(gene = genestest2,number = numbertest,tfs = tfs)

tflinks.results=tflinks(numbertest,genes = genestest,genes2=genestest2)

d[,1]=genes[[i]]

tflinks<-function(number,genes,genes2){
  dtotal=matrix((NA),ncol=3)
  for(i in number){
    if(length(genes[[i]])==1 && !is.null(genes2[[i]])){
      d=matrix(NA,ncol = 3,nrow = length(genes2[[i]]))
      d[,1]=genes[[i]]
      d[,2]=genes2[[i]]
      dtotal=rbind(dtotal,d)
    }else if (length(genes[[i]])>1 && !is.null(genes2[[i]])){ for (t in 1:length(genes[[i]])){
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

dtotal


tflinks.results=as.data.frame(tflinks.results)
dtotal.unique=distinct(tflinks.results, .keep_all = TRUE)
dtotal.unique[,3]=1
plot_network(dtotal.unique)
####


dtotal.unique


write.csv(dtotal.unique, "network1.8.csv")

## add trrust database


trrust_ <- read.table('trrust_rawdata.human.tsv', sep='\t', header=F)
colnames(trrust_)=c('Source','Target','Interaction','n')

trrust_tf=trrust_[trrust_[,1]%in%tfs,]
trrust_tf_2=trrust_tf[trrust_tf[,2]%in%tfs,]  
trrust_tf_3=trrust_tf_2[!(trrust_tf_2[,1]==trrust_tf_2[,2]),]
interactions=unique(trrust_tf_3[,1:2])
inter=cbind(interactions,1)
colnames(inter)=colnames(dtotal.unique)


### NetAct
load("~/Desktop/8.17/Singlecell6.25/7.27/fisher exact test/hDB2.RData")
TF_hDB2=hDB[tfs]
TF_hDB=hDB[tfs]

for (i in 1:length(TF_hDB)){
  TF_hDB2[[i]]=TF_hDB[[i]][TF_hDB[[i]]%in%tfs]
}

inter_NATACT=matrix(NA,nrow=1,ncol=3)
inter_NATACT[1,]=c('Source','Target','Interaction')
for (i in 1:length(TF_hDB2)){
  if (length(TF_hDB2[[i]])>0){
    interaction=matrix(NA,nrow=length(TF_hDB2[[i]]),ncol=3)
    interaction[,1]=names(TF_hDB2[i])
    interaction[,2]=TF_hDB2[[i]]
    interaction[,3]=1
    inter_NATACT=rbind(inter_NATACT,interaction)
  }
}


colnames(inter_NATACT)=inter_NATACT[1,]
inter_NATACT=inter_NATACT[-1,]
inter_NATACT=as.data.frame(inter_NATACT)

dtotal.unique2=unique(rbind(dtotal.unique,inter_NATACT))

### NetAct

df=dtotal.unique2

plot_network(dtotal.unique2)

only_source <- setdiff(unique(df$Source), unique(df$Target))

# Exclude SOX17 and ETV2 from this list
only_source <- setdiff(only_source, c("SOX17", "ETV2"))

# Remove rows where the Source is in the only_source list
df_filtered <- df[!df$Source %in% only_source, ]

plot_network(df_filtered)

plot_network2(df_filtered)
