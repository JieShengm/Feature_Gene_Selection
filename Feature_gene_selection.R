rm(list=ls())
library("DuoClustering2018")
library("EnsDb.Hsapiens.v86")
library("Seurat")
library("scran")
library("scuttle")
library("EDGE")
library("scry")
library("BASiCS")
library("SLICER")
library("M3Drop")
library("HIPPO")
library("calibrate")
library("ggplot2")
library("forcats")
library("mclust")


## Dataset
sce=sce_full_Zhengmix8eq(metadata = FALSE) 
counts_mat=assay(sce,"counts")
annotation=sce$phenoid

## Main function
FeatureGeneSelection=function(sce, method, seed=2021, n_feature=NULL, p_thresh=0.05){
  set.seed(seed)
  result=list()
  counts_mat=assay(sce,"counts")
  
  start_time = Sys.time()
  
  if(method%in%c("disp","mvp","vst","sctransform")){
    data <- CreateSeuratObject(counts = counts_mat)
    if(method!="sctransform"){
      data = NormalizeData(data,verbose = FALSE)
      if(is.null(n_feature)){n_feature=2000}
      data = FindVariableFeatures(data, selection.method = method, nfeatures=n_feature, verbose = FALSE)
      genelist = data@assays[["RNA"]]@var.features
    }else{
      if(is.null(n_feature)){n_feature=3000}
      data = SCTransform(data,variable.features.n = n_feature, verbose = FALSE)
      genelist = data@assays[["SCT"]]@var.features
      }
  }else if(method%in%c("scran")){
    scran_sce=sce
    scran_sce = logNormCounts(scran_sce)
    dec = modelGeneVar(scran_sce)
    top.hvgs = getTopHVGs(dec) ## keep all HVGs without filtering in this step, filter when calculating correlation.
    cor.pairs = correlatePairs(scran_sce, subset.row=top.hvgs)
    if(is.null(n_feature)){
      sig.cor = cor.pairs$FDR <= p_thresh
      genelist = unique(c(cor.pairs$gene1[sig.cor], cor.pairs$gene2[sig.cor]))  ## fdr.threshold=0.05
    }else{
      cor.genes = correlateGenes(cor.pairs)
      genelist = cor.genes$gene[1:n_feature] 
      }
  }else if(method%in%c("scVEGs")){
    if(is.null(n_feature)){
    scVEGs_mat=counts_mat
    ensembl.genes = rownames(sce)
    geneIDs = ensembldb::select(EnsDb.Hsapiens.v86, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
    sub_ind=c()
    for(i in 1:dim(geneIDs)[1]){
      sub_ind=c(sub_ind,which(rownames(scVEGs_mat)==geneIDs$GENEID[i]))
    }
    scVEGs_mat = scVEGs_mat[sub_ind,]
    rownames(scVEGs_mat) = geneIDs$SYMBOL
    pVal = 0.05
    pFlag = 1
    species = 'hsa'
    cellSum = apply(scVEGs_mat, 2, sum)
    scaleNum = mean(cellSum)/cellSum
    scaleFactor = t(kronecker( matrix(1, 1, dim(scVEGs_mat)[1]), scaleNum))
    processed_mat = scaleFactor * scVEGs_mat
    outputName="result"
    source("scVEGs.r")
    sig = scVEGs(processed_mat, pVal, pFlag, species, outputName)
    genelist = geneIDs$GENEID[geneIDs$SYMBOL%in%rownames(sig)]
    }
  }else if(method%in%c("NBDisp","NBDrop")){
    if(is.null(n_feature)){n_feature=2000}
    NBumi_mat <- NBumiConvertData(counts_mat, is.counts=TRUE)
    DANB_fit <- NBumiFitModel(NBumi_mat)
    if(method=="NBDisp"){
      genelist = names(NBumiFeatureSelectionHighVar(DANB_fit)[1:n_feature])
    }else if(method=="NBDrop"){
      genelist = NBumiFeatureSelectionCombinedDrop(DANB_fit, qval.thresh = p_thresh)$Gene
      }
  }else if(method%in%c("HIPPO")){
    hippo_sce=sce
    hippo_sce = hippo(hippo_sce, K=8, feature_method = "zero_inflation",
                      clustering_method = "kmeans")
    k_length=length(hippo_sce@int_metadata$hippo$features)
    genelist = hippo_sce@int_metadata$hippo$features[[k_length]]$gene
  }else if(method%in%c("Townes")){
    if(is.null(n_feature)){n_feature=2000}
    Townes_sce = sce
    Townes_sce = devianceFeatureSelection(sce, assay="counts", nkeep=n_feature, sorted=TRUE)
    genelist = rownames(Townes_sce)
  }else if(method%in%c("SLICER")){
    data = CreateSeuratObject(counts = counts_mat) 
    data = NormalizeData(data, verbose = FALSE)
    mat = as.matrix(data@assays[["RNA"]]@data)
    mat = t(mat)
    genes = select_genes(mat)
    genelist = rownames(counts_mat)[genes]
  }else{
    print("Method is not available in R.")
    genelist = NA
  }
  end_time = Sys.time()
  running_time = end_time - start_time
  result$time = running_time
  result$genelist = genelist
  return(result)
}

Cluster=function(sce, gene_list){
  if(length(gene_list)==1){
    cluster_ident="NA"
    return(cluster_ident)}
  else{
    counts_mat=assay(sce,"counts")
    SeuObj = CreateSeuratObject(counts=counts_mat)
    SeuObj <- NormalizeData(SeuObj,verbose = FALSE)
    all.gene = rownames(SeuObj)
    SeuObj <- ScaleData(SeuObj, features = all.gene, verbose = FALSE)
    SeuObj <- RunPCA(SeuObj, features = gene_list,verbose = FALSE)
    SeuObj <- FindNeighbors(SeuObj, verbose = FALSE)
    SeuObj <- FindClusters(SeuObj, verbose = FALSE)
    cluster_ident=as.numeric(Idents(SeuObj))
    return(cluster_ident)
  }
}

## Get feature gene list
methods=c("disp","mvp","vst","SLICER","scVEGs","scran","scGEAToolbox",
          "Brennecke","BASiCS","NBDisp","NBDrop","sctransform","Townes","HIPPO",
          "GiniClust","M3Drop","EDGE")
GeneList=list()
for(i in 1:length(methods)){
  GeneList[[i]]=FeatureGeneSelection(sce,method = methods[i])
}
names(GeneList)=methods

## Clustering
Cluster_result=list()
for(i in 1:length(methods)){
  Cluster_result[[i]]=Cluster(sce, GeneList[[i]]$genelist)
}
names(Cluster_result)=methods

## Visualization
annotation=sce$phenoid

### Remove the methods that clustering results are not available
Cluster_result$Brennecke=NULL
Cluster_result$BASiCS=NULL
Cluster_result$EDGE=NULL
Cluster_result$M3Drop=NULL

EDB = sum(names(Cluster_result)%in% c("disp","mvp","vst","SLICER","scVEGs","scran","scGEAToolbox"))
GMB = sum(names(Cluster_result)%in% c("Brennecke","BASiCS","NBDisp","NBDrop","sctransform","Townes","HIPPO"))
DF = sum(names(Cluster_result)%in% c("GiniClust","M3Drop","EDGE"))

if(DF!=0){
  Category=factor(c(rep("Empirical-distribution--based methods",EDB),
                    rep("Generative-model--based methods",GMB),
                    "Distribution free methods",DF),
                  levels = c("Empirical-distribution--based methods",
                             "Generative-model--based methods", 
                             "Distribution free methods"))}else{
  Category=factor(c(rep("Empirical-distribution--based methods",EDB),
                    rep("Generative-model--based methods",GMB)),
                  levels = c("Empirical-distribution--based methods",
                             "Generative-model--based methods"))
                             }

AdjustedRandIndex_res=unlist(lapply(Cluster_result, function(x){adjustedRandIndex(x,annotation)}))
AdjustedRandIndex_res=data.frame(AdjustedRandIndex=AdjustedRandIndex_res,methods=as.factor(names(Cluster_result)),Category=Category)
reorder=order(AdjustedRandIndex_res$AdjustedRandIndex)
AdjustedRandIndex_res=AdjustedRandIndex_res[reorder,]
AdjustedRandIndex_res$methods=fct_inorder(AdjustedRandIndex_res$methods)
pdf("ADI_figure.pdf")
ggplot(data = AdjustedRandIndex_res, mapping = aes(x = methods, y = AdjustedRandIndex, fill=Category)) + 
  geom_bar(stat = 'identity', show.legend = FALSE,width=0.9) +
  geom_text(aes(label=round(AdjustedRandIndex,3)), hjust=-0.1, size=6)+
  scale_y_continuous(limits = c(0,1)) +
  scale_fill_brewer(palette = "Dark2") +
  coord_flip()  + theme_bw()+theme(panel.grid = element_blank(),
                                   panel.border = element_blank(),
                                   axis.line.x=element_line(colour="black"),
                                   axis.line.y=element_line(colour="black"),
                                   axis.text=element_text(size=20),
                                   axis.title=element_text(size=20)) +
  
  xlab(NULL) + ylab(NULL)
dev.off()