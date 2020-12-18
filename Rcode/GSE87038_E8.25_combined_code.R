# This code ran BioTIP analysis on the published E8.25 developing mesoderming cells, following the guildline given by the two tutorials:
# http://bioconductor.org/books/release/OSCA/pijuan-sala-chimeric-mouse-embryo-10x-genomics.html#ref-pijuansala2019single
# https://bioconductor.org/packages/release/data/experiment/vignettes/MouseGastrulationData/inst/doc/MouseGastrulationData.html
# Last update 12/15/2020
# by Holly Yang xyang2_at_uchicago.edu

## setting working directory
setwd('F:/projects/scRNA/results/GSE87038_gastrulation/')
FigPath <- 'F:/projects/BioTIP/doc/2020_/Applications/ROutputFigs/GSE87038_gastrulation/'
DataPath <- 'F:/projects/BioTIP/doc/2020_/Applications/outPut.Rdata/GSE87038_gastrulation/'

library(MouseGastrulationData)
## load in only E8.25 samples from the atlas to reduce memory consumption when compiling this vignette.
subset(AtlasSampleMetadata, stage=='E8.25')
#     sample stage pool_index seq_batch ncells
# 23     24 E8.25         19         2   6707
# 24     25 E8.25         19         2   7289
# 27     28 E8.25         21         2   4646
sampleID <- subset(AtlasSampleMetadata, stage=='E8.25')$sample

##########################################
## download E8.25 datasets and size-factor normalized one-by-one
##############################################
# pca.corrected  is the output of author's fastMNN correction over the whole atlas. 
library(scater)
sce.list = list()

for(i in 1:length(sampleID))
{
  sce <- EmbryoAtlasData(samples = sampleID[i])
  rownames(sce) <- uniquifyFeatureNames(
    rowData(sce)$ENSEMBL, rowData(sce)$SYMBOL)
  
  ## Quality control on the cells has already been performed by the authors, 
  # We additionally remove cells that are labelled as stripped nuclei or doublets.
  # and check if any cell has high mitochondrial gene content
  table(grepl("^MT-", rownames(sce)))  # FALSE
  drop <- (sce$doublet | sce$stripped)
  sce <- sce[,!drop]
  dim(sce) #29452  5717
  sce.list[[i]] <- sce
  
  #We do log-normalization using the pre-computed size factors in sce.chimera.
  sce.list[[i]] <- logNormCounts(sce.list[[i]])
  
}
names(sce.list) <- sampleID


##########################################
##  merge three E8.25 datasets(10x sampels)
##############################################
#http://bioconductor.org/books/release/OSCA/integrating-datasets.html#setting-up-the-data

library(scran)
dec24 <- modelGeneVar(sce.list[[1]])
dec25 <- modelGeneVar(sce.list[[2]])
dec28 <- modelGeneVar(sce.list[[3]])

# The multiBatchNorm() function recomputes log-normalized expression values 
# after adjusting the size factors for systematic differences in coverage between SingleCellExperiment objects. 
library(batchelor)
sce <- multiBatchNorm(sce.24, sce.25, sce.28)

# Feature selection.
library(scran)
sce.24 <- sce[[1]]
sce.25 <- sce[[2]]
sce.28 <- sce[[3]]
chosen.hvgs <- dec24$bio > 0 & dec25$bio>0 & dec28$bio>0
sum(chosen.hvgs)  # 11313

# author did UMAP analysis based on globally batch-corrected PCA based on all samples
reducedDimNames(sce.24)[2] <- 'UMAP' # rename to call scater::runPCA later 
reducedDimNames(sce.25)[2] <- 'UMAP' # rename to call scater::runPCA later
reducedDimNames(sce.28)[2] <- 'UMAP' # rename to call scater::runPCA later

#################################################
# Diagnosing batch effects
# Synchronizing the metadata for cbinding. 
#################################################
uncorrected <- cbind(sce.list[[1]],sce.list[[2]],sce.list[[3]])
uncorrected$batch =rep(1, ncol(uncorrected))
uncorrected$batch[which(uncorrected$sample==25)] <- 2
uncorrected$batch[which(uncorrected$sample==28)] <- 3
table(uncorrected$batch)
#    1    2    3 
# 5717 6244 3974
uncorrected$batch <- factor(uncorrected$batch)
uncorrected$sample <- factor(uncorrected$sample)

# Using RandomParam() as it is more efficient for file-backed matrices.
library(scater)
set.seed(0010101010)
uncorrected <- runPCA(uncorrected, subset_row=chosen.hvgs,
                      BSPARAM=BiocSingular::RandomParam())

library(scran)
snn.gr <- buildSNNGraph(uncorrected, use.dimred="pca.corrected")
clusters <- igraph::cluster_walktrap(snn.gr)$membership

rm(sce.list)
# save(uncorrected, clusters, 
#      file='uncorrected/E8.25_uncorrected.RData', compress=TRUE) #####################

plotUMAP(uncorrected, colour_by="sample",point_size=0.5) + 
  scale_colour_manual(values = c("green", "blue", "orange"))
dev.copy2pdf(file=paste0(FigPath, 'UMAP_author.corrected_E8.25.pdf'))


# We use the adjusted Rand index (Section 10.6.2) to quantify the agreement between the clusterings before and after batch correction
library(bluster)
pairwiseRand(clusters[uncorrected$sample==24], uncorrected$celltype[uncorrected$sample==24], mode="index")  #0.7004609
pairwiseRand(clusters[uncorrected$sample==25], uncorrected$celltype[uncorrected$sample==25], mode="index")  #0.69
pairwiseRand(clusters[uncorrected$sample==28], uncorrected$celltype[uncorrected$sample==28], mode="index")  #0.6446281

pairwiseRand(clusters[uncorrected$sample==24], uncorrected$cluster.stage[uncorrected$sample==24], mode="index") ## 0.7168464
pairwiseRand(clusters[uncorrected$sample==25], uncorrected$cluster.stage[uncorrected$sample==25], mode="index") ## 0.7166239
pairwiseRand(clusters[uncorrected$sample==28], uncorrected$cluster.stage[uncorrected$sample==28], mode="index") ## 0.7371658

#myplot_TSNE_UMAP(uncorrected, path='uncorrected', filename='uncorrected')
  

###########################################
## Focus on 7240 developing mesoderm cels
###########################################

int <- c('Somitic mesoderm','Intermediate mesoderm','ExE mesoderm','Paraxial mesoderm','Allantois',
         'Pharyngeal mesoderm','Cardiomyocytes', 'Mesenchyme',
         'Haematoendothelial progenitors','Endothelium','Blood progenitors 1','Blood progenitors 2'
)
sce <- uncorrected[,uncorrected$celltype %in% int]
dim(sce) #29452 7240 
rm(uncorrected)

## Feature selection 
set.seed(0010101)
dec.pois <- modelGeneVarByPoisson(sce,  block=sce$batch)
### have positive biological components in all combined E8.25 cells while accounting for batch effects and 
hvg.var <- rownames(dec.pois)[ dec.pois$bio>0]
length(hvg.var)  #[1] 16052
### have positive biological components in each dataset
chosen.hvgs <- rownames(sce)[dec24$bio > 0 & dec25$bio>0 & dec28$bio>0]
length(chosen.hvgs)  # 11313

hvg.var <- intersect(chosen.hvgs, hvg.var)
length(hvg.var) # [1] 10938

sce <- sce[hvg.var,]
dim(sce) # 10938 7240 
reducedDimNames(sce)[2] <- 'umap'  # recover the name as author assigned


##########################################
## cluster cells, using k=10 and the downloaded pca.corrected, regardless of feature selection
##########################################
library(scran)
g.10 <- buildSNNGraph(sce, k=10, use.dimred = 'pca.corrected')
clust.10 <- igraph::cluster_walktrap(g.10)$membership
table(clust.10)

colLabels(sce) <- factor(clust.10)

set.seed(11000)
reducedDim(sce, "force") <- igraph::layout_with_fr(g)
  
gridExtra::grid.arrange(
    plotReducedDim(sce, "umap", colour_by="celltype",text_by='celltype', text_size = 3, add_legend=FALSE),
    plotReducedDim(sce, "umap", colour_by="label",text_by='label', text_size = 3),
    ncol=2)
  
dev.copy2pdf(file=paste0(FigPath,'ReduceDim_cluster.k10.pdf'), height=5)

  
####################################################################
## annotate the cluster 
####################################################################


library(scran)
library(pheatmap)
markers.up <- findMarkers(sce, test="wilcox", # if wilcox test rather than t-test, get AUC rather than lfc
                          groups=sce$label,  #lfc=1,
                          direction="up") #, block=mnn.all$sample)
save(markers.up, file=paste0(DataPath,'markers.up.RData'))

pdf(file=paste0(FigPath,'heatmap_pval_upgene.pdf'))
for(chosen in 1:length(markers.z))
{
  interesting.up <- markers.up[[chosen]]
  best.set <- interesting.up[interesting.up$Top <= 5,]
  AUCs <- getMarkerEffects(best.set, prefix="AUC")
  pheatmap(AUCs, breaks=seq(-5, 5, length.out=101),
           fontsize = 10, main= paste('S',chosen),
           scale='row')
  
}
dev.off()




##########################################
## trajectory analysis
##########################################

## Minimum spanning tree constructed trajectory

library(scater)

## by clusters ------------------
by.cluster <- aggregateAcrossCells(sce, ids=colData(sce)$label)

centroids <- reducedDim(by.cluster, "pca.corrected")
dmat <- dist(centroids)
dmat <- as.matrix(dmat)
set.seed(1000)
g <- igraph::graph.adjacency(dmat, mode = "undirected", weighted = TRUE)
mst <- igraph::minimum.spanning.tree(g)
plot(mst)

pairs <- Matrix::which(mst[] > 0, arr.ind=TRUE)
coords <- reducedDim(by.cluster, "umap")
group <- rep(seq_len(nrow(pairs)), 2)
stuff <- data.frame(rbind(coords[pairs[,1],], coords[pairs[,2],]), group)

plotReducedDim(sce, 'umap', colour_by="label", 
               text_by="label", text_size = 8) + 
  geom_line(data=stuff, mapping=aes(x=x, y=y, group=group))

dev.copy2pdf(file=paste0(FigPath,"trajectory_libsf_umap.pdf"))




#####################################################
#######  BioTIP analysis    ######################
#####################################################
# Pre-selection Transcript
cut.preselect = 0.1
cut.fdr = 0.2   
cut.minsize = 60  

samplesL <- split(rownames(colData(sce)),f = colData(sce)$label)
lengths(samplesL)

##libsf.var is the top 4000 HVGs per cell from whole E8.25 dataset
library(scran)

libsf.var <- getTopHVGs(dec.pois, n=4000)
table(libsf.var %in% rownames(sce))
# FALSE  TRUE
#  927  3073
libsf.var <- intersect(libsf.var, rownames(sce))
length(libsf.var) # 3073

dat <- sce[libsf.var,]

#logcounts normalized by library size
logmat <- as.matrix(logcounts(dat))
dim(logmat) # [1] 3073 7240

rm(dat)

#this optimizes sd selection
set.seed(2020)
testres <- optimize.sd_selection(logmat[,unlist(samplesL)], samplesL, B=100, cutoff=cut.preselect,
                                 times=.75, percent=0.8)
#save(testres, file=paste0(DataPath,"optimized_test_sd_selection_E8.25.RData"), compress=TRUE)


# Network Partition
igraphL <- getNetwork(testres, fdr = cut.fdr)
# 1:282 nodes
# 2:288 nodes
# 3:281 nodes
# 4:268 nodes
# 5:286 nodes
# 6:294 nodes
# 7:229 nodes
# 8:273 nodes
# 9:268 nodes
# 10:284 nodes
# 11:173 nodes
# 12:249 nodes
# 13:279 nodes
# 14:205 nodes
# 15:255 nodes
# 16:247 nodes
# 17:222 nodes
# 18:64 nodes
# 19:113 nodes

cluster <- getCluster_methods(igraphL)

# 2.1) Identifying CTS, new MCI score #########
membersL <- getMCI(cluster,testres, adjust.size = FALSE, fun='BioTIP')
names(membersL)
#[1] "members" "MCI"     "sd"      "PCC"     "PCCo"  
save(membersL, file=paste0(DataPath,"membersL.RData"))

pdf(file=paste0(FigPath,"MCIBar_E8.25_", cut.preselect,
                "_fdr",cut.fdr,"_minsize",cut.minsize,".pdf"),
    width=40, height=5)
plotBar_MCI(membersL, ylim=c(0,10), minsize = cut.minsize)
dev.off()
#the numbers in the parenthesis: they are total number of clusters, no control of the cell (sample) size per  cluster


# get the statistics using the MCI system
n.states = 4
topMCI = getTopMCI(membersL[["members"]], membersL[["MCI"]], membersL[["MCI"]], min=cut.minsize, n=n.states)
names(topMCI)
#[1] "15" "6"  "16" "13" 

# get the state ID and MCI statistics for the two leading MCI scores per state
maxMCIms <- getMaxMCImember(membersL[["members"]],
                            membersL[["MCI"]],min =cut.minsize, n=2)
names(maxMCIms)
#[1] "idx"             "members"         "2topest.members"

maxMCI = getMaxStats(membersL[['MCI']], maxMCIms[['idx']])
unlist(maxMCI)
# 1.527133 1.438823 1.564383 1.618711 

names(maxMCIms[["members"]][names(topMCI)])
#[1] "15" "6"  "16" "13" 

CTS.Lib = getCTS(maxMCI[names(topMCI)], maxMCIms[["members"]][names(topMCI)])
# Length: 67
# Length: 90
# Length: 79
# Length: 60


# translate to ID, already SYMBOL
CTS.Lib.Symbol <- CTS.Lib
for(i in 1:length(CTS.Lib))
{
  CTS.Lib[[i]] <- rowData(sce)[CTS.Lib.Symbol[[i]],]$ENSEMBL
}
'Etv2' %in% unlist(CTS.Lib.Symbol) #[1] TRUE for cut.fdr=0.2
save(CTS.Lib,  CTS.Lib.Symbol, maxMCI, file=paste0(DataPath,"CTS_Lib_E8.25.RData"))
rm(sce)

# estimate significance NOT repeat (takes a while), no need to rerun
dim(logmat[,unlist(samplesL)])
#[1]  3073 7240
# M is precalculated correlation matrix, will be reused in the downstream simulation analysis
M <- cor.shrink(logmat[,unlist(samplesL)], Y = NULL, MARGIN = 1, shrink = TRUE)
#save(M, file=paste0(outputpath,"CTS_ShrinkM_E8.25.RData"), compress=TRUE) # save the file as its calculation runs a while
#load(file=paste0(outputpath,"CTS_ShrinkM_E8.25.RData"))
dim(M)  #3073  3073

C = 1000
simuMCI = list()
set.seed(2020)
for (i in 1:length(CTS.Lib)){
  n <- length(CTS.Lib[[i]])
  simuMCI[[i]] <- simulationMCI(n, samplesL, logmat,  B=C, fun="BioTIP",M=M)
}
names(simuMCI) = names(CTS.Lib)

#save(simuMCI, file=paste0(outputpath,"SimuMCI_",C,"_E8.25", cut.preselect, "_fdr",cut.fdr,"_minsize",cut.minsize,".RData"))

par(mfrow=c(2,2))  # plot of 17
for (i in 1:4){
  plot_MCI_Simulation(maxMCI[i], simuMCI[[i]], las=2, ylim=c(0, 3.6),
                      main=paste("S", names(maxMCI)[i], ";", length(CTS.Lib[[i]]), "genes",
                                 "\n","vs. ",C, "times of gene-permutation"),
                      which2point=names(maxMCI)[i])
}
dev.copy2pdf(file=paste0(FigPath,"plot_MCI_Sim_RandomGene_E8.25_top4.pdf"))


######## 3)  Finding Tipping Point and verify using IC score  #################
##################################################
## test PCC_s -> no shrinkage; PCC_g -> zero, permutation both genes and sampels ---------------------------------
C= 1000
SimResults_s <- list()

#set.seed(2020)
set.seed(101010)
for(i in 1:length(CTS.Lib)){
  ######  BioTIP score, shuffling both genes and samples ##############
  CTS <- CTS.Lib.Symbol[[i]]
  n <- length(CTS)
  
  SimResults_s[[i]] <- matrix(nrow=length(samplesL), ncol=C)
  rownames(SimResults_s[[i]]) = names(samplesL)
  for(j in 1:length(samplesL)) {
    ns <- length(samplesL[[j]])  # of each state
    CTS.sim <- sample(rownames(logmat), n) # !!!!
    SimResults_s[[i]][j,] <- simulation_Ic_sample(logmat[,unlist(samplesL)], ns, #Ic=BioTIP_score[i],
                                                  genes= CTS.sim, B=C,
                                                  fun="BioTIP", shrink=TRUE, PCC_sample.target = 'none') 
  }
}

names(SimResults_s) <- names(CTS.Lib)

# save( SimResults_s, 
#       file=paste0(outputpath,"LibSF_IC_sim_BioTIP_E8.25.PermutateBoth.RData"), compress=TRUE)
#load('uncorrected/E8.25_mesoderm/BioTIP_top0.1FDR0.2/PCCsNoShrink/LibSF_IC_sim_BioTIP_E8.25.PermutateBoth.RData')

pdf(file=paste0(FigPath,"IC_SS_SimresultBoth_E8.25.pdf"),width=10, height=10) 
ylim=0.6
par(mfrow=c(2,2)) # 16 plots
for(i in 1:4){
  n <- length(CTS.Lib[[i]])
  interesting = which(names(samplesL) == names(BioTIP_scores[i]))
  plot_Ic_Simulation(BioTIP_scores[[i]], SimResults_s[[i]], las = 2, ylab="Ic*", ylim=c(0, ylim),
                     main=paste(names(CTS.Lib)[i],"_",n, "genes", "\n","vs. ",C, "sample-permutations"),
                     fun="matplot", which2point= interesting)
  plot_SS_Simulation(BioTIP_scores[[i]], SimResults_s[[i]], 
                     main = "Delta Ic*, gene & label shuffling", 
                     ylab=NULL,
                     xlim=range(c(BioTIP_scores[[i]][names(BioTIP_scores)[i]],SimResults_s[[i]])))
}

dev.off()


#############################################
### 3) motif analysis using Homer2, do not repeat #######
### see bioTIP_librarysize_motif_xy.R #########






sessionInfo( )
# R version 4.0.2 (2020-06-22)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 18363)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
# [3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.1252    
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods  
# [9] base     
# 
# other attached packages:
#   [1] bluster_1.0.0               rstudioapi_0.13            
# [3] MouseGastrulationData_1.4.0 scater_1.18.3              
# [5] ggplot2_3.3.2               scran_1.18.1               
# [7] SingleCellExperiment_1.11.9 SummarizedExperiment_1.19.9
# [9] Biobase_2.49.1              GenomicRanges_1.41.6       
# [11] GenomeInfoDb_1.25.11        IRanges_2.23.10            
# [13] S4Vectors_0.27.14           BiocGenerics_0.35.4        
# [15] MatrixGenerics_1.1.8        matrixStats_0.57.0         
# 
# loaded via a namespace (and not attached):
#   [1] bitops_1.0-6                  bit64_4.0.5                  
# [3] RColorBrewer_1.1-2            RcppAnnoy_0.0.16             
# [5] httr_1.4.2                    tools_4.0.2                  
# [7] R6_2.5.0                      irlba_2.3.3                  
# [9] vipor_0.4.5                   uwot_0.1.8                   
# [11] DBI_1.1.0                     colorspace_1.4-1             
# [13] withr_2.3.0                   tidyselect_1.1.0             
# [15] gridExtra_2.3                 bit_4.0.4                    
# [17] curl_4.3                      compiler_4.0.2               
# [19] BiocNeighbors_1.8.1           DelayedArray_0.15.17         
# [21] labeling_0.4.2                scales_1.1.1                 
# [23] rappdirs_0.3.1                digest_0.6.27                
# [25] XVector_0.29.3                htmltools_0.5.0              
# [27] pkgconfig_2.0.3               sparseMatrixStats_1.2.0      
# [29] fastmap_1.0.1                 dbplyr_2.0.0                 
# [31] limma_3.46.0                  rlang_0.4.8                  
# [33] RSQLite_2.2.1                 shiny_1.5.0                  
# [35] DelayedMatrixStats_1.12.0     farver_2.0.3                 
# [37] generics_0.1.0                BiocParallel_1.23.3          
# [39] dplyr_1.0.2                   RCurl_1.98-1.2               
# [41] magrittr_1.5                  BiocSingular_1.6.0           
# [43] GenomeInfoDbData_1.2.4        scuttle_1.0.0                
# [45] Matrix_1.2-18                 Rcpp_1.0.5                   
# [47] ggbeeswarm_0.6.0              munsell_0.5.0                
# [49] viridis_0.5.1                 lifecycle_0.2.0              
# [51] yaml_2.2.1                    edgeR_3.32.0                 
# [53] zlibbioc_1.35.0               Rtsne_0.15                   
# [55] BiocFileCache_1.14.0          AnnotationHub_2.22.0         
# [57] grid_4.0.2                    blob_1.2.1                   
# [59] promises_1.1.1                dqrng_0.2.1                  
# [61] ExperimentHub_1.16.0          crayon_1.3.4                 
# [63] lattice_0.20-41               cowplot_1.1.0                
# [65] beachmat_2.6.1                locfit_1.5-9.4               
# [67] pillar_1.4.6                  igraph_1.2.6                 
# [69] codetools_0.2-18              glue_1.4.2                   
# [71] BiocVersion_3.12.0            BiocManager_1.30.10          
# [73] httpuv_1.5.4                  vctrs_0.3.4                  
# [75] gtable_0.3.0                  purrr_0.3.4                  
# [77] assertthat_0.2.1              mime_0.9                     
# [79] rsvd_1.0.3                    xtable_1.8-4                 
# [81] RSpectra_0.16-0               later_1.1.0.1                
# [83] viridisLite_0.3.0             tibble_3.0.4                 
# [85] AnnotationDbi_1.52.0          beeswarm_0.2.3               
# [87] memoise_1.1.0                 statmod_1.4.35               
# [89] ellipsis_0.3.1                interactiveDisplayBase_1.28.0

