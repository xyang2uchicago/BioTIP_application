# BioTIP_Application: a collection of transcriptomic data analysis using the package BioTIP
## An Overview of BioTIP
In short, BioTIP is an R-package designated for the characterization of biological tipping points. A more detailed overview - (1) what are tipping points, (2) how BioTIP improves on existing methods, (3) what datasets can BioTIP be applied to, and (4) how to install BioTIP - can be found at [BioTIP](https://github.com/xyang2uchicago/BioTIP).

<a name='WORKFLOW'></a>
## BioTIP workflow  
The general workflow of BioTIP is detailed in the following chart: 

   <img src="https://github.com/xyang2uchicago/BioTIP_application/blob/master/FigS1_BioTIP_pipeline_detailed_v7.jpg"
     alt="BioTIP workflow"
     width="600" type="application/pdf"/>

Examples of BioTIP applied to various datasets can be found in the folder [Rcode](https://github.com/xyang2uchicago/BioTIP_application/tree/master/Rcode). 

We here present two cases in the following sections, one with [a single-cell RNA dataset](https://github.com/xyang2uchicago/BioTIP_application/blob/master/README.md#example-with-single-cell-rna-dataset-bloodnet-nestorowa-2016) and the other with [a bulk RNA dataset](https://github.com/xyang2uchicago/BioTIP_application/blob/master/README.md#example-with-bulk-rna-dataset-gse6136).

<a name='scRNA'></a>
## Example with Single-Cell RNA Dataset: BloodNet (Nestorowa 2016)
The following chart shows how BioTIP performs in the single-cell downstream analysis and how BioTIP differs from the existing methods - instead of discovering marker genes from the differential-expression analysis between stable states, BioTIP introduces the tipping-point analysis that focuses on unstable transition states.

<img src="https://github.com/xyang2uchicago/BioTIP_application/blob/master/Fig1_scRNARNAseq_pipeline_2021_xy.jpg"
     alt="BioTIP workflow"
     width="600" type="application/pdf"/> 

We will perform downstream analysis of the existing the BloodNet dataset, the hematopoietic stem and progenitor cells (HSPCs) from ten female mice. The detail of the experiment can be found in details on our [Vignette](https://bioconductor.org/packages/release/bioc/vignettes/BioTIP/inst/doc/BioTIP.html#Data%20Preprocessing%20with%20Trajectory%20Building%20Tools). 

<a name='bulkRNA'></a>
## Example with Bulk RNA dataset: GSE6136
As indicated in the [BioTIP workflow chart](https://github.com/xyang2uchicago/BioTIP_application/blob/master/README.md#biotip-workflow), the analysis can be divided into the following two sections:
1. [Data preprocessing before running BioTIP](https://github.com/xyang2uchicago/BioTIP_application/blob/master/README.md#data-preprocessing-before-running-biotip),
2. [Standard Identification in 6 steps](https://github.com/xyang2uchicago/BioTIP_application/blob/master/README.md#standard-identification-in-6-steps), detailed in the [BioTIP workflow chart](https://github.com/xyang2uchicago/BioTIP_application/blob/master/README.md#biotip-workflow).
###
Data preprocessing before running BioTIP 
------

<a name="Data Preprocessing"></a>
 __Data Preprocessing__

An existing dataset, GSE6136, is used to demonstrate how our functions are applied. Samples were collected from transgenic mouse lymphomas and divided into five groups based on clinical presentation, pathology, and flow cytometry ([Lenburg 2007](https://www.ncbi.nlm.nih.gov/pubmed/?term=17166848)), thus belonging to cross-sectional profiles. Noticing these five group represent a control stage similar to non-transgenic B cells and four major periods of B-cell lymphomagenesis, Dr. Chen and coauthors used the DNB method to identify the pre-disease state exists around the normal activated period (P2), i.e., the system transitions to the disease state around the transitional lymphoma period (Figure S12 in publication ([Chen 2012)](https://www.ncbi.nlm.nih.gov/pubmed/22461973)). 
Start by installing the package 'BioTIP' and other dependent packages such as stringr, psych, and igraph if not. Below are some examples.

```{r}
# load package
library(BioTIP)
# first row of the GEO-downloaded matrix should be column-name.
data(GSE6136_matrix)
GSE6136 = GSE6136_matrix[,-1]  
rownames(GSE6136) <- GSE6136_matrix[,1]
dim(GSE6136)               #[1] 22690 rows and 26 columns
rm(GSE6136_matrix)
```

  The summary of GEO-downloaded clinical information GSE6136_cli is shown below. 

```{r}
#requires library(stringr)
library(BioTIP)
data(GSE6136_cli)

#dim(GSE6136_cli) #check dimension
cli = t(GSE6136_cli)

library(stringr)
colnames(cli) = str_split_fixed(cli[1,],'_',2)[,2]
cli = cli[-1,]
cli = data.frame(cli)
cli[,"cell-type:ch1"] = str_split_fixed(cli$characteristics_ch1.1,": ",2)[,2]
cli[,"Ig clonality:ch1"] = str_split_fixed(cli$characteristics_ch1.3,": ",2)[,2]

colnames(cli)[colnames(cli) == "cell-type:ch1"] = "group"
cli$Row.names = cli[,1]
head(cli[,1:3])
```

  We normalized the expression of genes using log2() scale. This normalization 
  will ensure a more accurate comparison of the variance between the expression
  groups (clusters).

```{r}
df <- log2(GSE6136+1)
 df[1:3,1:5]
#           GSM142398 GSM142399 GSM142400 GSM142401 GSM142402
#1415670_at  10.21917  10.29048  9.987548  10.07682  9.827343
#1415671_at  10.90358  11.15993 10.776186  10.93000 11.468268
#1415672_at  11.11530  10.89209 11.091303  11.04029 11.109504

```

The sample assembles are previously published.
```
tmp <- names(table(cli$group))
cli$group = factor(cli$group, levels=c('resting','activated','lymphoma_marginal','lymphoma_transitional','lymphoma_aggressive'))
samplesL <- split(cli[,1],f = cli$group)
names(samplesL)
#[1] "resting"               "activated"             "lymphoma_marginal"    
#[4] "lymphoma_transitional" "lymphoma_aggressive"  
```

[Go to Top](https://github.com/xyang2uchicago/BioTIP_application/blob/master/README.md#example-with-bulk-rna-dataset-gse6136) 


Standard Identification in 6 steps 
------
<a name="Finding Tipping Point"></a>
__S1: Finding Tipping Point__

The first step is to calculate an Index of Critical transition (Ic score) of the dataset. One can use the function 'getIc' to calculate the Ic score based on randomly selected genes (20-1000). 'GetIc' has a key parameter fun which gives the options to call the existing Ic score or the proposed new Ic* score.
  
```{r}
RandomG <- sample(rownames(df), 100)
IC.new <- getIc(df, samplesL, RandomG, fun='BioTIP')
plotIc(IC.new, las = 2)
IC.old <- getIc(df, samplesL, RandomG, fun='cor')
plotIc(IC.old, las = 2)
```

To estimate the distribution of random Ic scores, one can call the function simulation_Ic. 

```{r,warning=FALSE}
simuIC <- simulation_Ic(length(RandomG),samplesL,df, fun='BioTIP')
par(mar = c(10,5,0,2))
plot_Ic_Simulation(IC,simuIC,las = 2)
```

<a name="Acknowledgements"></a>

<a name="Pre-selection Transcript"></a>
 __S2: Pre-selection Transcript__

  Once pre-processed, we can analyze the clinical-stage ensembles. Here, we
  demonstrate a pre-selection of 226 (1%) genes for each state.

```{r, warning=FALSE}
test <- sd_selection(df, samplesL,0.01)
class(test)  #[1] "list"
lapply(test, dim)
#$resting
#[1] 226   5
#
#$activated
#[1] 226   3
#
#$lymphoma_marginal
#[1] 226   6
#
#$lymphoma_transitional
#[1] 226   5
#
#$lymphoma_aggressive
#[1] 226   7
```

<a name="Network Partition"></a>
 __S3: Network Partition__

  A graphical presentation of genes of interest can be achieved using the following functions. The `getNetwork` function will obtain an igraph object based on Pearson Correlation of `test`. This `igraphL` object is then run using the `getCluster_methods` function classify nodes.

```{r,echo=TRUE, warning=FALSE}
library(igraph)
library(psych)
igraphL <- getNetwork(test, fdr = 1)

cluster <- getCluster_methods(igraphL)
```
```{r,echo=TRUE, warning=FALSE}
names(cluster)
class(cluster)  #[1] "list"
class(cluster[[1]])  #[1] "communities"
```

<a name="Identifying Dynamic Network Biomodule"></a>
 __S4: Identifying Dynamic Network Biomodule__

Here, ‘module’ refers to a cluster of network nodes (e.g. transcripts) highly linked (e.g. by correlation). “Biomodule” refers to the module resenting significantly highest score called “Module-Criticality Index (MCI)” per state.

  The following step shows a graph of classified clustered samples for five different states. MCI score is calculated for each module using the `getMCI`
  function. The `getMaxMCImember` function will obtain a list of modules with the highest
  MCI at each stage. Use `"head(maxCIms)"` to view the MCI scores calculated. Use
  `plotMaxMCI` function to view the line linking the highest MCI score at each state.

```{r,echo=TRUE, warning=FALSE}
membersL_noweight <- getMCI(cluster,test,adjust.size = FALSE)
plotBar_MCI(membersL_noweight,ylim = c(0,6))
```
```{r,echo=TRUE, warning=FALSE}
maxMCIms <- getMaxMCImember(membersL_noweight[[1]],membersL_noweight[[2]],min =10)
names(maxMCIms)
names(maxMCIms[[1]])
names(maxMCIms[[2]])
```
```{r,echo=TRUE, warning=FALSE}
head(maxMCIms[['idx']])
head(maxMCIms[['members']][['lymphoma_aggressive']])
```


  To get the selected statistics of biomodule (i.e., the gene module that presents the highest MCI score) of each state, please run the following 

```{r}
biomodules = getMaxStats(membersL_noweight[['members']],maxMCIms[[1]])
maxMCI = getMaxStats(membersL_noweight[['MCI']],maxMCIms[[1]])
head(maxMCI)
```
```{r}
maxSD = getMaxStats(membersL_noweight[['sd']],maxMCIms[[1]])
head(maxSD)
```

  To get the biomodule with the highest MCI score among all states, as we call it CTS (Critical Transition Signals), please run the following.
  
```{r}
CTS = getCTS(maxMCI, maxMCIms[[2]])
```

 Run the following to visualize the trendence of every state represented by the cluster with the highest MCI scores.
```{r,echo=TRUE, warning=FALSE}
par(mar = c(10,5,0,2))
plotMaxMCI(maxMCIms,membersL_noweight[[2]],states = names(samplesL),las = 2)
```

  We then perform simulation for MCI scores based on identified signature size
  (length(CTS) ) using the `simulationMCI` function.Use `plot_MCI_simulation`
  function to visualize the result. This step usually takes 20-30mins, so here to save time, we picked a small number 3 as the length of the CTS.
  
```{r,echo=TRUE, warning=FALSE}
simuMCI <- simulationMCI(3,samplesL,df)
plot_MCI_Simulation(maxMCI,simuMCI)
```

<a name="Evaluate Tipping Point"></a>
 __S5: Evaluate Tipping Point__

The next step is to calculate an Index of Critical transition (Ic score) of the dataset. First, use the getIc function to calculate the Ic score based on the biomodule previously identified. We use the function plotIc to draw a line plot of the Ic score. 
  
```{r}
IC <- getIc(df,samplesL,CTS)
par(mar = c(10,5,0,2))
plotIc(IC,las = 2)
```

Then use the two functions to evaluate two types of empirical significance,
respectively. The first function simulation_Ic calculates random Ic-scores by
shuffling features (transcripts). Showing in the plot is Ic scores of the
identification (red) against its corresponding size-controlled random scores
(grey).

```{r,warning=FALSE}
simuIC <- simulation_Ic(length(CTS),samplesL,df)
par(mar = c(10,5,0,2))
plot_Ic_Simulation(IC, simuIC,las = 2)
```
  
Importantly, the function plot_simulation_sample calculates random Ic-scores by shuffling
samples and visualizes the results. Showing in the plot is observed Ic-score
(red vertical line) comparing to the density of random scores (grey), at the
tipping point identified.

```{r}
RandomIc = simulation_Ic_sample(df, sampleNo=length(samplesL[['lymphoma_aggressive']]), 
                                 Ic=IC[['lymphoma_aggressive']], genes=CTS, B=100)
plot_Ic_Simulation(IC, RandomIc, las = 2)                                 
```

<a name="Evaluate Delta Score"></a>
 __S6: Evaluate Delta Score__

The final step is to calculate the Delta Score, which evaluates the significance of our identified CTS. 
```{r}
plot_SS_Simulation(Ic,  Ic_simuG, 
                   main = paste("Delta of Ic*", length(CTS),"genes"), 
                   ylab='BioTIP 35 genes') # [1] 0 

plot_Ic_Simulation(Ic,  Ic_simuG,  las = 0,  ylim =  c(0,1.5),  
                   order = NULL,  main = 'BioTIP 35 genes',  
                   ylab = "Ic*",  fun = 'boxplot', 
                   which2point = NULL) 

plot_SS_Simulation(Ic,  Ic_simuG, 
                   main = paste("Delta of Ic*", length(CTS),"genes"), 
                   ylab='BioTIP 35 genes')
```

[Go to Top](https://github.com/xyang2uchicago/BioTIP_application/blob/master/README.md#example-with-bulk-rna-dataset-gse6136) 





 




