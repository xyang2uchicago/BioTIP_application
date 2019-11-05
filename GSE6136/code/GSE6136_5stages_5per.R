library(stringr)
library(psych)
library(igraph)

####################### load in the downloaded expressinal matrix of GSE6136 #####################
 setwd("F:\\ZheZhen\\NB\\DNS\\result")
 GSE6136 = read.table("F:\\ZheZhen\\NB\\DNS\\data\\GSE6136_matrix.txt",head = T)
 dim(GSE6136) #[1] 22690    27
 row.names(GSE6136) = GSE6136$ID_REF
 GSE6136 = GSE6136[,-1]
 dim(GSE6136) #[1] 22690    26

 ## read the published sample states
 cli = read.delim("F:\\ZheZhen\\NB\\DNS\\data\\GSE6136_cli.txt",head = F)
 cli = t(cli)
 library(stringr)
 colnames(cli) = str_split_fixed(cli[1,],'_',2)[,2]
 cli = cli[-1,]
 cli = data.frame(cli)
 cli[,'cell-type:ch1'] = str_split_fixed(cli$characteristics_ch1.1,': ',2)[,2]
 cli[,'Ig clonality:ch1'] = str_split_fixed(cli$characteristics_ch1.3,': ',2)[,2]

 table(cli[,'cell-type:ch1'])
#            activated (P2)  lymphoma_aggressive (P5)    lymphoma_marginal (P3)
#                    3                     7                     6
#lymphoma_transitional (P4)              resting (P1)
#                   5                     5
#table(cli[,'genotype:ch1'])
# E-mu-BRD2  wildtype
#       21         5
table(cli[,'Ig clonality:ch1'])
# monoclonal oligoclonal  polyclonal
#          7          11           8

####################### BioTIP application  #####################
## as long as the paackage BioTIP is availabel from Bioconductor, you will use
# library(BioTIP)
# instead of the following one command
source('F:\\ZheZhen\\NB\\DNS\\BioTIP\\R\\BioTIP.R')

cli$group = cli$'cell-type:ch1'
df = log2(GSE6136+1)
tmp = names(table(cli$group))
samplesL = split(cli$geo_accession,f = cli$group)

## step 2.1) pre-select the transcrips
## optimation is ignored in this trial due to small sample size
table(cli$group)
#            activated   lymphoma_aggressive     lymphoma_marginal lymphoma_transitional               resting
#                    3                     7                     6                     5                     5

## Instead, preselect the top 1% transcripts with the highest relative standard deviation
## in each state compared to transcript in the other states
subcounts = sd_selection(df,samplesL,cutoff= 0.05)
sapply(subcounts,dim)
#     activated lymphoma_aggressive lymphoma_marginal lymphoma_transitional resting
#[1,]      1134                1134              1134                  1134    1134
#[2,]         3                   7                 6                     5       5

pre5s_5per = row.names(subcounts[['activated']])
save(pre5s_5per, file = 'pre5s_5per_preselect_GSM142418.rData')

## step 3.1) Build gene regulatory network based on expresional PCC
igraphL = getNetwork(subcounts,fdr=1)
#activated:1134 nodes
#lymphoma_aggressive:1134 nodes
#lymphoma_marginal:1134 nodes
#lymphoma_transitional:1134 nodes
#resting:1134 nodes

## step 3.2) Network partition using random walk
## this step constructs modules (geoups of transcriopts) based on their expression levels
cluster = getCluster_methods(igraphL)

## step 4.1-4.3) calculate MCI score for each module
membersL_noweight = getMCI(cluster,subcounts,adjust.size = F)
par(mfrow = c(1,4))
plotBar_MCI(membersL_noweight,ylim = c(0,4),min = 30)
dev.copy2pdf(file = 'bar_GSE6136_5states_5per.pdf')

## step 4.4) identify biomodule per state
maxMCIms <- getMaxMCImember(membersL_noweight[[1]],membersL_noweight[[2]],min =30)

maxMCI = getMaxStats(membersL_noweight[['MCI']],maxMCIms[[1]])

## get biomodule with the highest MCI, CTS
CTS = getCTS(maxMCI, maxMCIms[[2]])
#Length: 141
save(CTS, file = 'CTS.GSE6136_5states_5per.rData')

par(mfrow = c(2,1))
simMCI = simulationMCI(105,samplesL, df)
plot_MCI_Simulation(maxMCI,simMCI,ylim =c(0,5),las = 2,order = c('resting','activated','lymphoma_marginal','lymphoma_transitional','lymphoma_aggressive'))

Ic = getIc(df,samplesL,CTS)
simIc = simulation_Ic(length(CTS),samplesL,df)
plot_Ic_Simulation(Ic,simIc,ylim = c(0,6),order = c('resting','activated','lymphoma_marginal','lymphoma_transitional','lymphoma_aggressive'))

dev.copy2pdf(file = 'MCI_Ic_GSE6136_5states_5per.pdf')
