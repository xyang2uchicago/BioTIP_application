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

cli$group = ifelse(cli$'cell-type:ch1' %in% c('resting','activated'),'stage1',
ifelse(cli$'cell-type:ch1' == 'lymphoma_marginal','stage2',
ifelse(cli$'cell-type:ch1' == 'lymphoma_transitional','stage3','stage4')))

## rename the clinially meaningful states to enable computationally in a expected order
#cli$group = cli$'cell-type:ch1'
#cli$group = ifelse(cli$'cell-type:ch1' %in% c('resting','activated'),'stage1',cli$group)

#table(cli$group)
#  lymphoma_aggressive     lymphoma_marginal lymphoma_transitional                stage1
#                    7                     6                     5                     8

## step 1) log transform and adjust to avoide negative inf
df = log2(GSE6136+1)
tmp = names(table(cli$group))
samplesL = split(cli$geo_accession,f = cli$group)

## step 2.1) pre-select the transcrips
## optimation is ignored in this trial due to small sample size
table(cli$group)
#stage1 stage2 stage3 stage4
#     8      6      5      7

## Instead, preselect the top 1% transcripts with the highest relative standard deviation
## in each state compared to transcript in the other states
subcounts = sd_selection(df,samplesL)
sapply(subcounts,dim)
#     stage1 stage2 stage3 stage4
#[1,]    226    226    226    226
#[2,]      8      6      5      7

pre4s = sapply(subcounts,row.names)[,1]
save(pre4s, file = 'pre4s.sdrownames_state1.rData')

#subcounts = sd_selection(df,samplesL,cutoff = 0.05)
#pre4s = row.names(subcounts[[1]])
#save(pre4s, file = 'pre4s.sdrownames_state1_5per.rData')

#subcounts = sd_selection(df,samplesL,cutoff = 0.1)
#pre4s = row.names(subcounts[[1]])
#save(pre4s, file = 'pre4s.sdrownames_state1_10per.rData')

pre4s = do.call(c, lapply(subcounts,row.names))
save(pre4s, file = 'pre4s.sdrownames_all.rData')

#subcounts = sd_selection(df,samplesL)
pre4s = row.names(subcounts[['stage4']])
save(pre4s, file = 'pre4s.sdrownames_p5.rData')

#subcounts = sd_selection(df,samplesL,cutoff = 0.05)
#pre4s = row.names(subcounts[['stage4']])
#save(pre4s, file = 'pre4s.sdrownames_p5_5per.rData')

#subcounts = sd_selection(df,samplesL,cutoff = 0.1)
#pre4s = row.names(subcounts[['stage4']])
#save(pre4s, file = 'pre4s.sdrownames_p5_10per.rData')


## step 3.1) Build gene regulatory network based on expresional PCC
igraphL = getNetwork(subcounts,fdr=1)
#stage1:220 nodes
#stage2:0 nodes
#stage3:181 nodes
#stage4:218 nodes

## step 3.2) Network partition using random walk
## this step constructs modules (geoups of transcriopts) based on their expression levels
cluster = getCluster_methods(igraphL)

## step 4.1-4.3) calculate MCI score for each module
membersL_noweight = getMCI(cluster,subcounts,adjust.size = F)
par(mfrow = c(1,4))
plotBar_MCI(membersL_noweight,ylim = c(0,4),min = 10)

## step 4.4) identify biomodule per state
maxMCIms <- getMaxMCImember(membersL_noweight[[1]],membersL_noweight[[2]],min =10)

maxMCI = getMaxStats(membersL_noweight[['MCI']],maxMCIms[[1]])

## get biomodule with the highest MCI, CTS
CTS = getCTS(maxMCI, maxMCIms[[2]])
#Length: 105
save(CTS, file = 'CTS.GSE6136_4states.rData')

simMCI = simulationMCI(105,samplesL, df)
plot_MCI_Simulation(maxMCI,simMCI,ylim =c(0,4))
dev.copy2pdf(file = 'box_GSE6136_MCIsim_4stages.pdf')

symbol = read.delim("F:\\ZheZhen\\NB\\DNS\\data\\GPL8321-17396.txt",comment = '#')
CTS_symbol = subset(symbol,ID %in% CTS)$'Gene.Symbol'
load("F:\\ZheZhen\\NB\\DNS\\result\\dnb_symbol.GSE6136_6thtry.rData")

library(gplots)
venn(list(CTS = CTS_symbol, oldCTS = dnb_symbol))

## step 5.1) Ic
Ic = getIc(df,samplesL,CTS)
simIc = simulation_Ic(length(CTS),samplesL,df)
#plot_Ic_Simulation(Ic,simIc)
#dev.copy2pdf(file = 'lines_GSE6136_Icsim_4stages.pdf')

plot_Ic_Simulation(Ic,simIc,ylim = c(0,2.5))
dev.copy2pdf(file = 'lines_GSE6136_Icsim_4stages.pdf')
#plot_Ic_Simulation(Ic,simIc,ylim = c(0,15))
#dev.copy2pdf(file = 'lines_GSE6136_Icsim_4stages_15ylim.pdf')

############# 




