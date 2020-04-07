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

cli$group = cli[,'Ig clonality:ch1']

## step 1) log transform and adjust to avoide negative inf
df = log2(GSE6136+1)
tmp = names(table(cli$group))
samplesL = split(cli$geo_accession,f = cli$group)

## step 2.1) pre-select the transcrips
## optimation is ignored in this trial due to small sample size
table(cli[,'Ig clonality:ch1'])
# monoclonal oligoclonal  polyclonal
#          7          11           8


## Instead, preselect the top 1% transcripts with the highest relative standard deviation
## in each state compared to transcript in the other states
subcounts = sd_selection(df,samplesL)
sapply(subcounts,dim)
#     monoclonal oligoclonal polyclonal
#[1,]        226         226        226
#[2,]          7          11          8

pre3s = sapply(subcounts,row.names)[,'monoclonal']
save(pre3s, file = 'pre3s.sdrownames_p5.rData')

pre3s = sapply(subcounts,row.names)[,'polyclonal']
save(pre3s, file = 'pre3s.sdrownames_p1p2.rData')

pre3s = do.call(c,lapply(subcounts,row.names))
save(pre3s, file = 'pre3s.sdrownames_all.rData')

subcounts = sd_selection(df,samplesL,cutoff = 0.05)
pre3s = row.names(subcounts[['polyclonal']])
save(pre3s, file = 'pre3s.sdrownames_p1p2_5per.rData')

#pre3s = row.names(subcounts[['monoclonal']])
#save(pre3s, file = 'pre3s.sdrownames_p5_5per.rData')

#subcounts = sd_selection(df,samplesL,cutoff = 0.1)
#pre3s = row.names(subcounts[['monoclonal']])
#save(pre3s, file = 'pre3s.sdrownames_p5_10per.rData')

## step 3.1) Build gene regulatory network based on expresional PCC
igraphL = getNetwork(subcounts,fdr=1)
#monoclonal:218 nodes
#oligoclonal:218 nodes
#polyclonal:220 nodes

## step 3.2) Network partition using random walk
## this step constructs modules (geoups of transcriopts) based on their expression levels
cluster = getCluster_methods(igraphL)
par(mfrow = c(1,3))

## step 4.1-4.3) calculate MCI score for each module
membersL_noweight = getMCI(cluster,subcounts,adjust.size = F)
plotBar_MCI(membersL_noweight,ylim = c(0,4),min = 10,order = c('polyclonal','oligoclonal','monoclonal'))

## step 4.4) identify biomodule per state
maxMCIms <- getMaxMCImember(membersL_noweight[[1]],membersL_noweight[[2]],min =10)

maxMCI = getMaxStats(membersL_noweight[['MCI']],maxMCIms[[1]])
## get biomodule with the highest MCI, CTS
CTS = getCTS(maxMCI, maxMCIms[[2]])
#Length: 105

save(CTS, file = 'CTS.GSE6136_3states.rData')

simMCI = simulationMCI(91,samplesL, df)
plot_MCI_Simulation(maxMCI,simMCI,ylim =c(0,4))
dev.copy2pdf(file = 'box_GSE6136_MCIsim_3stages.pdf')

symbol = read.delim("F:\\ZheZhen\\NB\\DNS\\data\\GPL8321-17396.txt",comment = '#')
CTS_symbol = subset(symbol,ID %in% CTS)$'Gene.Symbol'
load("F:\\ZheZhen\\NB\\DNS\\result\\dnb_symbol.GSE6136_6thtry.rData")

library(gplots)
venn(list(CTS = CTS_symbol, oldCTS = dnb_symbol))

## step 5.1) Ic 
par(mfrow = c(1,1))
Ic = getIc(df,samplesL,CTS)
simIc = simulation_Ic(length(CTS),samplesL,df)
#plot_Ic_Simulation(Ic,simIc,order = c('polyclonal','oligoclonal','monoclonal'))
plot_Ic_Simulation(Ic,simIc,order = c('polyclonal','oligoclonal','monoclonal'),ylim = c(0,2.5))
dev.copy2pdf(file = 'lines_GSE6136_Icsim_3stages.pdf')
#plot_Ic_Simulation(Ic,simIc,order = c('polyclonal','oligoclonal','monoclonal'),ylim = c(0,15))
#dev.copy2pdf(file = 'lines_GSE6136_Icsim_3stages_15ylim.pdf')
