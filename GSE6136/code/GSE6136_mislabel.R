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
 #activated   lymphoma_aggressive     lymphoma_marginal lymphoma_transitional               resting
 #3                     7                     6                     5                     5

 ### F:\ZheZhen\NB\DNS\result\mislabel\PCA_vsn_category_wname.pdf 1 sample from lymphoma_marginal to activated
 cli[which(cli$geo_accession == 'GSM142409'),'cell-type:ch1'] = 'activated'
 table(cli[,'cell-type:ch1'])
 #activated   lymphoma_aggressive     lymphoma_marginal lymphoma_transitional               resting
 #4                     7                     5                     5                     5


source('F:\\ZheZhen\\NB\\DNS\\BioTIP\\R\\BioTIP.R')

cli$group = cli[,'cell-type:ch1']

## step 1) log transform and adjust to avoide negative inf
df = log2(GSE6136+1)
tmp = names(table(cli$group))
samplesL = split(cli$geo_accession,f = cli$group)
subcounts = sd_selection(df,samplesL)

mislabel = row.names(subcounts[['activate']])
save(mislabel, file = 'mislabel_preselect_GSM142409.rData')

sapply(subcounts,dim)

igraphL = getNetwork(subcounts,fdr=1)
#monoclonal:218 nodes
#oligoclonal:218 nodes
#polyclonal:220 nodes

## step 3.2) Network partition using random walk
## this step constructs modules (geoups of transcriopts) based on their expression levels
cluster = getCluster_methods(igraphL)

## step 4.1-4.3) calculate MCI score for each module
membersL_noweight = getMCI(cluster,subcounts,adjust.size = F)
plotBar_MCI(membersL_noweight,ylim = c(0,4),min = 30)
dev.copy2pdf(file = 'bar_GSE6136_mislabel_GSM142409.pdf')

maxMCIms <- getMaxMCImember(membersL_noweight[[1]],membersL_noweight[[2]],min =30)

maxMCI = getMaxStats(membersL_noweight[['MCI']],maxMCIms[[1]])
## get biomodule with the highest MCI, CTS
CTS = getCTS(maxMCI, maxMCIms[[2]])
#Length: 120

save(CTS, file = 'CTS.GSE6136_mislabel_GSM142409.rData')

par(mfrow = c(2,1))
simMCI = simulationMCI(105,samplesL, df)
plot_MCI_Simulation(maxMCI,simMCI,ylim =c(0,4),las = 2,order = c('resting','activated','lymphoma_marginal','lymphoma_transitional','lymphoma_aggressive'))

Ic = getIc(df,samplesL,CTS)
simIc = simulation_Ic(length(CTS),samplesL,df)
plot_Ic_Simulation(Ic,simIc,ylim = c(0,2.5),order = c('resting','activated','lymphoma_marginal','lymphoma_transitional','lymphoma_aggressive'))

dev.copy2pdf(file = 'MCI_Ic_GSE6136_mislabel_GSM142409.pdf')

############################### related1per: still identifiying 'activated' as tipping point
cli[which(cli$geo_accession == 'GSM142409'),'cell-type:ch1'] = 'lymphoma_marginal'
cli[which(cli$geo_accession == 'GSM142402'),'cell-type:ch1'] = 'activated'

table(cli[,'cell-type:ch1'])
cli$group = cli[,'cell-type:ch1']
samplesL = split(cli$geo_accession,f = cli$group)
subcounts = sd_selection(df,samplesL)
sapply(subcounts,dim)

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
plotBar_MCI(membersL_noweight,ylim = c(0,5),min = 30)

maxMCIms <- getMaxMCImember(membersL_noweight[[1]],membersL_noweight[[2]],min =30)

maxMCI = getMaxStats(membersL_noweight[['MCI']],maxMCIms[[1]])
## get biomodule with the highest MCI, CTS
CTS = getCTS(maxMCI, maxMCIms[[2]])
#Length: 86

save(CTS, file = 'CTS.GSE6136_mislabel_GSM142402.rData')

par(mfrow = c(2,1))
simMCI = simulationMCI(105,samplesL, df)
plot_MCI_Simulation(maxMCI,simMCI,ylim =c(0,4),las = 2,order = c('resting','activated','lymphoma_marginal','lymphoma_transitional','lymphoma_aggressive'))

Ic = getIc(df,samplesL,CTS)
simIc = simulation_Ic(length(CTS),samplesL,df)
plot_Ic_Simulation(Ic,simIc,ylim = c(0,2.5),order = c('resting','activated','lymphoma_marginal','lymphoma_transitional','lymphoma_aggressive'))

dev.copy2pdf(file = 'MCI_Ic_GSE6136_mislabel_GSM142402.pdf')

############################### misclassified1per: identifiying 'lymphoma_agressive' as tipping point
cli[which(cli$geo_accession == 'GSM142402'),'cell-type:ch1'] = 'lymphoma_marginal'
cli[which(cli$geo_accession == 'GSM142418'),'cell-type:ch1'] = 'activated'

table(cli[,'cell-type:ch1'])
cli$group = cli[,'cell-type:ch1']
samplesL = split(cli$geo_accession,f = cli$group)
subcounts = sd_selection(df,samplesL)

mislabel = row.names(subcounts[['activate']])
save(mislabel, file = 'mislabel_preselect_GSM142418.rData')

sapply(subcounts,dim)

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
plotBar_MCI(membersL_noweight,ylim = c(0,5),min = 30)

maxMCIms <- getMaxMCImember(membersL_noweight[[1]],membersL_noweight[[2]],min =30)

maxMCI = getMaxStats(membersL_noweight[['MCI']],maxMCIms[[1]])
## get biomodule with the highest MCI, CTS
CTS = getCTS(maxMCI, maxMCIms[[2]])
#Length: 51

save(CTS, file = 'CTS.GSE6136_mislabel_GSM142418.rData')

par(mfrow = c(2,1))
simMCI = simulationMCI(105,samplesL, df)
plot_MCI_Simulation(maxMCI,simMCI,ylim =c(0,5),las = 2,order = c('resting','activated','lymphoma_marginal','lymphoma_transitional','lymphoma_aggressive'))

Ic = getIc(df,samplesL,CTS)
simIc = simulation_Ic(length(CTS),samplesL,df)
plot_Ic_Simulation(Ic,simIc,ylim = c(0,2.5),order = c('resting','activated','lymphoma_marginal','lymphoma_transitional','lymphoma_aggressive'))

dev.copy2pdf(file = 'MCI_Ic_GSE6136_mislabel_GSM142418.pdf')

