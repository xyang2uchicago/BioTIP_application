### This code recalculate the Ic*
### and the the p-value for Delta(Ic*)
### This is based on Zhezhen's code: E:\Zhezhen2019_bk\NB\DNS\code\GSE6136_chen_Ic_CI.R


setwd('F:\\projects\\BioTIP\\result\\GSE6136_Zhezhen\\Holly_desktop\\')
#setwd('F:\\ZheZhen\\NB\\DNS\\result')


GSE6136 = read.table("E:\\Zhezhen2019_bk\\NB\\DNS\\data\\GSE6136_matrix.txt",head = T,
                     sep = '\t',comment = '!', row.names=1)
# row.names(GSE6136) = GSE6136[,1]
# GSE6136 = GSE6136[,-1]
dim(GSE6136) # 22690   26
class(GSE6136)
#[1] "data.frame"

max(GSE6136) # 80843.2

# df = log2(as.matrix(GSE6136))
# summary(df)
df = log2(as.matrix(GSE6136+1))
summary(df)
df[1:3,1:5]
#           GSM142398 GSM142399 GSM142400 GSM142401 GSM142402
# 1415670_at  10.21917  10.29048  9.987548  10.07682  9.827343
# 1415671_at  10.90358  11.15993 10.776186  10.93000 11.468268
# 1415672_at  11.11530  10.89209 11.091303  11.04029 11.109504
dim(df)
#[1] 22690    26

library(stringr)
cli = t(read.delim("E:\\Zhezhen2019_bk\\NB\\DNS\\data\\GSE6136_cli.txt"))
dim(cli)
#[1]  26 38
colnames(cli) = str_split_fixed(cli[1,],'_',2)[,2]
cli = cli[-1,]
cli = data.frame(cli)
cli[,'cell-type:ch1'] = str_split_fixed(cli$characteristics_ch1.1,': ',2)[,2]
cli[,'Ig clonality:ch1'] = str_split_fixed(cli$characteristics_ch1.3,': ',2)[,2]
colnames(cli)[colnames(cli)=='cell-type:ch1'] = 'group'
cli$Row.names = cli[,1]
table(cli$group)
#activated   lymphoma_aggressive     lymphoma_marginal lymphoma_transitional 
# 3                     7                     6                     5 
# resting 
# 5 
cli$group <- factor(cli$group, 
                    levels=c('resting','activated','lymphoma_marginal','lymphoma_transitional', 'lymphoma_aggressive'))

table(cli[,'group'])
# resting             activated     lymphoma_marginal lymphoma_transitional 
#   5                     3                     6                     5 
# lymphoma_aggressive 
#   7 

load("F:/projects/BioTIP/doc/2020_/Applications/outPut.Rdata/GSE6136_robust_analysis/CTS.GSE6136_5states_1per.rData")
length(CTS)# [1] 35

  
sampleL <- split(rownames(cli),f = cli$group)


library(BioTIP)
packageVersion('BioTIP') #[1] '1.5.0'

names(sampleL)
#[1] "resting"               "activated"             "lymphoma_marginal"    
#[4] "lymphoma_transitional" "lymphoma_aggressive"  

Ic  <- getIc(df, sampleL = sampleL, genes = CTS, 
            output = 'Ic',
            fun = "BioTIP", shrink = TRUE, PCC_sample.target = 'none') 
Ic 
#              resting             activated     lymphoma_marginal lymphoma_transitional 
# 0.12516835            1.57747395            0.11492339            0.12532860 
# lymphoma_aggressive 
# 0.08390762 
set.seed(101010)
Ic_simuG <- simulation_Ic(length(CTS), sampleL,  df,  B = 1000,  fun = "BioTIP",   
                          shrink = TRUE, 
                          PCC_sample.target = 'none') 
  


save(Ic, Ic_simuG,
     file='IC_new_1kgenesimulation_BioTIP35G.RData')




pdf(file = 'lines.box_GSE6136_BioTIP35G_simulation.pdf')
par(mfrow = c(2,2))
plot_Ic_Simulation(Ic,  Ic_simuG,  las = 0,  ylim =  c(0,1.5),  
                                order = NULL,  main = 'DNB 35 genes',  
                                ylab = "Ic*",  fun = 'matplot', 
                                which2point = NULL) 
######## calculate delta score
plot_SS_Simulation(Ic,  Ic_simuG, 
                   main = paste("Delta of Ic*", length(CTS),"genes"), 
                   ylab='BioTIP 35 genes') # [1] 0 

plot_Ic_Simulation(Ic,  Ic_simuG,  las = 0,  ylim =  c(0,1.5),  
                   order = NULL,  main = 'BioTIP 35 genes',  
                   ylab = "Ic*",  fun = 'boxplot', 
                   which2point = NULL) 
######## calculate delta score
plot_SS_Simulation(Ic,  Ic_simuG, 
                   main = paste("Delta of Ic*", length(CTS),"genes"), 
                   ylab='BioTIP 35 genes') # [1] 0 

dev.off()

