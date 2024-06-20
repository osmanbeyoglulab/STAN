setwd("~/Desktop/STAN_project/notebooks/breast")
library(ggplot2)
library(ggrepel)
library(Seurat)
library(dplyr)
library(RColorBrewer)

FONTSIZE = 25
RELSIZE = 0.8

sample_list = c('ER_0', 'ER_1', 'TNBC_0', 'TNBC_1', 'TNBC_2',  'TNBC_3')

data_list = list()
for (sample in sample_list){
  data_list[[sample]] = read.csv(sprintf('../df2plot_inR/breast_tf_by_pathology_%s.csv', sample), row.names = 1)
}
var = 'diff'

threshold_list = list('ER_0' = c(1, 1.1), 
                      'ER_1' = c(1.5, 1.65), 
                      'TNBC_0' = c(2, 2.2), 
                      'TNBC_1' = c(2, 2.2), 
                      'TNBC_2' = c(3, 3.2), 
                      'TNBC_3' = c(1.5, 1.65))

for (sample in sample_list){
  data_list[[sample]]$Score = pmax(pmin(data_list[[sample]][[var]], 
                                        threshold_list[[sample]][2]),
                                   -threshold_list[[sample]][2])
}

data_list$ER_0$TF = factor(x = data_list$ER_0$TF,
                           levels = c('ESRRA', 'NFE2', 'MXI1', 'BCL6', 'KLF10', 'HOXA4', 'ELK3', 'STAT1',
                                      'FOSL1', 'NR3C1', 'SPIB', 'RFX3', 'RUNX2', 'MAZ', 'NFYB', 'GATA3',
                                      'MEF2C', 'E2F3', 'HSF2', 'ZNF250', 'SPDEF', 'BATF', 'NR5A2', 'FOXM1',
                                      'RARG'))
data_list$ER_0$Pathology = factor(x = data_list$ER_0$Pathology,
                                  levels = c('Stroma', 'Invasive cancer + stroma',
                                             'Invasive cancer + stroma + lymphocytes'))


data_list$ER_1$TF = factor(x = data_list$ER_1$TF,
                           levels = c('LEF1', 'SP2', 'AR', 'E2F2', 'FOXP3', 'E2F7', 'ELK1', 'BACH2', 'EBF1',
                                      'SOX2', 'NR3C1', 'TCF3', 'TEAD4', 'CEBPA', 'SNAI2', 'SPI1', 'GATA3',
                                      'VEZF1', 'AHR', 'KLF5', 'HOXA9', 'NFE2L2', 'ZNF250', 'GTF2B', 'TCF7L2',
                                      'ETV6', 'VDR', 'SREBF2', 'CREB1', 'GABPA', 'GMEB2', 'RXRA', 'FLI1',
                                      'SETDB1', 'ZNF143'))
data_list$ER_1$Pathology = factor(x = data_list$ER_1$Pathology,
                                  levels = c('Invasive cancer + adipose tissue + lymphocytes', 'Adipose tissue',
                                             'Invasive cancer', 'Invasive cancer + lymphocytes', 'Lymphocytes',
                                             'Stroma'))

data_list$TNBC_0$TF = factor(x = data_list$TNBC_0$TF,
                             levels = c('STAT2', 'GMEB1', 'SPIB', 'CEBPA', 'RFX5', 'FOXP3', 'RFX3', 'TFAP2C',
                                        'CLOCK', 'ZBTB7B', 'GABPA', 'ARNT', 'ASCL2', 'E2F1', 'FLI1', 'MAX',
                                        'EPAS1', 'NFE2L2', 'ZBTB7A', 'FOXM1', 'ZKSCAN1', 'ELK1', 'HSF2', 'ETV4',
                                        'KLF4'))
data_list$TNBC_0$Pathology = factor(x = data_list$TNBC_0$Pathology,
                                    levels = c('Lymphocytes', 'TLS', 'Necrosis',
                                               'Invasive cancer + stroma + lymphocytes', 'Stroma'))

data_list$TNBC_1$TF = factor(x = data_list$TNBC_1$TF,
                             levels = c('TEAD4', 'ELF1', 'FOXP1', 'FOXA1', 'GATA3', 'ELK3', 'PBX3', 'PGR',
                                        'ELK4', 'MEIS1', 'TP63', 'IRF1', 'IRF4', 'NFATC1', 'FOXP3', 'TFAP4',
                                        'HINFP', 'HBP1', 'STAT2', 'GLI2', 'HNF4G', 'THAP1', 'FLI1', 'MAX',
                                        'CREB1', 'ETS2', 'MZF1', 'ETV5', 'KLF5', 'SMAD3', 'ETV4', 'EOMES',
                                        'FOXD2'))
data_list$TNBC_1$Pathology = factor(x = data_list$TNBC_1$Pathology,
                                    levels = c('Cancer trapped in lymphocyte aggregation', 'DCIS', 'Lymphocytes',
                                               'Invasive cancer + stroma + lymphocytes', 'Normal glands + lymphocytes',
                                               'Adipose tissue', 'Stroma'))

data_list$TNBC_2$TF = factor(x = data_list$TNBC_2$TF,
                             levels = c('LHX2', 'PGR', 'MEF2C', 'NFATC1', 'STAT5A', 'EGR2', 'ELK1', 'ATF1',
                                        'ELF2', 'ELF3', 'NR2F2', 'JUN', 'NFIC', 'GTF2I', 'NR3C1', 'ZNF263',
                                        'HEY1', 'BCL11A', 'RFX5', 'ESR1', 'MAFB', 'ZNF236'))
data_list$TNBC_2$Pathology = factor(x = data_list$TNBC_2$Pathology,
                                    levels = c('Normal duct', 'Invasive cancer + stroma + lymphocytes', 'Stroma'))

data_list$TNBC_3$TF = factor(x = data_list$TNBC_3$TF,
                             levels = c('STAT3', 'TP73', 'BACH1', 'HSF2', 'MEF2C', 'ARID3A', 'ETV6', 'FOXA1',
                                        'LHX2', 'E2F8', 'KLF5', 'GABPA', 'ZFX', 'EHF', 'JUND', 'BACH2', 'IRF1',
                                        'STAT2', 'IRF2', 'RFX3', 'FOXP3', 'EBF1', 'YY1', 'PAX6', 'VDR', 'MAZ',
                                        'RELA', 'TFDP1', 'E2F1', 'VEZF1', 'MAX', 'MNT'))
data_list$TNBC_3$Pathology = factor(x = data_list$TNBC_3$Pathology,
                                    levels = c('Normal + stroma + lymphocytes', 'Stroma', 'Stroma + adipose tissue',
                                               'DCIS', 'Invasive cancer + lymphocytes', 'Lymphocytes'))

dotplot = function(sample, legend.key.height = 0.6, legend.break = 1){
  threshold = threshold_list[[sample]]
  data_copy = data_list[[sample]]
  p = ggplot(data_copy, aes(y=Pathology, x=TF)) + 
    geom_point(aes(size=X.log.p_adj., color=Score)) +
    scale_color_gradientn('TFa (scaled)',
                          colors=brewer.pal(11,'PiYG') %>% rev(),
                          breaks=seq(-threshold[1],threshold[1],legend.break), 
                          limits=c(-threshold[2],threshold[2])) +
    theme_classic() + ggtitle(sample) + 
    labs(size='Significance') +
    theme(panel.grid = element_blank(),
          strip.placement = "outside",
          strip.background = element_rect(colour="grey", fill="white"),
          strip.text.y.left = element_text(face='bold', size = rel(RELSIZE), angle = 0),
          plot.title = element_text(hjust = 0.5, size = rel(RELSIZE)),
          legend.title = element_text(size = rel(RELSIZE^3)),
          legend.text = element_text(size = rel(RELSIZE^3)),
          legend.key.height = unit(legend.key.height,"cm"),
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          text = element_text(size=FONTSIZE)) +
    scale_size_area(breaks=c(1, 2, 3, 4),
                    labels = c("ns","*","**", '***')) 
  return(p)
} 

sample = 'ER_0'
p = dotplot(sample, legend.key.height = 0.4)
pdf(sprintf('../df2plot_inR/breast_dotplot_%s.pdf', sample), width=12.8, height=3.2)
print(p)
dev.off()

sample = 'ER_1'
p = dotplot(sample, legend.break = 0.5)
pdf(sprintf('../df2plot_inR/breast_dotplot_%s.pdf', sample), width=16.8, height=4)
print(p)
dev.off()

sample = 'TNBC_0'
p = dotplot(sample)
pdf(sprintf('../df2plot_inR/breast_dotplot_%s.pdf', sample), width=14, height=4)
print(p)
dev.off()

sample = 'TNBC_1'
p = dotplot(sample)
pdf(sprintf('../df2plot_inR/breast_dotplot_%s.pdf', sample), width=16, height=4)
print(p)
dev.off()

sample = 'TNBC_2'
p = dotplot(sample, legend.key.height = 0.45, legend.break = 1.5)
pdf(sprintf('../df2plot_inR/breast_dotplot_%s.pdf', sample), width=12.8, height=3.2)
print(p)
dev.off()

sample = 'TNBC_3'
p = dotplot(sample, legend.break = 0.5)
pdf(sprintf('../df2plot_inR/breast_dotplot_%s.pdf', sample), width=14, height=4)
print(p)
dev.off()

