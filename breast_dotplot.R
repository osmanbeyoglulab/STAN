setwd("~/Desktop/Pitt/Project_sagan_stan/Code_Linan/stan_v3")
library(ggplot2)
library(ggrepel)
library(Seurat)
library(dplyr)
library(RColorBrewer)

FONTSIZE = 18
RELSIZE = 0.8
var = 'diff'
id_list = c('CID4290'='ER_1', 'CID4535'='ER_2', 'CID44971'='TNBC_4',
  'CID4465'='TNBC_3', '1142243F'='TNBC_1', '1160920F'='TNBC_2')

data = read.csv('df2plot_inR/breast_Stroma.csv', row.names = 1)
data$Score = pmax(pmin(data[[var]], 1.75), -1.75)
data$TF = factor(x = data$TF, 
                      levels = c('ZNF76', 'IKZF1', 'CREB3', 'NR2C2', 'GATA6', 'MAFB', 'HINFP', 'CREB1',
                                 'GABPA', 'SOX17', 'NR5A2', 'SPDEF', 'STAT3', 'TCF4', 'ARNTL', 'HSF2',
                                 'ZKSCAN1', 'HOXC6', 'GMEB1', 'MEIS1', 'VDR', 'ETV6', 'ELK1', 'TEAD4',
                                 'MEF2C', 'DBP', 'NFATC1'))
data$sample = factor(x = id_list[ data$sample] %>% unname(), 
                      levels = id_list[c('CID4465', 'CID4290', '1142243F', 'CID44971', '1160920F', 'CID4535')])

p = ggplot(data, aes(y=sample, x=TF)) + 
  geom_point(aes(size=X.log.p_adj., color=Score)) +
  scale_color_gradientn('TFa (scaled)',
                        colors=brewer.pal(11,'PiYG') %>% rev(),
                        breaks=seq(-1.5, 1.5, 1.5), limits=c(-1.75,1.75)) +
  theme_classic() + ggtitle('Stroma') + 
  labs(size='-log(p_adj)') +
  theme(panel.grid = element_blank(),
        strip.placement = "outside",
        strip.background = element_rect(colour="grey", fill="white"),
        strip.text.y.left = element_text(face='bold', size = rel(RELSIZE), angle = 0),
        plot.title = element_text(hjust = 0.5, size = rel(RELSIZE)),
        legend.title = element_text(size = rel(RELSIZE^3)),
        legend.text = element_text(size = rel(RELSIZE^3)),
        legend.key.height = unit(0.6,"cm"),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        text = element_text(size=FONTSIZE)) +
  scale_size_area(breaks=c(1, 2, 3, 4),
                  labels = c("1","2","3", '4+')) 
p
pdf('df2plot_inR/breast_dotplot_Stroma.pdf', width=9, height=3.6)
print(p)
dev.off()

# ========== ========== ========== ========== ========== ==========  

data = read.csv('df2plot_inR/breast_Lymphocytes.csv', row.names = 1)
data$Score = pmax(pmin(data[[var]], 1.75), -1.75)
data$TF = factor(x = data$TF, 
                      levels = c('RFX5', 'STAT2', 'FOXP3', 'GATA6', 'BACH2', 'STAT5B', 'HOXC6', 'NR2C2',
                                 'ZKSCAN1', 'NFAT5', 'POU2F1', 'IKZF1', 'MAFB', 'HES1', 'STAT3', 'CREB1',
                                 'HINFP', 'SPDEF', 'ZNF76', 'GABPA', 'NR5A2', 'SOX17', 'TCF4'))
data$sample = factor(x = id_list[ data$sample] %>% unname(),  
                          levels = id_list[c('1160920F', '1142243F', 'CID44971', 'CID4535')])

p = ggplot(data, aes(y=sample, x=TF)) + 
  geom_point(aes(size=X.log.p_adj., color=Score)) +
  scale_color_gradientn('TFa (scaled)',
                        colors=brewer.pal(11,'PiYG') %>% rev(),
                        breaks=seq(-1.5, 1.5, 1.5), limits=c(-1.75,1.75)) +
  theme_classic() + ggtitle('Lymphocytes') + 
  labs(size='-log(p_adj)') +
  theme(panel.grid = element_blank(),
        strip.placement = "outside",
        strip.background = element_rect(colour="grey", fill="white"),
        strip.text.y.left = element_text(face='bold', size = rel(RELSIZE), angle = 0),
        plot.title = element_text(hjust = 0.5, size = rel(RELSIZE)),
        legend.title = element_text(size = rel(RELSIZE^3)),
        legend.text = element_text(size = rel(RELSIZE^3)),
        legend.key.height = unit(0.4,"cm"),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        text = element_text(size=FONTSIZE)) +
  scale_size_area(breaks=c(1, 2, 3, 4),
                  labels = c("1","2","3", '4+')) 

pdf('df2plot_inR/breast_dotplot_Lymphocytes.pdf', width=9, height=3)
print(p)
dev.off()

# ========== ========== ========== ========== ========== ========== 

data = read.csv('df2plot_inR/breast_Invasive cancer + stroma + lymphocytes.csv', row.names = 1)
data$Score = pmax(pmin(data[[var]], 1.75), -1.75)
data$TF = factor(x = data$TF, 
                      levels = c('RFX3', 'ZNF76', 'MAFB', 'NR3C1', 'GMEB1', 'GATA6', 'HOXC6', 'ZKSCAN1',
                                 'IKZF1', 'CREB1', 'STAT3', 'TCF4', 'SOX17', 'TAL1', 'NR2C2', 'NR5A2',
                                 'GABPA', 'SPDEF', 'HINFP', 'E2F1', 'FOXM1'))
data$sample = factor(x = id_list[ data$sample] %>% unname(),  
                          levels = id_list[c('CID4290', '1160920F', '1142243F', 'CID4465')])

p = ggplot(data, aes(y=sample, x=TF)) + 
  geom_point(aes(size=X.log.p_adj., color=Score)) +
  scale_color_gradientn('TFa (scaled)',
                        colors=brewer.pal(11,'PiYG') %>% rev(),
                        breaks=seq(-1.5, 1.5, 1), limits=c(-1.75,1.75)) +
  theme_classic() + ggtitle('Invasive cancer + stroma + lymphocytes') + 
  labs(size='-log(p_adj)') +
  theme(panel.grid = element_blank(),
        strip.placement = "outside",
        strip.background = element_rect(colour="grey", fill="white"),
        strip.text.y.left = element_text(face='bold', size = rel(RELSIZE), angle = 0),
        plot.title = element_text(hjust = 0.5, size = rel(RELSIZE)),
        legend.title = element_text(size = rel(RELSIZE^3)),
        legend.text = element_text(size = rel(RELSIZE^3)),
        legend.key.height = unit(0.4,"cm"),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        text = element_text(size=FONTSIZE)) +
  scale_size_area(breaks=c(1, 2, 3, 4),
                  labels = c("1","2","3", '4+')) 

pdf('df2plot_inR/breast_dotplot_Invasive.pdf', width=9, height=3)
print(p)
dev.off()

# ========== ========== ========== ========== ========== ========== 

data = read.csv('df2plot_inR/breast_tf_by_pathology_CID44971.csv', row.names = 1)
data$Score = pmax(pmin(data[[var]], 1.75), -1.75)
data$TF = factor(x = data$TF, 
                 levels = c('STAT3', 'ARNTL', 'TP73', 'IKZF1', 'NR2C2', 'SREBF2', 'FOXP3', 'STAT5B',
                            'MITF', 'BACH1', 'HSF2', 'KLF5', 'ZNF250', 'MBD2', 'FOXA1', 'MAF',
                            'KMT2A', 'ZFX', 'ETV6', 'LEF1', 'RELA', 'MAZ', 'NR5A2', 'STAT2', 'IRF2',
                            'RFX3', 'BACH2', 'POU2F1', 'VDR', 'MNT', 'MAX', 'SREBF1', 'HBP1',
                            'PAX6'))
data$Pathology = factor(x = data$Pathology %>% unname(),  
                     levels = c('Normal + stroma + lymphocytes', 'Stroma + adipose tissue',
                                'Invasive cancer + lymphocytes', 'Lymphocytes', 'DCIS', 'Stroma'))

p = ggplot(data, aes(y=Pathology, x=TF)) + 
  geom_point(aes(size=X.log.p_adj., color=Score)) +
  scale_color_gradientn('TFa (scaled)',
                        colors=brewer.pal(11,'PiYG') %>% rev(),
                        breaks=seq(-1.5, 1.5, 1.5), limits=c(-1.75,1.75)) +
  theme_classic() + ggtitle('TNBC_4') + 
  labs(size='-log(p_adj)') +
  theme(panel.grid = element_blank(),
        strip.placement = "outside",
        strip.background = element_rect(colour="grey", fill="white"),
        strip.text.y.left = element_text(face='bold', size = rel(RELSIZE), angle = 0),
        plot.title = element_text(hjust = 0.5, size = rel(RELSIZE)),
        legend.title = element_text(size = rel(RELSIZE^3)),
        legend.text = element_text(size = rel(RELSIZE^3)),
        legend.key.height = unit(0.4,"cm"),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        text = element_text(size=FONTSIZE)) +
  scale_size_area(breaks=c(1, 2, 3, 4),
                  labels = c("1","2","3", '4+')) 

pdf('df2plot_inR/breast_dotplot_TNBC_4.pdf', width=12, height=3)
print(p)
dev.off()
