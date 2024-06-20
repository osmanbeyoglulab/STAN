setwd("~/Desktop/Pitt/Project_sagan_stan/Code_Linan/April/notebooks/breast")
library(ggplot2)
library(ggrepel)
library(Seurat)
library(dplyr)
library(RColorBrewer)

FONTSIZE = 25
RELSIZE = 0.8
var = 'diff'

data = read.csv('../df2plot_inR/breast_Stroma.csv', row.names = 1)
data$Score = pmax(pmin(data[[var]], 1.25), -1.25)
data$TF = factor(x = data$TF, 
                      levels = c('ZKSCAN1', 'BCL6', 'KLF10', 'MXI1', 'ELK1', 'JUN', 'VDR', 'NFATC1',
                                 'TEAD4', 'HSF2', 'MEF2C', 'KLF4', 'TP63', 'ELK4', 'PBX3', 'EBF1',
                                 'GTF2I', 'NR2F2', 'SPDEF', 'MAX', 'FLI1', 'KLF5', 'CREB1', 'GABPA',
                                 'NR3C1', 'MAZ', 'NFYB', 'SPI1', 'JUND', 'YY1'))
data$sample = factor(x = data$sample, 
                      levels = c('ER_0', 'TNBC_2', 'ER_1', 'TNBC_1', 'TNBC_0', 'TNBC_3'))

p = ggplot(data, aes(y=sample, x=TF)) + 
  geom_point(aes(size=X.log.p_adj., color=Score)) +
  scale_color_gradientn('TFa (scaled)',
                        colors=brewer.pal(11,'PiYG') %>% rev(),
                        breaks=seq(-1, 1, 1), limits=c(-1.25,1.25)) +
  theme_classic() + ggtitle('Stroma') + 
  labs(size='Significance') +
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
                  labels = c("ns","*","**", '***')) 

pdf('../df2plot_inR/breast_dotplot_Stroma.pdf', width=10, height=4)
print(p)
dev.off()

# ========== ========== ========== ========== ========== ==========  

data = read.csv('../df2plot_inR/breast_diff_Lymphocytes.csv', row.names = 1)
data$Score = pmax(pmin(data[[var]], 2.25), -2.25)
data$TF = factor(x = data$TF, 
                      levels = c('ETV6', 'GATA3', 'KLF5', 'SMAD3', 'ETV5', 'ARNT', 'ASCL2', 'ARID3A',
                                 'GABPA', 'VEZF1', 'FOXP3', 'IRF1', 'IRF4', 'TFAP4', 'BACH2', 'EBF1',
                                 'RFX3', 'STAT2'))
data$sample = factor(x = data$sample, 
                          levels = c('TNBC_0', 'TNBC_1', 'ER_1', 'TNBC_3'))

p = ggplot(data, aes(y=sample, x=TF)) + 
  geom_point(aes(size=X.log.p_adj., color=Score)) +
  scale_color_gradientn('TFa (scaled)',
                        colors=brewer.pal(11,'PiYG') %>% rev(),
                        breaks=seq(-2, 2, 1), limits=c(-2.25,2.25)) +
  theme_classic() + ggtitle('Lymphocytes') + 
  labs(size='Significance') +
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
                  labels = c("ns","*","**", '***')) 

pdf('../df2plot_inR/breast_dotplot_Lymphocytes.pdf', width=8, height=3)
print(p)
dev.off()

# ========== ========== ========== ========== ========== ========== 

data = read.csv('../df2plot_inR/breast_Invasive cancer + stroma + lymphocytes.csv', row.names = 1)
data$Score = pmax(pmin(data[[var]], 1.25), -1.25)
data$TF = factor(x = data$TF, 
                      levels = c('NR3C1', 'FOSL1', 'GTF2I', 'SPIB', 'MEIS1', 'NR2F2', 'ETV4', 'EPAS1',
                                 'ZBTB7A', 'ZKSCAN1', 'TFAP4', 'TP63', 'E2F3', 'GATA3', 'MEF2C',
                                 'NFATC1', 'FLI1', 'STAT5A', 'FOXM1', 'CREB1', 'E2F1', 'MAX'))
data$sample = factor(x = data$sample, 
                          levels = c('ER_0', 'TNBC_1', 'TNBC_0', 'TNBC_2'))

p = ggplot(data, aes(y=sample, x=TF)) + 
  geom_point(aes(size=X.log.p_adj., color=Score)) +
  scale_color_gradientn('TFa (scaled)',
                        colors=brewer.pal(11,'PiYG') %>% rev(),
                        breaks=seq(-1, 1, 1), limits=c(-1.25,1.25)) +
  theme_classic() + ggtitle('Invasive cancer + stroma + lymphocytes') + 
  labs(size='Significance') +
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
                  labels = c("ns","*","**", '***')) 

pdf('../df2plot_inR/breast_dotplot_Invasive.pdf', width=8.4, height=3)
print(p)
dev.off()
