setwd("~/Desktop/STAN_project/notebooks/lymphnode")
library(ggplot2)
library(ggrepel)
library(Seurat)
library(dplyr)
library(RColorBrewer)

FONTSIZE = 25
RELSIZE = 0.8

# ========== ========== ========== ========== ========== ==========  
# STAN
# ========== ========== ========== ========== ========== ==========  

data = read.csv('../df2plot_inR/lymphnode_tf_by_celltype.csv', row.names = 1)
threshold = 10
data$p_adj = -log10(data$p_adj+1e-10)
data$Score = pmax(pmin(data$coef, threshold+2), -threshold-2)
data$ct = factor(x = data$ct, 
                      levels = c('NKT', 'Monocytes', 'VSMC', 'B_Cycling', 'B_IFN', 'T_TfR', 'T_Treg',
                                 'B_plasma', 'Macrophages', 'B_mem', 'DC', 'FDC', 'B_GC', 'B_preGC',
                                 'T_TIM3+', 'B_activated', 'NK', 'T_CD8+', 'ILC', 'T_CD4+'))
data$tf = factor(x = data$tf, 
                      levels = c('PBX1', 'PPARG', 'TCF21', 'STAT1', 'FOXP2', 'PGR', 'IRF1', 'IRF2',
                                 'STAT2', 'E2F7', 'ETS1', 'GTF2B', 'NR5A2', 'BCL11A', 'RFX3', 'SPIB',
                                 'EBF1', 'CREB1', 'PAX5', 'BACH2', 'POU2F2', 'STAT3', 'STAT4', 'KLF1',
                                 'CDX2', 'FOXP3', 'TBX21', 'MYB', 'STAT5B', 'ETV6', 'MAX', 'KMT2A',
                                 'E2F1', 'FOXM1'))

p = ggplot(data, aes(y=tf, x=ct)) + 
  geom_point(aes(size=p_adj, color=Score)) +
  scale_color_gradientn('STAN\nTF\nActivity\nScore',
                        colors=brewer.pal(11,'PiYG') %>% rev(),
                        breaks=seq(-threshold,threshold,5), limits=c(-threshold-2,threshold+2)) +
  theme_classic() +
  labs(size='Significance') +
  theme(panel.grid = element_blank(),
        strip.placement = "outside",
        strip.background = element_rect(colour="grey", fill="white"),
        strip.text.y.left = element_text(face='bold', size = rel(RELSIZE), angle = 0),
        plot.title = element_text(hjust = 0.5, size = rel(RELSIZE)),
        legend.title = element_text(size = rel(RELSIZE^2)),
        legend.text = element_text(size = rel(RELSIZE^2)),
        legend.key.height = unit(0.8,"cm"),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        text = element_text(size=FONTSIZE)) +
  scale_size_area(breaks=c(1, 2, 3, 4),
                  labels = c("ns","*","**", '***')) 

p
pdf('../df2plot_inR/lymphnode_dotplot_tf.pdf', width=9, height=12)
print(p)
dev.off()

# ========== ========== ========== ========== ========== ==========  
# gene
# ========== ========== ========== ========== ========== ==========  

data = read.csv('../df2plot_inR/lymphnode_gene_by_celltype.csv', row.names = 1)
threshold = 5
data$p_adj = -log10(data$p_adj+1e-10)
data$Score = pmax(pmin(data$coef, threshold+1), -threshold-1)
data$ct = factor(x = data$ct, 
                      levels = c('NKT', 'Monocytes', 'VSMC', 'B_Cycling', 'B_IFN', 'T_TfR', 'T_Treg',
                                 'B_plasma', 'Macrophages', 'B_mem', 'DC', 'FDC', 'B_GC', 'B_preGC',
                                 'T_TIM3+', 'B_activated', 'NK', 'T_CD8+', 'ILC', 'T_CD4+'))
data$tf = factor(x = data$tf, 
                      levels = c('PBX1', 'PPARG', 'TCF21', 'STAT1', 'FOXP2', 'PGR', 'IRF1', 'IRF2',
                                 'STAT2', 'E2F7', 'ETS1', 'GTF2B', 'NR5A2', 'BCL11A', 'RFX3', 'SPIB',
                                 'EBF1', 'CREB1', 'PAX5', 'BACH2', 'POU2F2', 'STAT3', 'STAT4', 'KLF1',
                                 'CDX2', 'FOXP3', 'TBX21', 'MYB', 'STAT5B', 'ETV6', 'MAX', 'KMT2A',
                                 'E2F1', 'FOXM1'))

p = ggplot(data, aes(y=tf, x=ct)) + 
  geom_point(aes(size=p_adj, color=Score)) +
  scale_color_gradientn('mRNA\nExpression\nScore',
                        colors=brewer.pal(11,'RdBu') %>% rev(),
                        breaks=seq(-threshold,threshold,2.5), limits=c(-threshold-1,threshold+1)) +
  theme_classic() +
  labs(size='Significance') +
  theme(panel.grid = element_blank(),
        strip.placement = "outside",
        strip.background = element_rect(colour="grey", fill="white"),
        strip.text.y.left = element_text(face='bold', size = rel(RELSIZE), angle = 0),
        plot.title = element_text(hjust = 0.5, size = rel(RELSIZE)),
        legend.title = element_text(size = rel(RELSIZE^2)),
        legend.text = element_text(size = rel(RELSIZE^2)),
        legend.key.height = unit(0.8,"cm"),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        text = element_text(size=FONTSIZE)) +
  scale_size_area(breaks=c(1, 2, 3, 4),
                  labels = c("ns","*","**", '***')) 

p
pdf('../df2plot_inR/lymphnode_dotplot_gene.pdf', width=9, height=12)
print(p)
dev.off()

# ========== ========== ========== ========== ========== ==========  
# decoupler
# ========== ========== ========== ========== ========== ==========   

data = read.csv('../df2plot_inR/lymphnode_tf_by_celltype_decoupler.csv', row.names = 1)
threshold = 10
data$p_adj = -log10(data$p_adj+1e-10)
data$Score = pmax(pmin(data$coef, threshold+2), -threshold-2)
data$ct = factor(x = data$ct, 
                      levels = c('NKT', 'Monocytes', 'VSMC', 'B_Cycling', 'B_IFN', 'T_TfR', 'T_Treg',
                                 'B_plasma', 'Macrophages', 'B_mem', 'DC', 'FDC', 'B_GC', 'B_preGC',
                                 'T_TIM3+', 'B_activated', 'NK', 'T_CD8+', 'ILC', 'T_CD4+'))
data$tf = factor(x = data$tf, 
                      levels = c('PBX1', 'PPARG', 'TCF21', 'STAT1', 'FOXP2', 'PGR', 'IRF1', 'IRF2',
                                 'STAT2', 'E2F7', 'ETS1', 'GTF2B', 'NR5A2', 'BCL11A', 'RFX3', 'SPIB',
                                 'EBF1', 'CREB1', 'PAX5', 'BACH2', 'POU2F2', 'STAT3', 'STAT4', 'KLF1',
                                 'CDX2', 'FOXP3', 'TBX21', 'MYB', 'STAT5B', 'ETV6', 'MAX', 'KMT2A',
                                 'E2F1', 'FOXM1'))

p = ggplot(data, aes(y=tf, x=ct)) + 
  geom_point(aes(size=p_adj, color=Score)) +
  scale_color_gradientn('decoupleR\nTF\nActivity\nScore',
                        colors=brewer.pal(11,'PiYG') %>% rev(),
                        breaks=seq(-threshold,threshold,5), limits=c(-threshold-2,threshold+2)) +
  theme_classic() +
  labs(size='Significance') +
  theme(panel.grid = element_blank(),
        strip.placement = "outside",
        strip.background = element_rect(colour="grey", fill="white"),
        strip.text.y.left = element_text(face='bold', size = rel(RELSIZE), angle = 0),
        plot.title = element_text(hjust = 0.5, size = rel(RELSIZE)),
        legend.title = element_text(size = rel(RELSIZE^2)),
        legend.text = element_text(size = rel(RELSIZE^2)),
        legend.key.height = unit(0.8,"cm"),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        text = element_text(size=FONTSIZE)) +
  scale_size_area(breaks=c(1, 2, 3, 4),
                  labels = c("ns","*","**", '***')) 

p
pdf('../df2plot_inR/lymphnode_dotplot_tf_decoupler.pdf', width=9, height=12)
print(p)
dev.off()
