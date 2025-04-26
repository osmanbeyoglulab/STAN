setwd("~/Desktop/Pitt/Project_sagan_stan/Code_Linan/stan_v3")
library(ggplot2)
library(ggrepel)
library(Seurat)
library(dplyr)
library(RColorBrewer)

FONTSIZE = 18
RELSIZE = 0.8

cts_stan = c('B_Cycling', 'B_GC', 'FDC', 'B_IFN', 'T_CD8+', 'T_Treg', 'T_CD4+',
             'T_TfR', 'B_naive', 'NK', 'B_activated', 'B_plasma', 'Macrophages',
             'Endo', 'Monocytes', 'NKT', 'VSMC', 'B_mem', 'B_preGC', 'DC')
tfs_stan = c('STAT3', 'FOXP2', 'STAT1', 'MAZ', 'LHX2', 'ELK3', 'TCF21', 'SOX2',
             'HNF1B', 'ETV6', 'NFYC', 'RFX3', 'ETV4', 'MAFB', 'ATF5', 'PBX1',
             'BACH2', 'CREB1', 'HSF1', 'FOXM1', 'KMT2A', 'E2F1', 'MAX', 'GTF2B',
             'E2F7', 'LEF1', 'KLF1', 'STAT4', 'STAT6', 'CREB3', 'FOXP3', 'MYB',
             'ZFHX3', 'STAT5B', 'TBX21', 'JUND', 'IRF2', 'STAT2')

# ========== ========== ========== ========== ========== ==========  
# STAN
# ========== ========== ========== ========== ========== ==========  

data = read.csv('df2plot_inR/lymphnode_tf_by_celltype.csv', row.names = 1)
cts = cts_stan
tfs = tfs_stan

threshold = 30
# data$p_adj = -log10(data$p_adj+1e-10)
data$Score = pmax(pmin(data$coef, threshold+2), -threshold-2)
data$ct = factor(x = data$ct, levels = cts)
data$tf = factor(x = data$tf, levels = tfs)

p = ggplot(data, aes(y=tf, x=ct)) +
  geom_point(aes(size=neg_log_p_adj, color=Score)) +
  scale_color_gradientn('Cell type\nspecific\nTF Score',
                        colors=brewer.pal(11,'PiYG') %>% rev(),
                        breaks=seq(-threshold,threshold,20), limits=c(-threshold-2,threshold+2)) +
  theme_classic() +
  labs(size='-log(p_adj)') +
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
  scale_size_area(breaks=c(10, 20, 30, 40),
                  labels = c("10","20","30", '40+')) 

p
pdf('df2plot_inR/lymphnode_dotplot_tf.pdf', width=7, height=10)
print(p)
dev.off()

data = read.csv('df2plot_inR/lymphnode_tf_by_celltype_ridge.csv', row.names = 1)
cts = cts_stan
tfs = tfs_stan

# data$p_adj = -log10(data$p_adj+1e-10)
data$Score = pmax(pmin(data$coef, threshold+2), -threshold-2)
data$ct = factor(x = data$ct, levels = cts)
data$tf = factor(x = data$tf, levels = tfs)

p = ggplot(data, aes(y=tf, x=ct)) +
  geom_point(aes(size=neg_log_p_adj, color=Score)) +
  scale_color_gradientn('Cell type\nspecific\nTF Score\n(Ridge)',
                        colors=brewer.pal(11,'PiYG') %>% rev(),
                        breaks=seq(-threshold,threshold,20), limits=c(-threshold-2,threshold+2)) +
  theme_classic() +
  labs(size='-log(p_adj)') +
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
  scale_size_area(breaks=c(10, 20, 30, 40),
                  labels = c("10","20","30", '40+')) 

p
pdf('df2plot_inR/lymphnode_dotplot_tf_ridge.pdf', width=7, height=10)
print(p)
dev.off()

data = read.csv('df2plot_inR/lymphnode_tf_by_celltype_decoupler.csv', row.names = 1)
cts = cts_stan
tfs = tfs_stan

# data$p_adj = -log10(data$p_adj+1e-10)
data$Score = pmax(pmin(data$coef, threshold+2), -threshold-2)
data$ct = factor(x = data$ct, levels = cts)
data$tf = factor(x = data$tf, levels = tfs)

p = ggplot(data, aes(y=tf, x=ct)) +
  geom_point(aes(size=neg_log_p_adj, color=Score)) +
  scale_color_gradientn('Cell type\nspecific\nTF Score\n(DecoupleR)',
                        colors=brewer.pal(11,'PiYG') %>% rev(),
                        breaks=seq(-threshold,threshold,20), limits=c(-threshold-2,threshold+2)) +
  theme_classic() +
  labs(size='-log(p_adj)') +
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
  scale_size_area(breaks=c(10, 20, 30, 40),
                  labels = c("10","20","30", '40+')) 

p
pdf('df2plot_inR/lymphnode_dotplot_tf_decoupler.pdf', width=7, height=10)
print(p)
dev.off()

# ========== ========== ========== ========== ========== ==========  
# gene
# ========== ========== ========== ========== ========== ==========  

data = read.csv('df2plot_inR/lymphnode_mrna_by_celltype.csv', row.names = 1)
cts = cts_stan %>% intersect(data$ct)
tfs = tfs_stan %>% intersect(data$tf)
threshold = 20
# data$p_adj = -log10(data$p_adj+1e-10)
data$Score = pmax(pmin(data$coef, threshold+2), -threshold-2)
data$ct = factor(x = data$ct, levels = cts)
data$tf = factor(x = data$tf, levels = tfs)

p = ggplot(data, aes(y=tf, x=ct)) + 
  geom_point(aes(size=neg_log_p_adj, color=Score)) +
  scale_color_gradientn('Cell type\nspecific\nmRNA\nScore',
                        colors=brewer.pal(11,'RdBu') %>% rev(),
                        breaks=seq(-threshold,threshold,10), limits=c(-threshold-2,threshold+2)) +
  theme_classic() +
  labs(size='-log(p_adj)') +
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
                  labels = c("1","2","3", '4+')) 

p
pdf('df2plot_inR/lymphnode_dotplot_gene.pdf', width=7, height=10)
print(p)
dev.off()

