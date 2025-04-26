setwd("~/Desktop/Pitt/Project_sagan_stan/Code_Linan/stan_v3")
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(matrixStats)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)

df_corr_tl = read.csv('df2plot_inR/gbm_corr_tl.csv', row.names=1)
df_corr_tr = read.csv('df2plot_inR/gbm_corr_tr.csv', row.names=1)
df_corr_lr = read.csv('df2plot_inR/gbm_corr_lr.csv', row.names=1)

L = 0.8
p1 = Heatmap(df_corr_tl,
             col = colorRamp2(c(-L, 0, L), c("#4D9221", "white", "#C51B7D")),
             cluster_columns = T,
             cluster_rows = T, 
             column_names_side ="top",
             row_names_side ="left",
             row_title = 'TFs',
             column_title = 'Ligands',
             name = "TF-\nLigand\nPearson r")
p1

p2 = Heatmap(df_corr_tr,
             col = colorRamp2(c(-L, 0, L), c("#4575B4", "white", "#D73027")),
             cluster_columns = T,
             # column_order = colnames(h2_cor_rt)[column_order(p1)],
             cluster_rows = T, 
             column_names_side ="top",
             row_names_side ="left",
             row_title = 'TFs',
             column_title = 'Receptors',
             name = "TF-\nReceptor\nPearson r")
p2

p1+p2

p3 = Heatmap(df_corr_lr %>% t(),
             col = colorRamp2(c(-L, 0, L), c("#8073AC", "white", "#E08214")),
             column_order = colnames(df_corr_tl)[column_order(p1)],
             cluster_rows = T, 
             column_names_side ="bottom",
             row_names_side ="left",
             row_title = 'Receptors',
             name = "Ligand-\nReceptor\nPearson r")
p3

# thre = 0.59
# pdf('df2plot_inR/gbm_heatmap_1.pdf', width=8, height=5.6)
# print(p1+p2)
# dev.off()
# 
# pdf('df2plot_inR/gbm_heatmap_2.pdf', width=5.5, height=3.5)
# print(p3)
# dev.off()

pdf('df2plot_inR/gbm_heatmap_1.pdf', width=10, height=7)
print(p1+p2)
dev.off()

pdf('df2plot_inR/gbm_heatmap_2.pdf', width=6.4, height=4.8)
print(p3)
dev.off()

