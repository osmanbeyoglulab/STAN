setwd("~/Desktop/STAN_project/notebooks/lymphnode")
library(ggplot2)
library(ggrepel)
library(Seurat)
library(dplyr)

FONTSIZE = 25
RELSIZE = 0.8

# ========== ========== ========== ========== ========== ==========  
# STAN
# ========== ========== ========== ========== ========== ==========  

dedf = read.csv('../df2plot_inR/lymphnode_dedf_tf.csv', row.names = 1)
tfs = rownames(dedf)
n = length(tfs)
dedf$order = seq(1,n)
dedf$label = ''
dedf$color = ''
for (i in seq(1,n)){
  if (dedf$pvals_adj[i]<1e-50 & abs(dedf$diff[i])>0.015){
    dedf$color[i] = 'diff'
    dedf$label[i] = tfs[i] 
  }
}

p = ggplot(dedf) + geom_point(aes(x = diff, y = -log10(pvals_adj), colour = color), size = 2.5) +
  geom_text_repel(aes(x = diff, y = -log10(pvals_adj), label = label), size = 6, max.overlaps = Inf) +
  scale_colour_manual(values = c("#BBBBBB" , "#882E72")) + theme_bw() +
  xlab("Mean Activity Difference\n<- Other                    Germinal Centers ->") +
  ylab("-log10(pvals_adj)") +
  xlim(-0.06, 0.06) + 
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = rel(RELSIZE)),
        legend.title = element_blank(),
        axis.title = element_text(size = rel(RELSIZE)),
        text = element_text(size=FONTSIZE)) + 
  ggtitle('TFs Activities by Germinal Centers\nInferred by STAN')  + NoLegend() 
p

pdf('../df2plot_inR/lymphnode_volcano_tf.pdf', width=9, height=7.2)
print(p)
dev.off()

# ========== ========== ========== ========== ========== ==========  
# gene
# ========== ========== ========== ========== ========== ==========  

tfs = read.csv('../df2plot_inR/lymphnode_dedf_tf.csv', row.names = 1) %>% rownames() 
dedf_gene = read.csv('../df2plot_inR/lymphnode_dedf_gene.csv', row.names = 1)
dedf_gene = dedf_gene[rownames(dedf_gene) %>% intersect(tfs),]

genes = rownames(dedf_gene)
n = length(genes)
dedf_gene$order = seq(1,n)
dedf_gene$label = ''
dedf_gene$color = ''
for (i in seq(1,n)){ 
  if (dedf_gene$pvals_adj[i]<1e-25 & abs(dedf_gene$diff[i])>1){
    dedf_gene$label[i] = genes[i] 
    dedf_gene$color[i] = 'diff'
  }}

q = ggplot(dedf_gene) + geom_point(aes(x = diff, y = -log10(pvals_adj), colour = color), size = 2.5) +
  geom_text_repel(aes(x = diff, y = -log10(pvals_adj), label = label), size = 6, max.overlaps = Inf) +
  scale_colour_manual(values = c("#BBBBBB" , "#f97b72")) + theme_bw() +
  xlab("Mean Expression Difference\n<- Other                    Germinal Centers ->") +
  ylab("-log10(pvals_adj)") +
  xlim(-3, 3) + 
  # ylim(-0.1, 120) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = rel(RELSIZE)),
        legend.title = element_blank(),
        axis.title = element_text(size = rel(RELSIZE)),
        text = element_text(size=FONTSIZE)) + 
  ggtitle('mRNA Expressions by Germinal Centers')  + NoLegend() 

q
pdf('../df2plot_inR/lymphnode_volcano_gene.pdf', width=9, height=7.2)
print(q)
dev.off()

# ========== ========== ========== ========== ========== ==========  
# decoupler
# ========== ========== ========== ========== ========== ==========  

dedf = read.csv('../df2plot_inR/lymphnode_dedf_tf_decoupler.csv', row.names = 1)
tfs = rownames(read.csv('../df2plot_inR/lymphnode_dedf_tf.csv', row.names = 1))
dedf = dedf[rownames(dedf) %>% intersect(tfs),]

tfs = rownames(dedf)
n = length(tfs)
dedf$order = seq(1,n)
dedf$label = ''
dedf$color = ''
for (i in seq(1,80)){ 
  if (dedf$pvals_adj[i]<1e-50 & abs(dedf$diff[i])>1){
    dedf$label[i] = tfs[i] 
    dedf$color[i] = 'diff'
  }}
for (i in seq(n-19,n)){
  if (dedf$pvals_adj[i]<1e-25 & abs(dedf$diff[i])>0.5){
    dedf$label[i] = tfs[i] 
    dedf$color[i] = 'diff'
  }}

p = ggplot(dedf) + geom_point(aes(x = diff, y = -log10(pvals_adj), colour = color), size = 2.5) +
  geom_text_repel(aes(x = diff, y = -log10(pvals_adj), label = label), size = 6, max.overlaps = Inf) +
  scale_colour_manual(values = c("#BBBBBB" , "#882E72")) + theme_bw() +
  xlab("Mean Activity Difference\n<- Other                    Germinal Centers ->") +
  ylab("-log10(pvals_adj)") +
  xlim(-3, 4.5) + ylim(-0.1, 200) +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = rel(RELSIZE)),
        legend.title = element_blank(),
        axis.title = element_text(size = rel(RELSIZE)),
        text = element_text(size=FONTSIZE)) + 
  ggtitle('TFs Activities by Germinal Centers\nInferred by decoupleR')  + NoLegend() 
p

pdf('../df2plot_inR/lymphnode_volcano_tf_decoupler.pdf', width=9, height=7.2)
print(p)
dev.off()

