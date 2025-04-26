setwd("~/Desktop/Pitt/Project_sagan_stan/Code_Linan/stan_v3")
library(ggplot2)
library(ggrepel)
library(dplyr)

FONTSIZE = 20
RELSIZE = 0.8

targets = c('PAX5', 'MYC', 'FOXP1', 'SPIB', 'AID', 'TBX21', 'IRF8')

load_data = function(fname, pval_cutoff, diff_cutoff){
  dedf = read.csv(fname, row.names = 1)
  tfs = rownames(dedf)
  n = length(tfs)
  dedf$order = seq(1,n)
  dedf$label = ''
  dedf$color = ''
  dedf$is_target = tfs %in% targets
  for (i in seq(1,n)){
    if (dedf$pvals_adj[i]<pval_cutoff & abs(dedf$diff[i])>diff_cutoff){
      dedf$color[i] = 'diff'
      dedf$label[i] = tfs[i] 
    }
    if (dedf$is_target[i]){
      dedf$color[i] = 'target'
      dedf$label[i] = tfs[i] 
    }
  }
  return(dedf)
}

[1] "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00"
[9] "#CAB2D6" "#6A3D9A" "#FFFF99" "#B15928"

# ========== ========== ========== ========== ========== ==========  
# STAN
# ========== ========== ========== ========== ========== ==========  

dedf = load_data('df2plot_inR/lymphnode_dat_of_gc_by_stan.csv', 
                 pval_cutoff=1e-40, diff_cutoff=0.75)

p = ggplot(dedf) + geom_point(aes(x = diff, y = -log10(pvals_adj), colour = color), size = 1.5) +
  geom_text_repel(aes(x = diff, y = -log10(pvals_adj), label = label), size = 6, max.overlaps = Inf) +
  scale_colour_manual(values = c("#BBBBBB", "#FF7F00", "#6A3D9A")) + theme_bw() +
  xlab("Mean Activity Difference\n<- Other                    Germinal Centers ->") +
  ylab("-log10(pvals_adj)") +
  xlim(-1.6, 1.6) + 
  theme(plot.title = element_text(hjust = 0.5, size = rel(RELSIZE)),
        legend.title = element_blank(),
        legend.position = 'None',
        axis.title = element_text(size = rel(RELSIZE)),
        text = element_text(size=FONTSIZE)) + 
  ggtitle('STAN-inferred TFs Activities by Germinal Centers')
p

pdf('df2plot_inR/lymphnode_volcano_tf.pdf', width=6, height=5.4)
print(p)
dev.off()

# ========== ========== ========== ========== ========== ==========  
# Ridge
# ========== ========== ========== ========== ========== ==========  

dedf_ridge = load_data('df2plot_inR/lymphnode_dat_of_gc_by_ridge.csv', 
                 pval_cutoff=1e-40, diff_cutoff=0.75)

p = ggplot(dedf_ridge) + geom_point(aes(x = diff, y = -log10(pvals_adj), colour = color), size = 1.5) +
  geom_text_repel(aes(x = diff, y = -log10(pvals_adj), label = label), size = 6, max.overlaps = Inf) +
  scale_colour_manual(values = c("#BBBBBB", "#FF7F00", "#6A3D9A")) + theme_bw() +
  xlab("Mean Activity Difference\n<- Other                    Germinal Centers ->") +
  ylab("-log10(pvals_adj)") +
  xlim(-1.6, 1.6) + 
  theme(plot.title = element_text(hjust = 0.5, size = rel(RELSIZE)),
        legend.title = element_blank(),
        legend.position = 'None',
        axis.title = element_text(size = rel(RELSIZE)),
        text = element_text(size=FONTSIZE)) + 
  ggtitle('Ridge-inferred TFs Activities by Germinal Centers') 
p

pdf('df2plot_inR/lymphnode_volcano_tf_ridge.pdf', width=6, height=5.4)
print(p)
dev.off()

# ========== ========== ========== ========== ========== ==========  
# DecoupleR
# ========== ========== ========== ========== ========== ==========  

dedf_dec = load_data('df2plot_inR/lymphnode_dat_of_gc_by_decoupler.csv', 
                       pval_cutoff=1e-40, diff_cutoff=0.75)

p = ggplot(dedf_dec) + geom_point(aes(x = diff, y = -log10(pvals_adj), colour = color), size = 1.5) +
  geom_text_repel(aes(x = diff, y = -log10(pvals_adj), label = label), size = 6, max.overlaps = Inf) +
  scale_colour_manual(values = c("#BBBBBB", "#FF7F00", "#6A3D9A")) + theme_bw() +
  xlab("Mean Activity Difference\n<- Other                    Germinal Centers ->") +
  ylab("-log10(pvals_adj)") +
  xlim(-1.6, 1.6) + 
  theme(plot.title = element_text(hjust = 0.5, size = rel(RELSIZE)),
        legend.title = element_blank(),
        legend.position = 'None',
        axis.title = element_text(size = rel(RELSIZE)),
        text = element_text(size=FONTSIZE)) + 
  ggtitle('DecoupleR-inferred TFs Activities by Germinal Centers') 
p

pdf('df2plot_inR/lymphnode_volcano_tf_dec.pdf', width=6, height=5.4)
print(p)
dev.off()

# ========== ========== ========== ========== ========== ==========  
# gene
# ========== ========== ========== ========== ========== ==========  

dedf_gene = load_data('df2plot_inR/lymphnode_deg_of_gc.csv', 
                     pval_cutoff=1e-40, diff_cutoff=0.75)

q = ggplot(dedf_gene) + geom_point(aes(x = diff, y = -log10(pvals_adj), colour = color), size = 1.5) +
  geom_text_repel(aes(x = diff, y = -log10(pvals_adj), label = label), size = 6, max.overlaps = Inf) +
  scale_colour_manual(values = c("#BBBBBB", "#E31A1C", "#1F78B4")) + theme_bw() +
  xlab("Mean Expression Difference\n<- Other                    Germinal Centers ->") +
  ylab("-log10(pvals_adj)") +
  xlim(-1.6, 1.6) + 
  theme(plot.title = element_text(hjust = 0.5, size = rel(RELSIZE)),
        legend.title = element_blank(),
        legend.position = 'None',
        axis.title = element_text(size = rel(RELSIZE)),
        text = element_text(size=FONTSIZE)) + 
  ggtitle('TF mRNA Expressions by Germinal Centers')

q
pdf('df2plot_inR/lymphnode_volcano_gene.pdf', width=6, height=5.4)
print(q)
dev.off()
