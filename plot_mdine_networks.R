# Plotting networks from mdine fit using igraph pacakge

library(igraph)
library(phyloseq)
library(mdine)

source("sim_eval_functions.R")
source('plot_networks.R')

load("cleaned_data_ileum.RData")

#Family analysis -- choosing top (named) families
n.fam = 15
named.fam = counts_ileum$family[,colnames(counts_ileum$family)!="f__"]
top.fam = names(head(sort(colMeans(named.fam),decreasing = TRUE), n.fam))
not.top.fam = colnames(counts_ileum$family)[!colnames(counts_ileum$family) %in% top.fam]
counts.top.fam = cbind(counts_ileum$family[,top.fam], rowSums(counts_ileum$family[,not.top.fam]))
colnames(counts.top.fam) = c(top.fam, "ref")

#Recoding diagnosis/antibiotics to drop zero categories
dat.ileum$diagnosis = factor(dat.ileum$diagnosis, levels=c("no", "CD"))
dat.ileum$antibiotics = factor(dat.ileum$antibiotics, levels=c("false", "true"))
covars = model.matrix(~diagnosis+age+sex+antibiotics, data=dat.ileum)
counts.top.fam = counts.top.fam[rownames(covars),]

# Running mdine
ci.probs <- c(0.025, 0.975) #c(0.05, 0.95)
mn.family = mdine(counts.top.fam, covars, covars[,"diagnosisCD"], quant=ci.probs)

# Calcuating predicted proportions for each individual
predict0 = exp(covars[covars[,"diagnosisCD"]==0,]%*%mn.family$post_mean$beta)
predict1 = exp(covars[covars[,"diagnosisCD"]==1,]%*%mn.family$post_mean$beta)
comp0 = predict0/(rowSums(predict0)+1)
comp1 = predict1/(rowSums(predict1)+1)
# Averaging proportions over all subjects to scale size of nodes
mean_comp0 = colMeans(comp0)
mean_comp1 = colMeans(comp1)

#Reorder based on phyla so they appear close together in network
taxa_table = tax_table(phyloseq_family)
fam_ind = match(colnames(counts.top.fam)[-NCOL(counts.top.fam)], taxa_table[,"Rank5"])
phyla_for_fam = taxa_table[fam_ind,"Rank2"]
reorder = order(phyla_for_fam)
phyla_names = substr(as.vector(phyla_for_fam), 4, nchar(phyla_for_fam))

# Plotting networks
cc_network = plotNetworks(counts.top.fam, mn.family$post_mean$invsigma0, mn.family$post_mean$invsigma1,
                          ci0=mn.family$ci$invsigma0, ci1=mn.family$ci$invsigma1, 
                          lay=layout.mds,
                          lab0="Controls", lab1="Crohn's cases",
                          vertex.size0 = 3.5*(log(mean_comp0)-min(log(mean_comp0))+0.1), 
                          vertex.size1 = 3.5*(log(mean_comp1)-min(log(mean_comp1))+0.1),
                          col = as.numeric(factor(phyla_for_fam)), seed=3984, vertex.label.cex = 1,
                          phyla_names=phyla_names, side.legend=TRUE, cex.main=2)

# Network differences
cor_diff = abs(cov2cor(mn.family$post_mean$invsigma1)) - abs(cov2cor(mn.family$post_mean$invsigma0))
adj_diff = ci_to_adj(mn.family$ci$invsigma_diff[[1]],
                     mn.family$ci$invsigma_diff[[2]])
sig.diffs = adj_diff*cor_diff

diff_graph = graph.adjacency(sig.diffs, mode="undirected", weighted=TRUE, diag = FALSE)
E(diff_graph)$color[E(diff_graph)$weight>0] = "black"
E(diff_graph)$color[E(diff_graph)$weight<=0] = "black"
E(diff_graph)$lty[E(diff_graph)$weight>0] = 1
E(diff_graph)$lty[E(diff_graph)$weight<=0] = 2

scale_line_width = 30
labs = substr(colnames(counts.top.fam)[-NCOL(counts.top.fam)], 4, nchar(colnames(counts.top.fam)[-NCOL(counts.top.fam)]))

# Colour to show abundance differences
abund_diff_col = ifelse(mean_comp1-mean_comp0>0, "green", "red")

# Plotting differences in networks
layout(matrix(c(1,2,3,3), nrow=2, byrow=TRUE), heights=c(4, 1, 1), widths=c(4,1,2))
par(mai=rep(0.2, 4))
plot.igraph(diff_graph, layout=cc_network, vertex.label.cex=1,
            edge.width=abs(E(diff_graph)$weight)*scale_line_width, main="Network differences (CD minus control)",
            vertex.color=abund_diff_col, vertex.size=15,#8*(abs(log(mean_comp1)-log(mean_comp0))), 
            vertex.shape="circle",
            vertex.label=labs, vertex.label.color="black")
par(mai=c(0,0,0,0))
plot.new()
legend("center", legend=c("Higher in CD", "Lower in CD"), title="Family abundance",
       col="black", pt.bg=c("green", "red"), pch=21, pt.cex = 3, box.col = "white", y.intersp = 2)
plot.new()
legend("center", legend=c("CD abs. association stronger", "Control abs. association stronger"), 
       lty=c(1,2), col=c("black", "black"), ncol=2, cex=1, lwd = 2, box.col = "white")

