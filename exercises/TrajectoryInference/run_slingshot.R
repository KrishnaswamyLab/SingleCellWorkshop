install.packages("BiocManager")
BiocManager::install("slingshot")

library(slingshot)
library(SingleCellExperiment)

# Load counts, PHATE, and cluster labels
counts <- as.data.frame(read.csv('~/scRNAseq/Treutlein.expression.csv'))
phate = as.matrix(read.csv('/home/dan/scRNAseq/Treutlein.phate.csv'))
pclusters = read.csv('/home/dan/scRNAseq/Treutlein.phate_clusters.csv')

# Create SingleCellExperiment 
# How am I actually supposed to have the column names not be passed as genes??
sim <- SingleCellExperiment(assays = List(counts = t(as.matrix(counts[,2:2001]))))

# Add dim red data and clusters to SCE
reducedDims(sim) <- SimpleList(PHATE=phate)
colData(sim)$pclusters <- pclusters[,2]

# Optional, plot data and clusters
library(RColorBrewer)
plot(phate, col = brewer.pal(5,"Set1")[sim$pclusters], pch=16, asp = 1)

# Do Slingshot
sce <- slingshot(sim, clusterLabels = 'pclusters', reducedDim = 'PHATE')

summary(sce$slingPseudotime_3)

# Plot Slingshot
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

# For some reason not all points are plotted here
plot(reducedDims(sce)$PHATE, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

# This plots the 'scaffold'
plot(reducedDims(sce)$PHATE, col = brewer.pal(5,"Set1")[sim$pclusters], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')

#You can get the orderings of the points from these variables
sce$slingPseudotime_1
sce$slingPseudotime_2
