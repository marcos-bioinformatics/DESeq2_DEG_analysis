## Author: Marcos Elizalde Horcada - marcos.elizaldeh@gmail.com
## Date: January 2021

#--------------------------------------------------------------------#
#------------- SETTING UP THE WORKSPACE AND DATA LOADING ------------# 
#--------------------------------------------------------------------#

# Load packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(DESeq2)
library("readxl")
library(dplyr)
library(tidyr)
library(stringr)

# Load the raw counts dataframe and the experimental design
raw_counts <- read.csv("010920_genes_annot.csv", 
                       stringsAsFactors = FALSE)

# Tidy some columns
genes <- raw_counts$ENSEMBL
rownames(raw_counts) <- genes
raw_counts <- raw_counts[,-1]

# Resolve problem with row duplicates
duplicates <- duplicated(raw_counts$gene_id)
table(duplicates)
nams <- rownames(counts)
rownames(counts) <- make.names(nams, unique=TRUE)

# Experimental design should contain sample names in the same order as in
# the raw counts data frame, as input for DESeq2
design <- read_excel("Experimental_design.xlsx")
design$Type <- as.factor(design$Type)
design$Sex <- as.factor(design$Sex)
design$Batch <- as.factor(design$Batch)
str(design)

keep <- colnames(raw_counts) %in% design$Sample
table(keep)
raw_counts <- raw_counts[,keep] #67 individuals

# Create the DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = design,
  design = ~ Sex + Type) #Check for design type

#--------------------------------------------------------------------#
#------------ EXPLORATORY DATA ANALYSIS AND VISUALIZATION -----------# 
#--------------------------------------------------------------------#

# First we filter removing rows with no counts or a selected value of counts.
# For RNA-seq counts, the expected variance grows with the mean, so there
# is no homokedasticity. vst transformation looks at the trend between
# variance and the mean in the data and transforms the data so this trend
# is removed. This transformation is only for visualization and has no 
# effect on the DEG analysis. 

# Plot dispersions
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# Pre-filtering the dataset
nrow(dds)
keep <- rowSums(counts(dds) >= 1) >= 4
dds <- dds[keep,]
nrow(dds)

# Compute VST transformation
vsd <- vst(dds, blind = FALSE) #Set to TRUE for unsupervised transformation
head(assay(vsd))

# Check for Euclidean distances between samples
sampleDists <- dist(t(assay(vsd)))
sampleDists     

# Heatmap of Euclidean sample distances
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(design$Type, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# PCA Plot with plotPCA
plotPCA(vsd, intgroup = c("Sample")) #CHECK FOR GROUPS

# Generalized PCA with glmpca
library("glmpca")
gpca <- glmpca(counts(dds), L=4)
gpca.dat <- gpca$factors #SPECIFY FACTORS FROM EXP DESIGN
gpca.dat$Type <- dds$Type
gpca.dat$Sex <- dds$Sex

ggplot(gpca.dat, aes(x = dim1, y = dim2, color = Type, shape = Sex)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")

ggplot(gpca.dat, aes(x = dim1, y = dim3, color = Type, shape = Sex)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")

ggplot(gpca.dat, aes(x = dim1, y = dim4, color = Type, shape = Sex)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")

ggplot(gpca.dat, aes(x = dim2, y = dim3, color = Type, shape = Sex)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")

ggplot(gpca.dat, aes(x = dim2, y = dim4, color = Type, shape = Sex)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")

ggplot(gpca.dat, aes(x = dim3, y = dim4, color = Type, shape = Sex)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")

# PCA Plot with FactoMineR and hierarchical clustering
library("FactoMineR")
library("factoextra")
# Increase max. overlaps
options(ggrepel.max.overlaps = Inf)

# Load the gene expression matrix and experimental design
# Principal Component Analysis (PCA)
res.pca <- PCA(t(raw_counts),
               graph = FALSE,
               scale.unit = TRUE)

get_eigenvalue(res.pca)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 70), main = "")

pdf("eigenvalues_patients.pdf")
fviz_pca_ind(res.pca, 
             col.ind = design$Type, 
             pointsize=2, 
             pointshape=19,
             palette="lancet",
             fill="black",
             label = "none",
             repel = TRUE, 
             addEllipses = TRUE,
             ellipse.type = "confidence",
             legend.title="Patients",
             title="Principal Components Analysis",
             show_legend=FALSE,show_guide=FALSE)
dev.off()

# Contributions of variables to PC1
print(fviz_contrib(res.pca, choice = "var", axes = 1, top = 20))
# Contributions of variables to PC2
print(fviz_contrib(res.pca, choice = "var", axes = 2, top = 20))
# Contributions of samples to PC1
print(fviz_contrib(res.pca, choice = "ind", axes = 1))
# Contributions of samples to PC2
print(fviz_contrib(res.pca, choice = "ind", axes = 2))

# Compute Hierarchical Clustering on Principal Components
res.hcpc <- HCPC(res.pca, graph=FALSE)    

pdf("cluster.pdf")
fviz_dend(res.hcpc,k=3,
          cex = 0.7,       # Label size
          palette = "lancet", # Color palette see ?ggpubr::ggpar
          rect = TRUE,     # Add rectangle around groups
          rect_fill = TRUE, 
          rect_border = "jco",           
          type="rectangle",
)
dev.off()

# Show sample clustering map
pdf("cluster_ind.pdf")
fviz_cluster(res.hcpc,ellipse = TRUE,
             repel = TRUE,               # Avoid label overlapping
             show.clust.cent = TRUE,     # Show cluster centers
             palette = "lancet",            # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Sample Clustering"
)
dev.off()

# Show cluster dendrogram with scale
pdf("hierarchical_clustering_patients.pdf")
plot(res.hcpc, choice ="tree")
dev.off()


#--------------------------------------------------------------------#
#----------------- DIFFERENTIAL EXPRESSION ANALYSIS -----------------# 
#--------------------------------------------------------------------#

# Independent filtering increases detection power for high-throughput
# experiments (Bourgon et al. 2010)
# Run the differential expression pipeline
dds <- DESeq(dds)
resultsNames(dds)

# Save normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(as.data.frame(normalized_counts),
          file="normalized_counts.csv")
#Calling results without any arguments will extract the estimated log2 
#fold changes and p values for the last variable in the design formula. 
#If there are more than 2 levels for this variable, results will extract 
#the results table for a comparison of the last level over the first level.
#res <- results(dds)

# Specify baseline, otherwise it does following alphabetical order
res <- results(dds, contrast = c("Type", "Treatment", "Control"),
               alpha = 0.05)

# Check number of genes specifying conditions
sum(res$padj < 0.05 & (res$log2FoldChange > abs(2)), na.rm=TRUE)

# Order dataframe of results by pvalue
resOrdered <- res[order(res$padj),]
write.csv(as.data.frame(resOrdered), 
          file="condition_treated_0.05.csv")

# Extract results for FDR 0.05 and LFC 0,59
res <- results(dds, contrast = c("Type", "Treatment", "Control"),
               alpha = 0.05,
               lfcThreshold = 0.59)

# Order dataframe of results by pvalue
resOrdered <- res[order(res$padj),]
write.csv(as.data.frame(resOrdered), 
          file="condition_treated_0.05_LFC059.csv")

# Extract results for FDR 0.05 and LFC 1
res <- results(dds, contrast = c("Type", "Treatment", "Control"),
               alpha = 0.05,
               lfcThreshold = 1)

# Order dataframe of results by pvalue
resOrdered <- res[order(res$padj),]
write.csv(as.data.frame(resOrdered), 
          file="condition_treated_0.05_LFC1.csv")

# Bonferroni correction
bonferroni <- read.csv("condition_treated_0.05.csv")
bonferroni$bonferroni <- p.adjust(bonferroni$pvalue, 
                                  method = "bonferroni", 
                                  n = length(bonferroni$pvalue))

write.csv(as.data.frame(bonferroni),
          file="condition_treated_0.05_Bonferroni.csv")

# Apply shrinkage for plotting
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="Type_Treatment_vs_Control", type="apeglm")

#MA Plots
plotMA(res, ylim=c(-4,4))
plotMA(resLFC, ylim=c(-4,4))


