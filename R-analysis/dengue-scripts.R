# Exploratory Analysis of Dengue Transcriptomic Data
# Date: Monday 23rd November 2020
# Author: AS

# Aim: To gain insights into molecular mechanisms underlying host response to DENV infection, through:
# (1) visualising relationships between gene expression profiles of four patient populations; and
# (2) Identifying the genes that show most significant difference in expression between those populations.


# Import libraries
library(GEOquery)
library(ggplot2)
library(dplyr)
library(dendextend)
library(limma)
library(gplots)
library(EnhancedVolcano)


# Import data from GEO ----

# Load data (DENV expression dataset)
gse <- getGEO("GDS5093", GSEMatrix = TRUE)

# Create data frame
X <- Table(gse)
View(X)
class(X) # check type of object

# Convert GSE to ExpressionSet object
eset <- GDS2eSet(gse, do.log2 = TRUE)

# View variable names
featureNames(eset)

# View sample names
sampleNames(eset)

# Extract phenotypic state of each sample
pDat <- pData(eset)
View(pDat)


# Basic exploration ----

# Check column and row names
colnames(X)
rownames(X)

# Summarise gene expression data
summary(X)

# Summarise phenotypic data
summary(pDat)

# Set row names of X as gene names
geneNames <- as.character(X$IDENTIFIER)
X <- exprs(eset)
rownames(X) <- geneNames
dim(X) # check dimensions

# Check for any duplicated gene (row) names
length(which(duplicated(rownames(X))))

# Take average of probes corresponding to same gene
X <- avereps(X, ID = rownames(X))
dim(X)

# Basic exploration
plot(X)
plot(pDat$disease.state)
summary(X)
boxplot(X)
pairs(X[1:6, 2:6])


# Hierarchical cluster analysis ----

# Calculate distance between samples
tX <- t(X)
dist.tX <- dist(tX, method = "euclidean")
plot(dist.tX)

# Perform hierarchical clustering
comp.hclust <- hclust(dist.tX, method = "complete")

# Plot cluster dendrograms
plot(comp.hclust)


# Plot dendrogram ----

# Plot dendrogram
comp.dend <- as.dendrogram(comp.hclust)
plot(comp.dend)

# Colour samples according to disease status
col.comp.dend <- set(comp.dend, what = "branches_k_color", value = pDat$disease.state)

# Plot coloured dendrograms
plot(col.comp.dend)


# Relationship between genes ----

# Calculate standard deviation of genes across all samples
X.SD <- transform(X, SD = apply(X, 1, sd, na.rm = TRUE))

# Sort rows by highest SD
sort.sdX <- X.SD[with(X.SD, order(-SD)),]

# Extract top genes with highest SD
genes.X <- head(sort.sdX, 100)
hclust.X <- head(sort.sdX, 7528) # number genes used in original paper
top5000.X <- head(sort.sdX, 5000)
top1000.X <- head(sort.sdX, 1000)

# Remove SD column
genes.X <- genes.X[1:(length(genes.X)-1)]
hclust.X <- hclust.X[1:(length(hclust.X)-1)]
top5000.X <- top5000.X[1:(length(top5000.X)-1)]
top1000.X <- top1000.X[1:(length(top1000.X)-1)]


# Visualise as heatmap ----

# Convert named num to numeric matrix
genes.X <- as.matrix(genes.X)
hclust.X <- as.matrix(hclust.X)
top5000.X <- as.matrix(top5000.X)
top1000.X <- as.matrix(top1000.X)

# Plot heatmap (base R)
heatmap(x = genes.X, col = pDat$disease.state)

# Plot heatmap (gplots) for top 100 genes
heatmap.2(x = genes.X, trace = "none", ColSideColors = my_cols[pDat$disease.state],
          xlab = "Patient Samples (n=56)", ylab = "Genes (n=100)",
          key.title = "SD from mean", key.xlab = NA, key.ylab = NA,
          lhei=c(1,6), lwid=c(1,3.5), keysize=0.35, key.par = list(cex=0.5),
          labRow = FALSE, labCol = FALSE, margins = c(2, 2), col = greenred(75))

# Add legend
legend("topleft",
       legend = c("Convalescent", "Dengue Hemorrhagic Fever", "Dengue Fever", "Healthy Control"),
       fill = c("#0096FF", "#F8766D", "#E76CF3", "#00BA38"),
       border = NA, bty = "y",
       cex=0.5)

# Plot heatmap for all genes
heatmap.2(x = X, trace = "none", ColSideColors = my_cols[pDat$disease.state],
          xlab = "Patient Samples (n=56)", ylab = "Genes (n=31,654)",
          key.title = "SD from mean", key.xlab = NA, key.ylab = NA,
          lhei=c(1,4), lwid=c(1,3.5), keysize=0.35, key.par = list(cex=0.5),
          labRow = FALSE, labCol = FALSE, margins = c(2, 2), col = greenred(75))


# Plot heatmap for 7528 genes (as done in paper)
new_cols <- c("#9ECEFC", "#F6F237", "#272727", "#BD2D2C")
names(new_cols) <- c("Convalescent", "healthy control", "Dengue Fever", "Dengue Hemorrhagic Fever")

heatmap.2(x = hclust.X, trace = "none", ColSideColors = new_cols[pDat$disease.state],
          xlab = "Patient Samples (n=56)", ylab = "Genes (n=7,528)",
          key.title = "SD from mean", key.xlab = NA, key.ylab = NA,
          lhei=c(1,4), lwid=c(1,3.5), keysize=0.35, key.par = list(cex=0.5),
          labRow = FALSE, labCol = FALSE, margins = c(2, 2), col = greenred(75))


# Heatmap for top 5000 genes
heatmap.2(x = top5000.X, trace = "none", ColSideColors = my_cols[pDat$disease.state],
          xlab = "Patient Samples (n=56)", ylab = "Genes (n=5,000)",
          key.title = "SD from mean", key.xlab = NA, key.ylab = NA,
          lhei=c(1,4), lwid=c(1,3.5), keysize=0.35, key.par = list(cex=0.5),
          labRow = FALSE, labCol = FALSE, margins = c(2, 2), col = greenred(75))

# Heatmap for top 1000 genes
heatmap.2(x = top1000.X, trace = "none", ColSideColors = my_cols[pDat$disease.state],
          xlab = "Patient Samples (n=56)", ylab = "Genes (n=1,000)",
          key.title = "SD from mean", key.xlab = NA, key.ylab = NA,
          lhei=c(1,4), lwid=c(1,3.5), keysize=0.35, key.par = list(cex=0.5),
          labRow = FALSE, labCol = FALSE, margins = c(2, 2), col = greenred(75))


# Perform PCA ----

# Perform PCA on full Dengue dataset
Xpca <- prcomp(t(X), scale = TRUE)

# Get summary of PCA results
summary(Xpca)

# Get percentage variance for each PC
summ = summary(Xpca)
exp_var = summ$importance[2,] * 100

# Get cumulative variance for each PC
cum_var = summ$importance[3,] * 100
head(cum_var)

# Plot bar charts showing relative importance of each PC
barplot(exp_var)
barplot(cum_var)


# Plot PCA ----

# Extract PC scores
Xscores = Xpca$x

# Colour samples by disease state
my_cols <- c("#0096FF", "#F8766D", "#E76CF3", "#00BA38")
names(my_cols) <- c("Convalescent", "Dengue Hemorrhagic Fever", "Dengue Fever", "healthy control")

# Produce 2D PC plot
plot(Xscores[,1], Xscores[,2], xlab = "PC1 (16.5%)", ylab = "PC2 (8.2%)", pch = 19, 
     col = my_cols[pDat$disease.state], main = "PCA scores (per sample)")

# Add legend
legend("topright",
       legend = c("Convalescent", "Dengue Hemorrhagic Fever", "Dengue Fever", "Healthy Control"),
       fill = c("#0096FF", "#F8766D", "#E76CF3", "#00BA38"),
       border = NA, bty = "n")


# Investigate PC loadings ----

# Plot a biplot for full dataset
biplot(Xpca$x, Xpca$x, var.axes = TRUE, col = my_cols[pDat$disease.state], cex = 0.5,
       xlab = "PC1 (16.5%)", ylab = "PC2 (8.2%)", main = "PCA loadings (per sample)")

# Make biplot of top 100 genes
genes.pca <- prcomp(t(genes.X), scale = TRUE)
raw <- genes.pca$x[,1:2]
raw2 <- genes.pca$rotation[,1:2]
raw <- as.data.frame(raw)
raw2 <- as.data.frame(raw2)

plot(genes.pca$x[,1], genes.pca$x[,2], col = my_cols[pDat$disease.state])

biplot(genes.pca, cex = 0.6, col = my_cols[pDat$disease.state], 
       xlab = "PC1 (16.5%)", ylab = "PC2 (8.2%)", main = "PCA loadings (per gene)")

p <- (ggplot(raw2, aes(PC1, PC2, label = rownames(raw2))) +
        geom_text() +
        theme_bw() +
        xlab("PC1 (16.5%)") + ylab("PC2 (8.2%)") + 
        ggtitle("PCA loadings (per gene)") +
        theme(plot.title = element_text(hjust = 0.5, size = 14)))

q <- ggplot(as.data.frame(raw), aes(PC1, PC2, colour = my_cols[pDat$disease.state])) +
  geom_point() +
  scale_fill_manual(name="Disease Status")

p + z$layers


# Task 9: Plot PCA scores without transforming or scaling ----

# Perform PCA
Xpca2 <- prcomp(X, scale = FALSE)

# Get percentage variance for each PC
summ2 = summary(Xpca2)
exp_var2 = summ2$importance[2,] * 100

# Get cumulative variance for each PC
cum_var2 = summ2$importance[3,] * 100
head(cum_var2)

# Plot bar charts showing relative importance of each PC
barplot(exp_var2) # differences a lot more extreme between PC1 and rest of PCs.
barplot(cum_var2) # less variation between PCs.

# Plot PCA
Xscores2 = Xpca2$x
plot(Xscores2[,1], Xscores2[,2], xlab="PC1 (94.3%)", ylab= "PC2 (1.2%)", pch=19)


# Volcano plot to show significant differentially expressed genes ----

# Create design matrix
f <- factor(pDat$disease.state, levels = c("Convalescent", "Dengue Fever", "Dengue Hemorrhagic Fever", "healthy control"))
design <- model.matrix(~0+f)
colnames(design) <- c("Convalescent", "Dengue.Fever", "Dengue.Hemorrhagic.Fever", "Healthy.Control")

# Fit the model
fit <- lmFit(eset, design)

# Make pair-wise comparisons between all disease state groups
contrast.matrix <- makeContrasts(Convalescent-Dengue.Fever,
                                 Convalescent-Dengue.Hemorrhagic.Fever,
                                 Convalescent-Healthy.Control,
                                 Dengue.Fever-Dengue.Hemorrhagic.Fever,
                                 Dengue.Fever-Healthy.Control,
                                 Dengue.Hemorrhagic.Fever-Healthy.Control,
                                 levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
plotSA(fit2, main = "Final model: Mean-variance trend") # each dot is a gene

# Get table of top-ranked genes
tt5 <- topTable(fit2, coef = 5, adjust = "BH")

# Get ranked genes for all disease states
tt.all <- topTable(fit2, adjust = "BH", number = 54715)

# Look for genes which vary between each different disease state (i.e. top 30 genes)
top30.tt <- topTable(fit2, number = 30)

# Get outcome of each hypothesis test
results <- decideTests(fit2)

# Get summary of number of significant genes in each comparison
summary(results)
results.df <- as.data.frame(summary(results))
results.df$Var2 <- as.factor(results.df$Var2)
results.df$Var1 <- as.factor(results.df$Var1)
str(results.df)

# Plot results
colnames(results.df) <- c("expression.levels", "disease.state", "frequency")
results.df <- subset(results.df,  ! paste(expression.levels) %in% c("NotSig") )

ggplot(results.df, aes(x=disease.state, y=frequency, fill=expression.levels)) +
  geom_bar(stat="identity", position="dodge") +
  theme_bw() +
  xlab("Disease State") + ylab("Frequency") +
  scale_x_discrete(labels=c("C-DF", "C-DHF", "C-HC", "DF-DHF", "DF-HC", "DHF-HC")) +
  ggtitle("Total expression level changes of genes across disease states") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14,face = "bold")) +
  scale_fill_manual(name="Expression Levels", values = c("Down" = "black", "Up" = "red"))

# Plot mean difference plot
plotMD(fit2, column = 6)

# Plot volcano plot
volcanoplot(fit2)


# Plot enhanced volcano ----

# Get ranked genes between healthy controls and patients with DHF
tt.1 <- topTable(fit2, coef = 1, adjust = "BH", number = 54715) # C vs DF
tt.2 <- topTable(fit2, coef = 2, adjust = "BH", number = 54715) # C vs DHF
tt.3 <- topTable(fit2, coef = 3, adjust = "BH", number = 54715) # C vs HC
tt.4 <- topTable(fit2, coef = 4, adjust = "BH", number = 54715) # DF vs DHC
tt.5 <- topTable(fit2, coef = 5, adjust = "BH", number = 54715) # DF vs HC
tt.6 <- topTable(fit2, coef = 6, adjust = "BH", number = 54715) # HC vs DHF

# Plot basic volcano plot for HC vs DHF
EnhancedVolcano(tt.6,
                lab = tt.6$Gene.symbol,
                x = "logFC",
                y = "adj.P.Val",
                title = "Healthy Controls vs Dengue Hemorrhagic Fever",
                labSize = 4,
                legendLabels=c("Not significant", bquote(~Log[2]~ "fold change"),"p-value", 
                               bquote("p-value & " ~Log[2]~ "fold change")),
                legendIconSize = 5,
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                widthConnectors = 0.8)


# Plot basic volcano plot for C vs DF
EnhancedVolcano(tt.1,
                lab = tt.1$Gene.symbol,
                x = "logFC",
                y = "adj.P.Val",
                title = "Convalescent vs Dengue Fever",
                labSize = 4,
                legendLabels=c("Not significant", bquote(~Log[2]~ "fold change"),"p-value", 
                               bquote("p-value & " ~Log[2]~ "fold change")),
                legendIconSize = 5,
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                widthConnectors = 0.8)

# Plot basic volcano plot for C vs DHF
EnhancedVolcano(tt.2,
                lab = tt.2$Gene.symbol,
                x = "logFC",
                y = "adj.P.Val",
                title = "Convalenscent vs Dengue Hemorrhagic Fever",
                labSize = 4,
                legendLabels=c("Not significant", bquote(~Log[2]~ "fold change"),"p-value", 
                               bquote("p-value & " ~Log[2]~ "fold change")),
                legendIconSize = 5,
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                widthConnectors = 0.8)

# Plot basic volcano plot for C vs HC
EnhancedVolcano(tt.3,
                lab = tt.3$Gene.symbol,
                x = "logFC",
                y = "adj.P.Val",
                title = "Convalescent vs Healthy Control",
                labSize = 4,
                legendLabels=c("Not significant", bquote(~Log[2]~ "fold change"),"p-value", 
                               bquote("p-value & " ~Log[2]~ "fold change")),
                legendIconSize = 5,
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                widthConnectors = 0.8)

# Plot basic volcano plot for DF vs DHC
EnhancedVolcano(tt.4,
                lab = tt.4$Gene.symbol,
                x = "logFC",
                y = "adj.P.Val",
                title = "Dengue Fever vs Dengue Hemorrhagic Fever",
                labSize = 4,
                legendLabels=c("Not significant", bquote(~Log[2]~ "fold change"),"p-value", 
                               bquote("p-value & " ~Log[2]~ "fold change")),
                legendIconSize = 5,
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                widthConnectors = 0.8)

# Plot basic volcano plot for DF vs HC
EnhancedVolcano(tt.5,
                lab = tt.5$Gene.symbol,
                x = "logFC",
                y = "adj.P.Val",
                title = "Dengue Fever vs Healthy Control",
                labSize = 4,
                legendLabels=c("Not significant", bquote(~Log[2]~ "fold change"),"p-value", 
                               bquote("p-value & " ~Log[2]~ "fold change")),
                legendIconSize = 5,
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                widthConnectors = 0.8)
