# ==============================
# STEP 1 — DATA OVERVIEW
# ==============================

library(tidyverse)

# Load dataset
bvals <- read.csv("../Data/GSE279988_bvals.csv.gz", row.names = 1)

# Check dimensions
dim(bvals)

# Preview data
head(bvals[,1:5])

# Summary statistics
summary(as.vector(as.matrix(bvals)))
# Check missing values
sum(is.na(bvals))
# ==============================
# HISTOGRAM — METHYLATION DISTRIBUTION
# ==============================

# Convert matrix to vector
all_values <- as.vector(as.matrix(bvals))

# Save plot
png("../Results/01_histogram.png", width = 800, height = 600)

hist(all_values,
     breaks = 50,
     main = "Distribution of DNA Methylation Values",
     xlab = "Beta values",
     col = "lightblue")

dev.off()
# ==============================
# PCA — SAMPLE SIMILARITY
# ==============================

# Transpose data (samples as rows)
bvals_t <- t(bvals)

# Run PCA
pca <- prcomp(bvals_t, scale. = TRUE)

# Save plot
png("../Results/02_PCA.png", width = 800, height = 600)

plot(pca$x[,1], pca$x[,2],
     main = "PCA of Samples",
     xlab = "PC1",
     ylab = "PC2",
     col = "blue",
     pch = 19)

# Close file
dev.off()
# ==============================
# HEATMAP — TOP VARIABLE CpGs
# ==============================

library(pheatmap)

# Compute variance of each CpG
var_sites <- apply(bvals, 1, var)

# Select top 100 most variable CpGs
top_sites <- names(sort(var_sites, decreasing = TRUE))[1:100]

# Save heatmap
png("../Results/03_heatmap.png", width = 800, height = 800)

pheatmap(bvals[top_sites, ],
         show_rownames = FALSE,
         main = "Top 100 Variable CpG Sites")

dev.off()
# ==============================
# DIFFERENTIAL METHYLATION
# ==============================

library(limma)

# Create groups (since no labels available)
group <- factor(c(rep("G1", 12), rep("G2", 12)))

# Design matrix
design <- model.matrix(~ group)

# Fit model
fit <- lmFit(bvals, design)
fit <- eBayes(fit)

# Extract all CpGs
results <- topTable(fit, coef = 2, number = Inf)

# Preview
head(results)
# ==============================
# VOLCANO PLOT — HIGHLIGHT CpGs
# ==============================

library(ggplot2)

# Define strong CpGs
results$strong <- -log10(results$P.Value) > 4 & abs(results$logFC) > 0.05

# Save plot
png("../Results/04_volcano.png", width = 800, height = 600)

ggplot(results, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(color = "grey70") +
  
  # Highlight important CpGs
  geom_point(data = subset(results, strong),
             color = "red", size = 1.2) +
  
  geom_hline(yintercept = 4, linetype = "dashed") +
  geom_vline(xintercept = c(-0.05, 0.05), linetype = "dashed") +
  
  labs(title = "Volcano Plot — Significant CpGs",
       x = "Log Fold Change",
       y = "-log10(P-value)") +
  theme_minimal()

dev.off()

# ==============================
# SUMMARY — CpG SIGNIFICANCE
# ==============================

# Define categories
results$category <- ifelse(results$adj.P.Val < 0.05, "Significant", "Not Significant")

# Count
table_counts <- table(results$category)

# Save plot
png("../Results/05_CpG_summary.png", width = 800, height = 600)

barplot(table_counts,
        col = c("grey", "red"),
        main = "Distribution of Significant CpG Sites",
        ylab = "Number of CpGs")

dev.off()
# ==============================
# BOXPLOT — METHYLATION PER GROUP
# ==============================

# Create groups (same as before)
group <- factor(c(rep("G1", 12), rep("G2", 12)))

# Convert data to long format
library(tidyr)
library(dplyr)

df_long <- as.data.frame(bvals) %>%
  mutate(CpG = rownames(bvals)) %>%
  pivot_longer(-CpG, names_to = "Sample", values_to = "Methylation")

# Add group information
df_long$Group <- rep(group, each = nrow(bvals))

# Save plot
png("../Results/06_boxplot_groups.png", width = 800, height = 600)

boxplot(Methylation ~ Group, data = df_long,
        col = c("lightblue", "lightgreen"),
        main = "Methylation Distribution per Group",
        ylab = "Beta values")

dev.off()
# ==============================
# DENSITY PLOT — METHYLATION
# ==============================

# Use previous vector
# all_values <- as.vector(as.matrix(bvals))

png("../Results/07_density.png", width = 800, height = 600)

plot(density(all_values),
     main = "Density of DNA Methylation Values",
     xlab = "Beta values",
     col = "blue",
     lwd = 2)

dev.off()