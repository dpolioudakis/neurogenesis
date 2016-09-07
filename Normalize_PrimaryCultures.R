# Damon Polioudakis
# 2016-07-21
# Normalize array expression data from Andreas’ VZ/CP primary cultures
# differentiation time course

### TODO
# Add bootstrapping regression for regressing out covariates

################################################################################

rm(list=ls())
sessionInfo()

require(ggplot2)
require(reshape2)
require(WGCNA)
# require(boot)

### Load data and assign variables

## Load data

# hNPs
# Primary cultures from human fetal VZ/CP
exDF <- read.csv("../data/primaryCultures_exprdata.csv")
annotDF <- read.csv("../metadata/primaryCultures_annot.csv")
metDF <- read.csv("../metadata/primaryCultures_metadata.csv")

# Out graphs
outGraphs <- "../analysis/graphs/Normalize_PrimaryCultures_"
graphsTitle <- "Normalize_PrimaryCultures.R"

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text=element_text(size=18)))
theme_update(plot.title = element_text(size = 16))
################################################################################

### Format

## Extract brain region and differentiation week from sample name

df <- metDF[match(row.names(metDF), colnames(exDF[ ,-(1:3)])), ]
prmWkZoneDF <- data.frame(WEEK = gsub("^T([0-9]).*", "\\1", row.names(df), perl = TRUE))
prmWkZoneDF$ZONE <- gsub("^T[0-9]_(.*)_.*", "\\1", row.names(df), perl = TRUE)
prmWkZoneDF <- prmWkZoneDF[-c(1:3), ]

### Process Primary Cultures array data

## Graph data pre-normalization

pdf(file = paste0(outGraphs, "QC.pdf"))
par(mfrow = c(4, 1))

boxplot(exDF, range = 0, col = as.factor(prmWkZoneDF$WEEK)
  , ylab = "Log2 expression ?"
  , main = "Boxplot log2 expression ?")

plot(density(exDF[ ,1]), col = as.factor(prmWkZoneDF$WEEK)[1]
  , main = "Density log2 expression ?")
for (i in 2:dim(exDF)[2]) {
  lines(density(exDF[ ,i]), col = as.factor(prmWkZoneDF$WEEK)[i])
}

# MDS colored by differentiation week
mdsG = cmdscale(dist(t(exDF)), eig = TRUE)
varExpl <- round(mdsG$eig / sum(mdsG$eig), 3)
plot(mdsG$points, col = as.numeric(as.factor(prmWkZoneDF$WEEK)), pch = 19
  , xlab = paste("Dimension 1 (", 100*varExpl[1], "%)", sep = "")
  , ylab = paste("Dimension 2 (", 100*varExpl[2], "%)", sep = "")
  , main = "MDS plot - colored by differentiation week")

# MDS colored by VZ/CP
mdsG = cmdscale(dist(t(exDF)), eig = TRUE)
varExpl <- round(mdsG$eig / sum(mdsG$eig), 3)
plot(mdsG$points, col = as.numeric(as.factor(prmWkZoneDF$ZONE)), pch = 19
  , xlab = paste("Dimension 1 (", 100*varExpl[1], "%)", sep = "")
  , ylab = paste("Dimension 2 (", 100*varExpl[2], "%)", sep = "")
  , main = "MDS plot - colored by VZ/CP")
dev.off()