# Damon Polioudakis
# 2016-09-06
# Normalize array expression data from Miller brain regions

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

# Miller
load("../orig.data/LCMDE/AllenLCM.Rdata")
exDF = AllenLCM$datExpr
metDF = AllenLCM$datTraits
Zones = read.csv("../orig.data/LCMDE/LCM_Zones_CPio.csv")

# Out graphs
outGraphs <- "../analysis/graphs/Normalize_Miller_"
graphsTitle <- "Normalize_Miller.R"

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text=element_text(size=18)))
theme_update(plot.title = element_text(size = 16))
################################################################################

### Format

## Select cortical layers only

ZonesMiller = metDF$structure_acronym

metDF$Labels = rep(NA, nrow(metDF))
for (Label in Zones$Label)
  metDF[grep(Label, ZonesMiller), "Labels"] = Label

Zoneindex = which(! is.na(metDF$Labels))

exDF = exDF[, Zoneindex]
metDF = metDF[Zoneindex, ]

## Metadata character class to factor
metDF[sapply(metDF, is.character)] <- lapply(metDF[sapply(metDF, is.character)]
  , as.factor)

## Format zone codes
millerZones <- metDF[match(metDF$well_id, colnames(exDF)), ]$Labels
millerZones <- gsub(".*SG.*", "SG", millerZones)
millerZones <- gsub(".*MZ.*", "MZ", millerZones)
millerZones <- gsub(".*CP.*", "CP", millerZones)
millerZones <- gsub(".*SP.*", "SP", millerZones)
millerZones <- gsub(".*IZ.*", "IZ", millerZones)
millerZones <- gsub(".*SZ.*", "SZ", millerZones)
millerZones <- gsub(".*VZ.*", "VZ", millerZones)
millerZones <- as.factor(millerZones)

metDF$Labels <- millerZones
################################################################################

## Graph data pre-normalization

pdf(file = paste0(outGraphs, "QC.pdf"))
par(mfrow = c(3, 1))

boxplot(exDF, range = 0, col = as.factor(metDF$Labels)
  , ylab = "Log2 expression ?"
  , main = "Boxplot log2 expression ?")

plot(density(log2(exDF[ ,1]+1)), col = as.factor(metDF$Labels)[1]
  , main = "Density log2 expression ?")
for (i in 2:dim(exDF)[2]) {
  lines(density(log2(exDF[ ,i]+1)), col = as.factor(metDF$Labels)[i])
}

mdsG = cmdscale(dist(t(exDF)), eig = TRUE)
varExpl <- round(mdsG$eig / sum(mdsG$eig), 3)
plot(mdsG$points, col = as.numeric(as.factor(metDF$Labels)), pch = 19
  , xlab = paste("Dimension 1 (", 100*varExpl[1], "%)", sep = "")
  , ylab = paste("Dimension 2 (", 100*varExpl[2], "%)", sep = "")
  , main = "MDS plot - colored by region")
dev.off()


## PCA of expression and technical covariates

# Centers the mean of all genes - this means the PCA gives us the eigenvectors
# of the geneXgene covariance matrix, allowing us to assess the proportion of
# variance each component contributes to the data
meanCenteredM <- t(scale(t(exDF), scale=F))
# Run PCA
pCdat <- prcomp(meanCenteredM, center=F);
topPCs <- pCdat$rotation[,1:5];
# Calculate variance explained by each PC
varExp <- (pCdat$sdev)^2 / sum(pCdat$sdev^2)
topVar <- varExp[1:5]
colnames(topPCs) <- paste("Expression\n", colnames(topPCs)
  , " (", signif(100 * topVar[1:5], 2), "%)", sep = "")


## Correlation matrix of Luis VZ CP RNAseq expression PCs, seq stats PCs, and
# metadata

pairsDat <- data.frame(metDF[c("well_id", "BrainLabel", "Sex", "Ancestry", "AvgRIN", "Hemisphere", "Labels")])

cond <- labels2colors(metDF$Labels)  ## colors

# Useful function for comparing multivariate data
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  if (class(x) == "numeric" & class(y) == "numeric") {
    r <- abs(cor(x, y, use = "pairwise.complete.obs", method = "spearman"))
  } else {
    lmout <- lm(y~x)
    r <- sqrt(summary(lmout)$adj.r.squared)
  }
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

pdf(paste0(outGraphs, "CovariatesMatrix.pdf"), height = 12, width = 12)
pairs(cbind(topPCs, pairsDat), col = cond, pch = 19
  , upper.panel = panel.cor
  , main = "Covariates and MaxQuant Comparison -- |Spearman's rho| correlation values")
dev.off()
################################################################################