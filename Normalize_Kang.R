# Damon Polioudakis
# 2016-09-06
# Normalize array expression data from Kang brain stages

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

# Kang
load("../orig.data/InVivoData/WGCNAinput_SestanBrain.RData")
exDF = t(datExpr)
rm(datExpr)
# sampleKey is ProcessedKangMetaData.csv
metDF <- sampleKey

# Out graphs
outGraphs <- "../analysis/graphs/Normalize_Kang_"
graphsTitle <- "Normalize_Kang.R"

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text=element_text(size=18)))
theme_update(plot.title = element_text(size = 16))
################################################################################

### Format

## Metadata character class to factor
metDF[sapply(metDF, is.character)] <- lapply(metDF[sapply(metDF, is.character)]
  , as.factor)

### Process Primary Cultures array data

## Graph data pre-normalization

pdf(file = paste0(outGraphs, "QC.pdf"))
par(mfrow = c(3, 1))

boxplot(exDF, range = 0, col = as.factor(metDF$Stage)
  , ylab = "Log2 expression ?"
  , main = "Boxplot log2 expression ?")

plot(density(log2(exDF[ ,1]+1)), col = as.factor(metDF$Stage)[1]
  , main = "Density log2 expression ?")
for (i in 2:dim(exDF)[2]) {
  lines(density(log2(exDF[ ,i]+1)), col = as.factor(metDF$Stage)[i])
}

mdsG = cmdscale(dist(t(exDF)), eig = TRUE)
varExpl <- round(mdsG$eig / sum(mdsG$eig), 3)
plot(mdsG$points, col = as.numeric(as.factor(metDF$Stage)), pch = 19
  , xlab = paste("Dimension 1 (", 100*varExpl[1], "%)", sep = "")
  , ylab = paste("Dimension 2 (", 100*varExpl[2], "%)", sep = "")
  , main = "MDS plot - colored by stage")
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

pairsDat <- metDF[c("BrainCode", "Stage", "PMI", "pH", "RIN"
  , "Hemisphere", "Sex")]

pairsDat <- data.frame(metDF[c("Stage", "PMI", "pH", "RIN")])

cond <- labels2colors(metDF$Stage)  ## colors

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

### Regress out covariates

covDF <- data.frame(metDF[c("Stage", "PMI", "RIN", "Sex")])
covDF$Sex <- as.numeric(covDF$Sex)

# Regress out confounding variables
RegressCovariates <- function (exM, covDF) {
  X = model.matrix(~ ., data = covDF)
  Y = exM
  beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
  b = as.data.frame(t(beta))
  to_regress = (as.matrix(X[,3:5]) %*% (as.matrix(beta[3:5,])))
  exRegCovM = exM - t(to_regress)
  return(exRegCovM)
}
nmExDF <- RegressCovariates(exDF, covDF)
quantile(exDF[ ,1], c(0,0.025,0.25,0.5,0.75,0.975,1))
quantile(nmExDF[ ,1], c(0,0.025,0.25,0.5,0.75,0.975,1))

## PCA of expression and technical covariates after regressing out

# Centers the mean of all genes - this means the PCA gives us the eigenvectors
# of the geneXgene covariance matrix, allowing us to assess the proportion of
# variance each component contributes to the data
meanCenteredM <- t(scale(t(nmExDF), scale=F))
# Run PCA
pCdat <- prcomp(meanCenteredM, center=F);
topPCs <- pCdat$rotation[,1:5];
# Calculate variance explained by each PC
varExp <- (pCdat$sdev)^2 / sum(pCdat$sdev^2)
topVar <- varExp[1:5]
colnames(topPCs) <- paste("Expression\n", colnames(topPCs)
  , " (", signif(100 * topVar[1:5], 2), "%)", sep = "")

pdf(paste0(outGraphs, "CovariatesMatrix_Regressed.pdf"), height = 12, width = 12)
pairs(cbind(topPCs, pairsDat), col = cond, pch = 19
  , upper.panel = panel.cor
  , main = "Covariates and MaxQuant Comparison -- |Spearman's rho| correlation values")
dev.off()

save(nmExDF, file = "../analysis/Kang_Expression_RgCv.RData")
