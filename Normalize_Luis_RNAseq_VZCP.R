# Damon Polioudakis
# 2016-07-21
# Normalize Luis' bulk RNAseq from experiments for ATAC work

### TODO
# Add bootstrapping regression for regressing out covariates

################################################################################

rm(list=ls())
sessionInfo()

require(ggplot2)
require(reshape2)
require(xlsx)
require(cqn)
require(WGCNA)
# require(boot)

### Load data and assign variables

## Load data

# Luis RNAseq fetal brain VZ CP
vzcpExDatDF <- read.csv("../data/bulk_VZ_CP_from_ATAC/Exprs_HTSCexon.csv"
                        , row.names = 1)
# Picard Sequencing Statistics - bulk RNAseq
vzcpPicStatsDF <- read.csv("../metadata/PicardToolsQC.csv", fill = TRUE
                           , header = TRUE)

# Luis RNAseq
metDatDF <- read.csv("../metadata/VZCP_sampleinfo.csv", header = TRUE)

# Gene lengths and GC content for Union Exon model
load("../source/ENSEMBLhg19_UnionAnno.rda")

# Out graphs
outGraphs <- "../analysis/graphs/Normalize_Luis_RNAseq_VZCP_"
graphsTitle <- "Normalize_Luis_RNAseq_VZCP.R"

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text=element_text(size=18)))
theme_update(plot.title = element_text(size = 16))
################################################################################

### Process Luis RNAseq


## filter htseqcounts 

pres = apply(vzcpExDatDF > 5, 1, sum) 
idx = which(pres > 0.4*dim(vzcpExDatDF)[2]) ## exp > 10 in 80% of samples
vzcpFtExDatDF = vzcpExDatDF[idx, ]
# vzcpFtExDatDF <- vzcpExDatDF


## Graph data pre-normalization

pdf(file = paste0(outGraphs, "QC_NoNorm.pdf"))
par(mfrow = c(3, 1))

boxplot(log2(vzcpFtExDatDF+1), range = 0, col = as.factor(metDatDF$ExpCondition)
        , ylab = "Log2 (expression + 1)"
        , main = "Boxplot HTSeq log2(Counts + 1)")

plot(density(log2(vzcpFtExDatDF[ ,1]+1)), col = as.factor(metDatDF$ExpCondition)[1]
     , main = "Density HTSeq log2(Counts + 1)")
for (i in 2:dim(vzcpFtExDatDF)[2]) {
  lines(density(log2(vzcpFtExDatDF[ ,i]+1)), col = as.factor(metDatDF$ExpCondition)[i])
}

mdsG = cmdscale(dist(t(vzcpFtExDatDF)), eig = TRUE)
varExpl <- round(mdsG$eig / sum(mdsG$eig), 3)
plot(mdsG$points, col = as.numeric(as.factor(metDatDF$ExpCondition)), pch = 19
     , xlab = paste("Dimension 1 (", 100*varExpl[1], "%)", sep = "")
     , ylab = paste("Dimension 2 (", 100*varExpl[2], "%)", sep = "")
     , main = "MDS plot post CQN for length and GC")
dev.off()


## CQN Normalization

# Use CQN to normalize for GC content and gene length
RunCQN <- function (exprDatDF, lengthGcDF) {
  # Remove genes not in GC and gene length annotation file
  keepV <- intersect(rownames(lengthGcDF), rownames(exprDatDF))
  geneAnno <- lengthGcDF[match(keepV, rownames(lengthGcDF)), ]
  # Set genes with length 0 to length 1 - why? - from Vivek's code
  geneAnno[geneAnno[,1] == 0] <- 1
  exprDatDF <- exprDatDF[match(keepV, rownames(exprDatDF)), ]
  
  # Run CQN with specified depths and no quantile normalization
  cqnDat <- cqn(exprDatDF, lengths = as.numeric(geneAnno[,1])
                , x = as.numeric(geneAnno[,2]), lengthMethod = c("smooth")
                , sqn = FALSE)
  # Get the log2(Normalized FPKM) values
  cqnDat <- cqnDat$y + cqnDat$offset
  cqnDat
}
vzcpCqnDatDF <- RunCQN(vzcpFtExDatDF, ENSEMBLhg19.70UnionAnno)

# Check correlation of expression level to GC content and gene length pre and
# post CQN normalization
CheckCQNnorm <- function (preCQN, postCQN, lengthGcDF) {
  # Remove genes not in GC and gene length annotation file
  keepV <- intersect(rownames(lengthGcDF), rownames(preCQN))
  geneAnno <- lengthGcDF[match(keepV, rownames(lengthGcDF)), ]
  # Set genes with length 0 to length 1 - why? - from Vivek's code
  geneAnno[geneAnno[,1] == 0] <- 1
  keepgenes <- intersect(rownames(preCQN),rownames(postCQN))
  preCQN <- preCQN[match(keepgenes,rownames(preCQN)),]
  postCQN <- postCQN[match(keepgenes,rownames(postCQN)),]
  geneAnno1 <- geneAnno[match(keepgenes,rownames(geneAnno)),]
  
  qcCorrCQNm <- matrix(NA,nrow=ncol(preCQN),ncol=4)
  colnames(qcCorrCQNm) <- c("preNormGCcorr", "preNormLengthCorr"
                            ,"postNormGCcorr", "postNormLengthCorr")
  for (i in 1:nrow(qcCorrCQNm)) {
    qcCorrCQNm[i,1] <- cor(preCQN[,i], geneAnno1[,2], method="spearman")
    qcCorrCQNm[i,2] <- cor(preCQN[,i], geneAnno1[,1], method="spearman")
    qcCorrCQNm[i,3] <- cor(postCQN[,i], geneAnno1[,2], method="spearman")
    qcCorrCQNm[i,4] <- cor(postCQN[,i], geneAnno1[,1], method="spearman")
  }
  qcCorrCQNm
}
qcCorrCQNm <- CheckCQNnorm(vzcpFtExDatDF, vzcpCqnDatDF, ENSEMBLhg19.70UnionAnno)
apply(qcCorrCQNm, 2, quantile)
qcCorrCQNm <- data.frame(qcCorrCQNm)
qcCorrCQNm <- melt(qcCorrCQNm)
colnames(qcCorrCQNm) <- c("CorrType", "Corr")

# Histogram of spearman's rho pre and post normalization for GC and gene length
ggplot(qcCorrCQNm, aes(x = Corr)) +
  facet_wrap(~CorrType, nrow = 2) +
  geom_histogram(binwidth = 0.01) +
  ylab("Counts") +
  xlab("Spearman's rho across samples") +
  ggtitle(paste0(graphsTitle
                 , "\nHistogram: Spearman's rho across samples - pre and post CQN"
                 , "\n"))
ggsave(file = paste0(outGraphs, "CQN_LengthGC_Correlations.pdf"), height = 6)


## Graph data post-CQN

pdf(file = paste0(outGraphs, "QC_CQNlenGC.pdf"))
par(mfrow = c(3, 1))

boxplot(vzcpCqnDatDF, range = 0, col = as.factor(metDatDF$ExpCondition)
        , ylab = "Normalized expression"
        , main = "Boxplot post CQN for length and GC")

plot(density(vzcpCqnDatDF[ ,1]), col = as.factor(metDatDF$ExpCondition)[1]
     , main = "Density post CQN for length and GC")
for (i in 2:dim(vzcpCqnDatDF)[2]) {
  lines(density(vzcpCqnDatDF[ ,i]), col = as.factor(metDatDF$ExpCondition)[i])
}

mdsG = cmdscale(dist(t(vzcpCqnDatDF)), eig = TRUE)
varExpl <- round(mdsG$eig / sum(mdsG$eig), 3)
plot(mdsG$points, col = as.numeric(as.factor(metDatDF$ExpCondition)), pch = 19
     , xlab = paste("Dimension 1 (", 100*varExpl[1], "%)", sep = "")
     , ylab = paste("Dimension 2 (", 100*varExpl[2], "%)", sep = "")
     , main = "MDS plot post CQN for length and GC")
dev.off()


## PCA of sequencing statistics Luis VZ CP RNA-seq

vzcpPicStatsDF <- vzcpPicStatsDF[vzcpPicStatsDF$X %in% gsub("X", "", colnames(vzcpFtExDatDF)), ]
# Mean center (column 12 IGNORED_READS skipped because all 0s)
mcVzcpPicDatDF <- data.frame(SAMPLE = vzcpPicStatsDF$X
                             , apply(vzcpPicStatsDF[ ,-c(1,12)], 2, function(x) scale(x, center = TRUE, scale = FALSE)))

# PCA
datSeqNorm <- t(scale(mcVzcpPicDatDF[ ,-1], scale = FALSE))
pcDatSeq <- prcomp(datSeqNorm);
varexp <- (pcDatSeq$sdev)^2 / sum(pcDatSeq$sdev^2)
print(varexp[1:5])
topPCdatSeq <- pcDatSeq$rotation[ ,1:5]
colnames(topPCdatSeq) <- paste("Seq", colnames(topPCdatSeq))

# Correlation matrix of sequencing stats PCs and sequencing stats

pairsDat <- data.frame(topPCdatSeq, mcVzcpPicDatDF[ ,-1])

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

pdf(paste0(outGraphs, "CovariatesMatrix_SeqStats.pdf"), height = 20, width = 24)
pairs(pairsDat[ ,-11], pch = 19, upper.panel = panel.cor
      , main = "Covariates and MaxQuant Comparison -- |Spearman's rho| correlation values")
dev.off()


## PCA of Luis VZ and CP RNAseq expression

# Centers the mean of all genes - this means the PCA gives us the eigenvectors
# of the geneXgene covariance matrix, allowing us to assess the proportion of
# variance each component contributes to the data
meanCenteredM <- t(scale(t(vzcpCqnDatDF), scale=F))
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

pairsDat <- data.frame(topPCdatSeq
                       , metDatDF[c("RIN.y", "X260.230", "X260.280"
                                    , "ExtractionDate.y", "ExpCondition")])

cond <- labels2colors(metDatDF$ExpCondition)  ## colors

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

pdf(paste0(outGraphs, "CovariatesMatrix_CQNlenGC.pdf"), height = 12, width = 12)
pairs(cbind(topPCs, pairsDat), col = cond, pch = 19
      , upper.panel = panel.cor
      , main = "Covariates and MaxQuant Comparison -- |Spearman's rho| correlation values")
dev.off()
################################################################################

save(vzcpCqnDatDF, topPCdatSeq, file = "../analysis/Expression_CQN_FtG5P4_Luis_RNAseq_VZCP.RData")
# save(vzcpCqnDatDF, file = "../analysis/Expression_CQN_Luis_RNAseq_VZCP.RData")
save(topPCdatSeq, file = "../analysis/SeqPC_Luis_RNAseq_VZCP.RData")

################################################################################

### Regress out covariates

covDF <- data.frame(metDatDF[ ,c("ExpCondition", "RIN.y", "X260.230"
      , "X260.280", "ExtractionDate.y")], topPCdatSeq)
covDF$ExpCondition <- as.numeric(covDF$ExpCondition)
covDF$ExtractionDate.y <- as.numeric(covDF$ExtractionDate.y)

# Regress out confounding variables
RegressCovariates <- function (exM, covDF) {
  X = model.matrix(~ ., data = covDF)
  Y = exM
  beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
  b = as.data.frame(t(beta))
  to_regress = (as.matrix(X[,3:11]) %*% (as.matrix(beta[3:11,])))
  exRegCovM = exM - t(to_regress)
  return(exRegCovM)
}
vzcpNmExDF <- RegressCovariates(vzcpCqnDatDF, covDF)
quantile(vzcpCqnDatDF, c(0,0.025,0.25,0.5,0.75,0.975,1))
quantile(vzcpNmExDF, c(0,0.025,0.25,0.5,0.75,0.975,1))
save(vzcpNmExDF, file = "../analysis/Expression_CQN_RgCv_Luis_RNAseq_VZCP.RData")


## Bootstrap regression
# 
# covDF <- data.frame(metDatDF[ ,c("ExpCondition", "RIN.y", "X260.230"
#   , "X260.280", "ExtractionDate.y")], topPCdatSeq)
# covDF$ExpCondition <- as.numeric(covDF$ExpCondition)
# covDF$ExtractionDate.y <- as.numeric(covDF$ExtractionDate.y)
# 
# # Factors must be coded as numeric for bootstrap
# RegressCovariatesBootstrap <- function (exM, covDF) {
#   # 1) lme using covariates on unnormalized data
#   options(stringsAsFactors=FALSE)
#   boot <- FALSE
#   numboot <- 10
#   bs <- function(formula, data, indices) {
#     d <- data[indices,] # allows bootstrap function to select samples
#     fit <- lm(formula, data=d)
#     return(coef(fit))
#   }
#   exBootRegM <- matrix(NA, nrow = nrow(exM), ncol = ncol(exM))
#   rownames(exBootRegM) <- rownames(exM)
#   colnames(exBootRegM) <- colnames(exM)
#   # ncol(covDF)+1 when condition has 2 levels
#   coefmat <- matrix(NA, nrow = nrow(exM), ncol = ncol(covDF) + 1)
#   set.seed(8675309)
#   for (i in 1:nrow(exM)) {
#     if (i%%1000 == 0) {print(i)}
#     thisexp <- as.numeric(exM[i, ])
#     #       bs.results <- boot(data = data.frame(thisexp, covDF), statistic = bs,
#     #         R = numboot, formula = thisexp~. + age + sex + PMI)
#     bs.results <- boot(data = data.frame(thisexp, covDF), statistic = bs,
#       R = numboot, formula = thisexp~ExpCondition+RIN.y+X260.230+X260.280+ExtractionDate.y+Seq.PC1+Seq.PC2+Seq.PC3+Seq.PC4+Seq.PC5)
#     ## get the median - we can sometimes get NA values here... so let's exclude
#     ## these - old code #bs.stats <- apply(bs.results$t,2,median) 
#     bs.stats <- rep(NA, ncol(bs.results$t)) ##ncol is 3 here (thisexp, construct and extracted)
#     for (n in 1:ncol(bs.results$t)) {
#       bs.stats[n] <- median(na.omit(bs.results$t[ ,n]))
#     }
#     coefmat[i,] <- bs.stats
#     exBootRegM[i,] <- thisexp - bs.stats[3]*covDF[ ,2] - bs.stats[4]*covDF[ ,3] - bs.stats[5]*covDF[ ,4] - bs.stats[6]*covDF[ ,5] - bs.stats[7]*covDF[ ,6] - bs.stats[8]*covDF[ ,7] - bs.stats[9]*covDF[ ,8] - bs.stats[10]*covDF[ ,9] - bs.stats[11]*covDF[ ,10]
#     # cat('Done for Gene',i,'\n')
#   }
# }
# 


