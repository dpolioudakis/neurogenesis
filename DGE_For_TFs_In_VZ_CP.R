# Damon Polioudakis
# 2016-07-21
# Run DGE for TFs identified by Andreas as enriched for targets sites in
# neurogenesis modules from Kang data
################################################################################

rm(list=ls())
sessionInfo()

require(ggplot2)
require(reshape2)
require(xlsx)
require(WGCNA)
require(biomaRt)

### Load data and assign variables

## Load data

# Luis RNAseq fetal brain VZ CP
load("../analysis/Expression_CQN_FtG5P4_Luis_RNAseq_VZCP.Rdata")

# Picard Sequencing Statistics - Luis RNAseq - Seq PC1-5
load("../analysis/SeqPC_Luis_RNAseq_VZCP.RData")

# Luis RNAseq
metDatDF <- read.csv("../metadata/VZCP_sampleinfo.csv", header = TRUE)

# Miller data
load("../orig.data/LCMDE/AllenLCM.Rdata")
MillerExprRAW = AllenLCM$datExpr
MillerMetaRAW = AllenLCM$datTraits
Zones = read.csv("../orig.data/LCMDE/LCM_Zones_CPio.csv")
MillerAnnotRAW = read.csv("../orig.data/LCMDE/annot.csv", row.names = 1)
# Kang
metKangDF <- read.csv("../orig.data/InVivoData/ProcessedKangMetaData.csv")

# Andreas TFs from VJs pipeline
tfDF <- read.xlsx2("../analysis/SummaryVJRegFacEnc.xls", 2)

# # Gene lengths and GC content for Union Exon model
# load("../source/ENSEMBLhg19_UnionAnno.rda")

# Out graphs
outGraphs <- "../analysis/graphs/DGE_For_TFs_In_VZ_CP_"
graphsTitle <- "DGE_For_TFs_In_VZ_CP.R"

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text=element_text(size=18)))
theme_update(plot.title = element_text(size = 16))
################################################################################

### Functions

## Function: DGE Linear model
DGE_Linear_Model <- function (exDatDF, termsDF, mod) {
  
  lmmod <- apply(as.matrix(exDatDF), 1
                 , function(y) {
                   mod <- as.formula(mod)
                   lm(mod, data = termsDF)})
  
  coefmat <- matrix(NA, nrow = nrow(exDatDF)
                    , ncol = length(coef(lmmod[[1]])))
  pvalmat <- matrix(NA, nrow = nrow(exDatDF)
                    , ncol = length(summary(lmmod[[1]])[[4]][ ,4]))
  colnames(coefmat) <- names(coef(lmmod[[1]]))
  rownames(coefmat) <- rownames(exDatDF)
  colnames(pvalmat) <- names(summary(lmmod[[1]])[[4]][ ,4])
  rownames(pvalmat) <- rownames(exDatDF)
  for (i in 1:nrow(exDatDF)) {
    if (i%%100 == 0) {cat(".")}
    coefmat[i, ] <- coef(lmmod[[i]])
    pvalmat[i, ] <- summary(lmmod[[i]])[[4]][ ,4]
  }
  dgeCoefPvalLM <- list(coefmat = coefmat, pvalmat = pvalmat)
  dgeCoefPvalLM
}

## Function: Histogram of DGE expression fold changes 
DGE_Histogram_Of_FoldChange <- function (foldChange, title, outFile) {
  ggDF <- data.frame(foldChange)
  ggplot(data = ggDF, aes(ggDF[ ,1])) +
    geom_histogram(fill = 4, col = "black") +
    xlab("Log2 fold change") +
    ylab("Count") +
    ggtitle(paste0(graphsTitle
                   , "\n", title
                   , "\n"
    ))
  ggsave(paste0(outGraphs, outFile))
}

## Function: Histogram of DGE expression pvals
DGE_Histogram_Of_FoldChange_Pvals <- function (foldChange, title, outFile) {
  ggDF <- data.frame(foldChange)
  ggplot(data = ggDF, aes(ggDF[ ,1])) +
    geom_histogram(fill = 4, col = "black") +
    xlab("P-value") +
    ylab("Count") +
    ggtitle(paste0(graphsTitle
                   , "\n", title
                   , "\n"
    ))
  ggsave(paste0(outGraphs, outFile))
}

## Function: Convert Ensembl IDs to Gene Symbols
ConvertEnsemblTranscriptToGeneSym <- function (ensemblList) {
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  moduleGenes <- data.frame(ensemblList)
  # bioMart manual:
  #http://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf
  # Attribute: values to retrieve
  # Filters: input query type
  # Values: input query
  #ensembl <- useMart("ensembl")
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
  ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)
  # Data frame of module Ensembl IDs and gene symbols
  ensemblGeneSymDF <- getBM(  attributes = c("ensembl_gene_id", "hgnc_symbol")
                              , filters = "hgnc_symbol"
                              , values = moduleGenes
                              , mart = ensembl
  )
  ensemblGeneSymDF
}

################################################################################

### DGE Luis RNAseq

# (CQN outputs log2 transformed values, these are input into linear model, and
#  then beta values output from linear mode are the log2 fold change)

## DGE Luis Linear model
termsDF <- data.frame(metDatDF[ ,c("ExpCondition", "RIN.y", "X260.230"
                                   , "X260.280", "ExtractionDate.y")], topPCdatSeq)

dgeLuVzcpLM <- DGE_Linear_Model(vzcpCqnDatDF, termsDF
  , "y~ExpCondition+RIN.y+X260.230+X260.280+ExtractionDate.y+Seq.PC1+Seq.PC2+Seq.PC3+Seq.PC4+Seq.PC5")

# Plot histogram of Luis VZ CP fold changes
DGE_Histogram_Of_FoldChange(dgeLuVzcpLM$coefmat[ ,2]
                            , "Fold change histogram: Luis RNAseq VZ CP"
                            , "LuisVzCp_FoldChange_Hist_CQNlenGC_Lm.pdf")
# Plot histogram of Luis VZ CP fold changes p-values
DGE_Histogram_Of_FoldChange_Pvals(dgeLuVzcpLM$pvalmat[ ,2]
                                  , "Fold change p-values histogram: Miller VZ CP"
                                  , "LuisVzCp_FoldChangePval_Hist_CQNlenGC_Lm.pdf")

# Save Luis DGE from linear model
save(dgeLuVzcpLM, file = "../analysis/DGE_Luis_VzCp_FtG5P4_CQNlenGc_Lm.rdata")
################################################################################

### Subset Luis VZ CP RNAseq DGE to Andreas TFs

## Convert Ensembl IDs to Gene Symbols

tfsEnsemblGeneDF <- ConvertEnsemblTranscriptToGeneSym(tfDF$GENENAME)

## Subset DGE by TFs

table(tfsEnsemblGeneDF$ensembl_gene_id %in% gsub("\\..*", "", rownames(dgeLuVzcpLM$coefmat)))
for(mName in names(dgeLuVzcpLM)) {
  m <- dgeLuVzcpLM[[mName]]
  dgeLuVzcpLM[[mName]] <- m[match(tfsEnsemblGeneDF$ensembl_gene_id, gsub("\\..*", "", rownames(m))), 2]
}

# ## FDR
# corrected <- fdrtool(vzcpTFs, statistic= "pvalue", plot = TRUE)

## Histogram of fold changes

DGE_Histogram_Of_FoldChange(dgeLuVzcpLM$coefmat
  , paste0(
    "Fold change histogram of Luis VZ CP RNAseq DGE subset to Andreas TFs"
    , "\nCQN length GC normalized"
    , "\nLinear module terms: Brain Region, RIN, 260/230, 260/280, "
    , "\nExtraction Date, Seq PC1-5")
  , "LuisVzCp_FoldChange_Hist_TFs_CQNlenGC_Lm.pdf")

## Histogram of fold change p-pvalues

DGE_Histogram_Of_FoldChange_Pvals(dgeLuVzcpLM$pvalmat
  , paste0("Fold change p-values histogram of Luis VZ CP RNAseq DGE subset to Andreas TFs"
    , "\nCQN length GC normalized"
    , "\nLinear module terms: Brain Region, RIN, 260/230, 260/280, "
    , "\nExtraction Date, Seq PC1-5")
  , "LuisVzCp_FoldChangePval_Hist_TFs_CQNlenGC_Lm.pdf")


## Make TF and log2 fold change table combining log2 fold changes, p-values
## for TFs and Andreas TF table

dgeDF <- data.frame(ENSEMBL = names(dgeLuVzcpLM$coefmat)
                    , LUIS_VZCP_LOG2_FC = dgeLuVzcpLM$coefmat
                    , LUIS_VZCP_PVAL = dgeLuVzcpLM$pvalmat)
dgeDF$ENSEMBL <- gsub("\\..*", "", dgeDF$ENSEMBL)
dgeDF <- merge(dgeDF, tfsEnsemblGeneDF, by.x = "ENSEMBL"
                     , by.y = "ensembl_gene_id", all.y = TRUE)
vzcpTFsDF <- merge(tfDF, dgeDF, by.x = "GENENAME", by.y = "hgnc_symbol", all.x = TRUE)

################################################################################

### DGE Miller

### Format and subset of data

ZonesMiller = MillerMetaRAW$structure_acronym

# Change name of Kang metadata sample column from X to SampleID
names(metKangDF)[1] = "SampleID"

## Select cortical layers only

MillerMetaRAW$Labels = rep(NA, nrow(MillerMetaRAW))
for (Label in Zones$Label)
  MillerMetaRAW[grep(Label, ZonesMiller), "Labels"] = Label

Zoneindex = which(! is.na(MillerMetaRAW$Labels))

MillerExpr = MillerExprRAW[, Zoneindex]
MillerMeta = MillerMetaRAW[Zoneindex, ]

## Only ENTREZ annotated probes

MillerExpr = MillerExpr[! is.na(MillerAnnotRAW$ENTREZ_ID), ]
MillerAnnot = MillerAnnotRAW[! is.na(MillerAnnotRAW$ENTREZ_ID), ]

## Get maximum expression probe

# Miller Dataset

Genes = unique(MillerAnnot$ENTREZ_ID)

keepind = matrix(nrow = 0, ncol = 0);

for (ii in 1:length(Genes)) {
  genematchind = which(MillerAnnot$ENTREZ_ID == Genes[ii]);
  if (length(genematchind) > 1) {
    themeans = rowMeans(MillerExpr[genematchind, ]);
    maxind = which(themeans == max(themeans))[1];
    keepind = c(keepind, genematchind[maxind]);
  } else {
    keepind = c(keepind, genematchind);
  }
}


MillerExpr = MillerExpr[keepind, ];
rownames(MillerExpr) = Genes;
MillerAnnot = MillerAnnot[keepind, ];

# Cleanup Miller zone labels
millerZones <- MillerMeta[match(MillerMeta$well_id, colnames(MillerExpr)), ]$Labels
millerZones <- gsub(".*SG.*", "SG", millerZones)
millerZones <- gsub(".*MZ.*", "MZ", millerZones)
millerZones <- gsub(".*CP.*", "CP", millerZones)
millerZones <- gsub(".*SP.*", "SP", millerZones)
millerZones <- gsub(".*IZ.*", "IZ", millerZones)
millerZones <- gsub(".*SZ.*", "SZ", millerZones)
millerZones <- gsub(".*VZ.*", "VZ", millerZones)
################################################################################

### DGE Miller VZ CP

## DGE Linear Model Miller VZ CP
# Subset expression matrix to CP and VZ
exDF <- MillerExpr[ ,millerZones == c("CP", "VZ")]
zones <- millerZones[millerZones == c("CP", "VZ")]
# Terms for LM
termsDF <- data.frame(Zone = zones)
# DGE Linear Model
dgeMlVzcpLM <- DGE_Linear_Model(exDF, termsDF, "y~Zone")
# Plot histogram of Miller VZ CP fold changes
DGE_Histogram_Of_FoldChange(dgeMlVzcpLM$coefmat[ ,2]
                            , "Fold change histogram: Miller VZ CP"
                            , "MillerVzCp_FoldChange_Hist_Lm.pdf")
# Plot histogram of Miller VZ CP fold changes p-values
DGE_Histogram_Of_FoldChange_Pvals(dgeMlVzcpLM$pvalmat[ ,2]
                                  , "Fold change p-values histogram: Miller VZ CP"
                                  , "MillerVzCp_FoldChangePval_Hist_Lm.pdf")

## DGE Linear Model Miller SVZ CP
# Subset expression matrix to SVZ CP
exDF <- MillerExpr[ ,millerZones == c("CP", "SZ")]
zones <- millerZones[millerZones == c("CP", "SZ")]
# Terms for LM
termsDF <- data.frame(Zone = zones)
# DGE Linear Model
dgeMlSvzcpLM <- DGE_Linear_Model(exDF, termsDF, "y~Zone")
# Plot histogram of Miller SVZ CP fold changes
DGE_Histogram_Of_FoldChange(dgeMlVzcpLM$coefmat[ ,2]
                            , "Fold change histogram: Miller SVZ CP"
                            , "MillerSvzCp_FoldChange_Hist_Lm.pdf")
# Plot histogram of Miller SVZ CP fold changes p-values
DGE_Histogram_Of_FoldChange_Pvals(dgeMlVzcpLM$pvalmat[ ,2]
                                  , "Fold change p-values histogram: Miller SVZ CP"
                                  , "MillerSvzCp_FoldChangePval_Hist_Lm.pdf")

## DGE Linear Model Miller SVZ SP
# Subset expression matrix to CP and VZ
exDF <- MillerExpr[ ,millerZones == c("SP", "SZ")]
zones <- millerZones[millerZones == c("SP", "SZ")]
# Terms for LM
termsDF <- data.frame(Zone = zones)
# DGE Linear Model
dgeMlSvzspLM <- DGE_Linear_Model(exDF, termsDF, "y~Zone")
# Plot histogram of Miller VZ CP fold changes
DGE_Histogram_Of_FoldChange(dgeMlSvzspLM$coefmat[ ,2]
                            , "Fold change histogram: Miller SVZ SP"
                            , "MillerSvzSp_FoldChange_Hist_Lm.pdf")
# Plot histogram of Miller VZ CP fold changes p-values
DGE_Histogram_Of_FoldChange_Pvals(dgeMlSvzspLM$pvalmat[ ,2]
                                  , "Fold change p-values histogram: Miller SVZ SP"
                                  , "MillerSvzSp_FoldChangePval_Hist_Lm.pdf")

# Save Miller DGE from linear model
save(dgeMlVzcpLM, dgeMlSvzcpLM, dgeMlSvzspLM,
     file = "../analysis/DGE_Miller_VzCpSvzCpSvzSp_LM.rdata")
################################################################################

### Subset Miller VZ SVZ CP RNAseq DGE to Andreas TFs

## Convert Ensembl IDs to Gene Symbols
tfsEnsemblGeneDF <- ConvertEnsemblTranscriptToGeneSym(tfDF$GENENAME)

## Subset DGE by TFs
# VZ CP
for(mName in names(dgeMlVzcpLM)) {
  m <- dgeMlVzcpLM[[mName]]
  row.names(m) <- MillerAnnot$SYMBOL[match(row.names(m), MillerAnnot$ENTREZ_ID)]
  dgeMlVzcpLM[[mName]] <- m[match(tfsEnsemblGeneDF$hgnc_symbol, rownames(m)), 2]
}
# SVZ CP
for(mName in names(dgeMlSvzcpLM)) {
  m <- dgeMlSvzcpLM[[mName]]
  row.names(m) <- MillerAnnot$SYMBOL[match(row.names(m), MillerAnnot$ENTREZ_ID)]
  dgeMlSvzcpLM[[mName]] <- m[match(tfsEnsemblGeneDF$hgnc_symbol, rownames(m)), 2]
}
# SVZ SP
for(mName in names(dgeMlSvzspLM)) {
  m <- dgeMlSvzspLM[[mName]]
  row.names(m) <- MillerAnnot$SYMBOL[match(row.names(m), MillerAnnot$ENTREZ_ID)]
  dgeMlSvzspLM[[mName]] <- m[match(tfsEnsemblGeneDF$hgnc_symbol, rownames(m)), 2]
}

# ## FDR
# corrected <- fdrtool(vzcpTFs, statistic= "pvalue", plot = TRUE)

## Histogram of VZ CP fold changes
DGE_Histogram_Of_FoldChange(dgeMlVzcpLM$coefmat
  , "Fold change histogram of Miller VZ CP DGE subset to Andreas TFs"
  , "MillerVzCp_FoldChange_Hist_TFs_Lm.pdf")
# And p-values
DGE_Histogram_Of_FoldChange_Pvals(dgeMlVzcpLM$pvalmat
 , "Fold change p-values histogram of Miller VZ CP DGE subset to Andreas TFs"
 , "MillerVzCp_FoldChangePval_Hist_TFs_Lm.pdf")

## Histogram of SVZ CP fold changes
DGE_Histogram_Of_FoldChange(dgeMlSvzcpLM$coefmat
  , "Fold change histogram of Miller SVZ CP DGE subset to Andreas TFs"
  , "MillerSvzCp_FoldChange_Hist_TFs_Lm.pdf")
# And p-values
DGE_Histogram_Of_FoldChange_Pvals(dgeMlSvzcpLM$pvalmat
  , "Fold change p-values histogram of Miller SVZ CP DGE subset to Andreas TFs"
  , "MillerSvzCp_FoldChangePval_Hist_TFs_Lm.pdf")

## Histogram of SVZ SP fold changes
DGE_Histogram_Of_FoldChange(dgeMlSvzspLM$coefmat
  , "Fold change histogram of Miller SVZ SP DGE subset to Andreas TFs"
  , "MillerSvzSp_FoldChange_Hist_TFs_Lm.pdf")
# And p-values
DGE_Histogram_Of_FoldChange_Pvals(dgeMlSvzspLM$pvalmat
  , "Fold change p-values histogram of Miller SVZ SP DGE subset to Andreas TFs"
  , "MillerSvzSp_FoldChangePval_Hist_TFs_Lm.pdf")


## Add to TF and log2 fold change table combining log2 fold changes, p-values
## for TFs and Andreas TF table

# Not filtering for significance
# VZ CP
dgeDF <- data.frame(GENE_SYM = names(dgeMlVzcpLM$coefmat)
                    , MILLER_VZCP_LOG2_FC = dgeMlVzcpLM$coefmat
                    , MILLER_VZCP_PVAL = dgeMlVzcpLM$pvalmat)
vzcpTFsDF <- merge(vzcpTFsDF, dgeDF, by.x = "GENENAME", by.y = "GENE_SYM", all.x = TRUE)
# SVZ CP
dgeDF <- data.frame(GENE_SYM = names(dgeMlSvzcpLM$coefmat)
                    , MILLER_SVZCP_LOG2_FC = dgeMlSvzcpLM$coefmat
                    , MILLER_SVZCP_PVAL = dgeMlSvzcpLM$pvalmat)
vzcpTFsDF <- merge(vzcpTFsDF, dgeDF, by.x = "GENENAME", by.y = "GENE_SYM", all.x = TRUE)
# SVZ SP
dgeDF <- data.frame(GENE_SYM = names(dgeMlSvzspLM$coefmat)
                    , MILLER_SVZSP_LOG2_FC = dgeMlSvzspLM$coefmat
                    , MILLER_SVZSP_PVAL = dgeMlSvzspLM$pvalmat)
vzcpTFsDF <- merge(vzcpTFsDF, dgeDF, by.x = "GENENAME", by.y = "GENE_SYM", all.x = TRUE)
################################################################################

### Heatmaps of expression for Luis VZ CP, and Miller VZ SVZ CP


# Select columns of interest and remove duplicate rows
df <- unique(vzcpTFsDF[ ,c("GENENAME"
                           , "LUIS_VZCP_LOG2_FC", "LUIS_VZCP_PVAL"
                           , "MILLER_VZCP_LOG2_FC", "MILLER_VZCP_PVAL"
                           , "MILLER_SVZCP_LOG2_FC", "MILLER_SVZCP_PVAL"
                           , "MILLER_SVZSP_LOG2_FC", "MILLER_SVZSP_PVAL"
                           , "module")])
df[order(df$module), ]

# Reverse fold change sign so that positive is higher expression in CP
df$LUIS_VZCP_LOG2_FC <- -df$LUIS_VZCP_LOG2_FC
df$MILLER_VZCP_LOG2_FC <- -df$MILLER_VZCP_LOG2_FC
df$MILLER_SVZCP_LOG2_FC <- -df$MILLER_SVZCP_LOG2_FC
df$MILLER_SVZSP_LOG2_FC <- -df$MILLER_SVZSP_LOG2_FC

## Heatmaps of expression for Luis VZ CP, and Miller VZ SVZ CP
# Split by Kang module

modL <- split(df, df$module)

pdf(paste0(outGraphs, "Heatmaps_DGE.pdf"))
for (modDF in modL) {
  print(head(modDF))
  
  # Matrix of log2 fold changes and p-values for labeling heatmap
  textMatrix <- paste0(round(as.matrix(modDF[ ,c(
    "LUIS_VZCP_LOG2_FC", "MILLER_VZCP_LOG2_FC", "MILLER_SVZCP_LOG2_FC", "MILLER_SVZSP_LOG2_FC")]), 2)
    , "\n"
    , round(as.matrix(modDF[ ,c(
      "LUIS_VZCP_PVAL", "MILLER_VZCP_PVAL", "MILLER_SVZCP_PVAL", "MILLER_SVZSP_PVAL")]), 4))
  
  # Set margins
  par(mar = c(6,10,5,1))
  # Plot
  labeledHeatmap(Matrix = as.matrix(modDF[ ,c("LUIS_VZCP_LOG2_FC"
                                              , "MILLER_VZCP_LOG2_FC"
                                              , "MILLER_SVZCP_LOG2_FC"
                                              , "MILLER_SVZSP_LOG2_FC")])
            , xLabels = c("Luis VZ / CP", "Miller VZ / CP", "Miller SVZ / CP", "Miller SVZ / SP")
            , yLabels = paste0("ME", modDF$module)
            , ySymbols = modDF$GENENAME
            , colorLabels = FALSE
            , colors = blueWhiteRed(50)
            , textMatrix = textMatrix
            , setStdMargins = FALSE
            , cex.text = 0.5
            , main = paste0(graphsTitle
              , "\nNormalized expression log2 fold change and p-value"
              , "\n", modDF$module[1], " module derived from Kang enriched TFs")
  )
}
dev.off()

# # Code for bargraphs instead of heatmaps
# pdf(paste0(outGraphs, "TEST.pdf"))
# modL <- split(df, df$module)
# for (modDF in modL) {
#   print(modDF)
#   ggDF <- melt(modDF)
#   print(head(ggDF))
#   ggDF[is.na(ggDF)] <- 0
#   head(ggDF)
#   plot <- ggplot(ggDF, aes(y = value, x = GENENAME)) +
#     geom_bar(aes(fill = variable), position = "dodge", stat = "identity")
#   print(plot)
# }
# dev.off()



################################################################################
