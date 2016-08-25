# Damon Polioudakis
# 2014-04-13
# Plot modules defined in Kang dataset across Kang, Miller, Luis VZ/CP, hNPCs,
# and Primary cultures
################################################################################

rm(list = ls())

require(WGCNA)
require(ggplot2)
require(reshape2)
require(biomaRt)
require(gridExtra)

options(stringsAsFactors = F)
args = commandArgs(trailingOnly = T)

## Load data and assign variables

# Kang
load("../orig.data/InVivoData/WGCNAinput_SestanBrain.RData")
KangExpr = t(datExpr)
rm(datExpr)
# sampleKey is ProcessedKangMetaData.csv
metKangDF <- sampleKey
# metKangDF <- read.csv("../orig.data/InVivoData/ProcessedKangMetaData.csv")
KangAnnotRAW = read.csv("../orig.data/InVivoData/annot.csv", row.names = 1)
KangAnnotnet2 = KangAnnotRAW[which(rownames(KangAnnotRAW) %in% rownames(KangExpr)), ] 
load("../orig.data/InVivoData/net2_power16_cutHeigt0.15.RData")

# Miller
load("../orig.data/LCMDE/AllenLCM.Rdata")
MillerExprRAW = AllenLCM$datExpr
MillerMetaRAW = AllenLCM$datTraits
Zones = read.csv("../orig.data/LCMDE/LCM_Zones_CPio.csv")
MillerAnnotRAW = read.csv("../orig.data/LCMDE/annot.csv", row.names = 1)

# Luis VZ/CP RNAseq
# Raw HTSeq counts
# vzcpExDF <- read.csv("../data/bulk_VZ_CP_from_ATAC/Exprs_HTSCexon.csv"
#                      , row.names = 1)
# CQN for GC and length + Regressed out covariates
load("../analysis/Expression_CQN_RgCv_Luis_RNAseq_VZCP.RData")
vzcpExDF <- vzcpNmExDF
rm(vzcpNmExDF)
# CQN for GC and length
# load("../analysis/Expression_CQN_Luis_RNAseq_VZCP.RData")
# vzcpExDF <- vzcpCqnDatDF
# rm(vzcpCqnDatDF)
metDatDF <- read.csv("../metadata/VZCP_sampleinfo.csv", header = TRUE)

# hNPs
HNPAnnot = read.csv("../orig.data/HNPData1.4.8.12/annot.csv", row.names = 1)
HNPExpr = read.csv("../orig.data/HNPData1.4.8.12/exprdata.csv", row.names = 1)
HNPMeta = read.csv("../orig.data/HNPData1.4.8.12/sampleinfo.csv", row.names = 1)

# Primary cultures from human fetal VZ/CP
prmExDF <- read.csv("../data/primaryCultures_exprdata.csv")
prmAnnotDF <- read.csv("../metadata/primaryCultures_annot.csv")
prmMetDF <- read.csv("../metadata/primaryCultures_metadata.csv")

## Variables
graphCodeTitle <- "Kang_ModuleEG_Expression.R"
outGraphPfx <- "../analysis/graphs/Kang_ModuleEG_Expression_"

## Output Directory
dir.create("../analysis/graphs", recursive = TRUE)

## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text=element_text(size=18)))
theme_update(plot.title = element_text(size = 16))
################################################################################

### Functions

## Function: Convert Ensembl IDs to Gene Symbols
ConvertEnsemblToEntrez <- function (ensemblList) {
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
  ensemblGeneSymDF <- getBM(  attributes = c("ensembl_gene_id", "entrezgene")
                              , filters = "ensembl_gene_id"
                              , values = moduleGenes
                              , mart = ensembl
  )
  ensemblGeneSymDF
}


################################################################################

### Format and subset of data

## Luis RNAseq VZ/CP Ensembl IDs to Entrez
rownames(vzcpExDF) <- gsub("\\..*", "", rownames(vzcpExDF))
ensemblEntrezDF <- ConvertEnsemblToEntrez(row.names(vzcpExDF))
# Remove ensembl IDs not in bioMart and add Entrez
vzcpExDF <- merge(ensemblEntrezDF, vzcpExDF
                  , by.x = "ensembl_gene_id", by.y = "row.names")
dim(vzcpExDF) # 53278
# Remove genes with no Entrez ID
vzcpExDF <- vzcpExDF[! is.na(vzcpExDF$entrezgene), ]
dim(vzcpExDF) # 24119


## Primary Cultures VZ/CP ILMN IDs to Entrez
prmExDF <- merge(prmAnnotDF, prmExDF, by.x = "row.names", by.y = "row.names")


## Format and subset Miller

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

HNPExpr = HNPExpr[! is.na(HNPAnnot$ENTREZ_GENE_ID), ]
HNPAnnot = HNPAnnot[! is.na(HNPAnnot$ENTREZ_GENE_ID), ]


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


# HNP Dataset

Genes = unique(HNPAnnot$ENTREZ_GENE_ID)

keepind = matrix(nrow = 0, ncol = 0);
for (ii in 1:length(Genes)) {
  genematchind = which(HNPAnnot$ENTREZ_GENE_ID == Genes[ii]);
  if (length(genematchind) > 1) {
    themeans = rowMeans(HNPExpr[genematchind, ]);
    maxind = which(themeans == max(themeans))[1];
    keepind = c(keepind, genematchind[maxind]);
  } else {
    keepind = c(keepind, genematchind);
  }
}

HNPExpr = HNPExpr[keepind, ];
rownames(HNPExpr) = Genes;
HNPAnnot = HNPAnnot[keepind, ];


## Select only ENTREZ genes evaluated in the Kang datasets

KangAnnot = KangAnnotnet2

Index = MillerAnnot$ENTREZ_ID %in% KangAnnot$ENTREZ_ID

MillerAnnot = MillerAnnot[Index, ]
MillerExpr = MillerExpr[Index, ]

nrow(MillerExpr)
# 16534

Index = HNPAnnot$ENTREZ_GENE_ID %in% KangAnnot$ENTREZ_ID
HNPAnnot = HNPAnnot[Index, ]
HNPExpr = HNPExpr[Index, ]

nrow(HNPExpr)
# 16679

## Limit to ENTREZ IDs present in all datasets
## Arbitrary for Module preservation analysis, but mandatory for follow up
## procedures

if (TRUE) {
  ENTREZset = intersect(
    intersect(
      intersect(
        intersect(KangAnnotnet2$ENTREZ_ID, MillerAnnot$ENTREZ_ID)
          , HNPAnnot$ENTREZ_GENE_ID)
        , vzcpExDF$entrezgene)
    , prmExDF$ENTREZ_GENE_ID)
  
  Index = match(ENTREZset, KangAnnotnet2$ENTREZ_ID)
  KangExpr = KangExpr[Index, ]
  KangAnnot = KangAnnotnet2[Index, ]
  
  Index = match(ENTREZset, MillerAnnot$ENTREZ_ID)
  MillerExpr = MillerExpr[Index, ]
  MillerAnnot = MillerAnnot[Index, ]
  
  Index = match(ENTREZset, HNPAnnot$ENTREZ_GENE_ID)
  HNPExpr = HNPExpr[Index, ]
  HNPAnnot = HNPAnnot[Index, ]
  
  Index = match(ENTREZset, vzcpExDF$entrezgene)
  vzcpExDF = vzcpExDF[Index, ]
  
  Index = match(ENTREZset, prmExDF$ENTREZ_GENE_ID)
  prmExDF = prmExDF[Index, ]
}

if (nrow(HNPExpr) == nrow(KangExpr) & nrow(HNPExpr) == nrow(MillerExpr) & nrow(MillerExpr) == nrow(vzcpExDF))
  cat('there are', nrow(HNPExpr), 'genes in all datasets\n')
# there are 16026 genes in all datasets
################################################################################

### Plot ME expression from Kang modules in datasets

## Calculate MEs
netColors <- net2$colors[match(ENTREZset, KangAnnotnet2$ENTREZ_ID)]
kangMEs <- moduleEigengenes(t(KangExpr), netColors)$eigengenes
hNPMEs <- moduleEigengenes(t(HNPExpr), netColors)$eigengenes
millerMEs <- moduleEigengenes(t(MillerExpr), netColors)$eigengenes
luisVzcpMEs <- moduleEigengenes(t(vzcpExDF[ ,-c(1:2)]), netColors)$eigengenes
prmMEs <- moduleEigengenes(t(prmExDF[ ,-c(1:3)]), netColors)$eigengenes

# Subset Kang metadata to Kang samples present in expression DF
sbMetKangDF <- metKangDF[match(colnames(KangExpr), metKangDF$SampleID), ]
sbMetKangDF <- sbMetKangDF[! is.na(sbMetKangDF[1]), ]

## Extract stage for each sample, ordered by expression DF - same order as
## module eigengene list
# Kang
kangStage <- sbMetKangDF[match(sbMetKangDF$SampleID, colnames(KangExpr)), ]$Stage
# hNPs
hNPwK <- HNPMeta[match(HNPMeta$RNAID, gsub("X", "", colnames(HNPExpr))), ]$DiffWk
# Miller
millerZones <- MillerMeta[match(MillerMeta$well_id, colnames(MillerExpr)), ]$Labels
millerZones <- gsub(".*SG.*", "SG", millerZones)
millerZones <- gsub(".*MZ.*", "MZ", millerZones)
millerZones <- gsub(".*CP.*", "CP", millerZones)
millerZones <- gsub(".*SP.*", "SP", millerZones)
millerZones <- gsub(".*IZ.*", "IZ", millerZones)
millerZones <- gsub(".*SZ.*", "SZ", millerZones)
millerZones <- gsub(".*VZ.*", "VZ", millerZones)
millerZones <- as.factor(millerZones)
# Primary cultures
df <- prmMetDF[match(row.names(prmMetDF), colnames(prmExDF[ ,-(1:3)])), ]
prmWkZoneDF <- data.frame(WEEK = gsub("^T([0-9]).*", "\\1", row.names(df), perl = TRUE))
prmWkZoneDF$ZONE <- gsub("^T[0-9]_(.*)_.*", "\\1", row.names(df), perl = TRUE)

# Luis RNAseq VZ/CP
df <- data.frame(Zone = metDatDF$ExpCondition, luisVzcpMEs)
df$Zone <- factor(df$Zone, levels = c("VZ", "CP"))
luisGgDF <- melt(df, id.vars = "Zone"
  , variable.name = "ME", value.name = "Expression")
luisGgLDF <- split(luisGgDF, luisGgDF$ME)
luisGGL <- lapply(names(luisGgLDF), function(name) {
  mnGgDF <- aggregate(luisGgLDF[[name]], list(luisGgLDF[[name]]$Zone), mean)
  mnGgDF$Zone <- mnGgDF$Group.1
  ggplot(luisGgLDF[[name]], aes(x = Zone, y = Expression, group = 1)) +
    geom_jitter(size = 0.5, alpha = 0.25, width = 0.2, height = 0.2) +
    geom_point(data = mnGgDF, aes(x = Zone, y = Expression, group = 1), color = "blue") +
    geom_line(data = mnGgDF, aes(x = Zone, y = Expression, group = 1), color = "blue") + 
    xlab("Zone\n") +
    ylab("") +
    ggtitle(paste0(name, "\n"))
})

# Kang
df <- data.frame(Stage = kangStage, kangMEs)
kangGgDF <- melt(df)
kangGgDF <- melt(df, id.vars = "Stage"
  , variable.name = "ME", value.name = "Expression")
kangGgLDF <- split(kangGgDF, kangGgDF$ME)
kangGGL <- lapply(names(kangGgLDF), function(name) {
  mnGgDF <- aggregate(kangGgLDF[[name]], list(kangGgLDF[[name]]$Stage), mean)
  ggplot(kangGgLDF[[name]], aes(x = Stage, y = Expression, group = 1)) +
    geom_jitter(size = 0.5, alpha = 0.25, width = 0.5, height = 0.5) +
    geom_point(data = mnGgDF, aes(x = Stage, y = Expression, group = 1), color = "red") +
    geom_line(data = mnGgDF, aes(x = Stage, y = Expression, group = 1), color = "red") + 
    xlab("Stage\n") +
    ylab("ME expression") +
    ggtitle(paste0(name, "\n"))
})

# Miller
df <- data.frame(Zone = millerZones, millerMEs)
df$Zone <- factor(df$Zone, levels = c("VZ", "SZ", "IZ", "SP", "CP", "MZ", "SG"))
millerGgDF <- melt(df)
millerGgDF <- melt(df, id.vars = "Zone"
  , variable.name = "ME", value.name = "Expression")
millerGgLDF <- split(millerGgDF, millerGgDF$ME)
millerGGL <- lapply(names(millerGgLDF), function(name) {
  mnGgDF <- aggregate(millerGgLDF[[name]], list(millerGgLDF[[name]]$Zone), mean)
  mnGgDF$Zone <- mnGgDF$Group.1
  ggplot(millerGgLDF[[name]], aes(x = Zone, y = Expression, group = 1)) +
    geom_jitter(size = 0.5, alpha = 0.25, width = 0.5, height = 0.5) +
    geom_point(data = mnGgDF, aes(x = Zone, y = Expression, group = 1), color = "purple") +
    geom_line(data = mnGgDF, aes(x = Zone, y = Expression, group = 1), color = "purple") + 
    xlab("Zone\n") +
    ylab("") +
    ggtitle(paste0(name, "\n"))
})

# hNPs
df <- data.frame(Week = hNPwK, hNPMEs)
hnpGgDF <- melt(df)
hnpGgDF <- melt(df, id.vars = "Week"
  , variable.name = "ME", value.name = "Expression")
hnpGgLDF <- split(hnpGgDF, hnpGgDF$ME)
hnpGGL <- lapply(names(hnpGgLDF), function(name) {
  mnGgDF <- aggregate(hnpGgLDF[[name]], list(hnpGgLDF[[name]]$Week), mean)
  mnGgDF$Week <- mnGgDF$Group.1
  ggplot(hnpGgLDF[[name]], aes(x = Week, y = Expression, group = 1)) +
    geom_jitter(size = 0.5, alpha = 0.25, width = 0.2, height = 0.2) +
    geom_point(data = mnGgDF, aes(x = Week, y = Expression, group = 1), color = "green") +
    geom_line(data = mnGgDF, aes(x = Week, y = Expression, group = 1), color = "green") + 
    xlab("Week\n") +
    ylab("") +
    ggtitle(paste0(name, "\n"))
})

# Primary cultures
df <- data.frame(WEEK = prmWkZoneDF$WEEK, ZONE = prmWkZoneDF$ZONE, prmMEs)
ggDF <- melt(df, id.vars = c("WEEK", "ZONE")
  , variable.name = "ME", value.name = "EXPRESSION")
ggLDF <- split(ggDF, ggDF$ME)
prmGGL <- lapply(names(ggLDF), function(name) {
  mnGgDF <- aggregate(EXPRESSION ~ WEEK + ZONE, ggLDF[[name]], mean)
  ggplot(ggLDF[[name]], aes(x = WEEK, y = EXPRESSION, color = ZONE)) +
    geom_jitter(size = 0.5, alpha = 0.5, width = 0.2, height = 0.2) +
    geom_point(data = mnGgDF, aes(x = WEEK, y = EXPRESSION, group = ZONE)) +
    geom_line(data = mnGgDF, aes(x = WEEK, y = EXPRESSION, linetype = ZONE, group = ZONE)) + 
    theme(legend.title = element_blank()) +
    xlab("Week\n") +
    ylab("") +
    ggtitle(paste0(name, "\n"))
})

# Combine to set layout of grid.arrange
ggLL <- mapply(list, kangGGL, luisGGL, millerGGL, hnpGGL, prmGGL)

pdf(paste0(outGraphPfx, "Split_Line_Graphs_LuisNm.pdf"), height = 120, width = 16)
do.call("grid.arrange", c(ggLL, ncol = 5))
dev.off()


# # Calculate mean ME expression at each time point or brain region
# Dev_Time_Point_Mean_MEexpr <- function (timePoint, MEs) {
#   df <- cbind(data.frame(Time_Point = timePoint), MEs)
#   tmpL <- split(df, df$Time_Point)
#   # Remove Time_Point column
#   tmpL <- lapply(tmpL, function(df) df[ ,-1])
#   outM <- sapply(tmpL, function(timePoint) {
#     apply(timePoint, 2, mean)
#   })
#   outM
# }
# 
# # Kang
# ggM <- Dev_Time_Point_Mean_MEexpr(kangStage, kangMEs)
# ggDF <- data.frame(t(ggM))
# ggDF$Time_Point <- row.names(ggDF)
# kangGgDF <- melt(ggDF, id.vars = "Time_Point"
#                  , variable.name = "ME", value.name = "Expression")
# kangGgDF$Data_Set <- "Kang"
# 
# # hNPs
# ggM <- Dev_Time_Point_Mean_MEexpr(hNPwK, hNPMEs)
# ggDF <- data.frame(t(ggM))
# ggDF$Time_Point <- row.names(ggDF)
# gghNPdF <- melt(ggDF, id.vars = "Time_Point"
#                 , variable.name = "ME", value.name = "Expression")
# gghNPdF$Data_Set <- "hNP"
# 
# # Miller
# ggM <- Dev_Time_Point_Mean_MEexpr(millerZones, millerMEs)
# ggDF <- data.frame(t(ggM))
# ggDF$Time_Point <- row.names(ggDF)
# # Convert zones to numeric for plotting on Kang time course
# ggDF$Time_Point <- as.numeric(factor(ggDF$Time_Point
#                                      , levels = c("VZ", "SZ", "IZ", "SP", "CP", "MZ", "SG")))
# ggMillerdF <- melt(ggDF, id.vars = "Time_Point"
#                    , variable.name = "ME", value.name = "Expression")
# ggMillerdF$Data_Set <- "Miller"
# 
# # Combine for ggplot2 plotting
# ggDF <- rbind(kangGgDF, gghNPdF, ggMillerdF)
# ggDF$Time_Point <- as.numeric(ggDF$Time_Point)
# 
# ggplot(ggDF, aes(x = Time_Point, y = Expression, color = Data_Set)) +
#   facet_wrap(~ME) +
#   geom_line() +
#   xlab("Stage / Week / Zone") +
#   ylab("ME Expression") +
#   theme(axis.text.x = element_blank()) +
#   ggtitle(paste0(graphCodeTitle
#                  ,"\nKang Module Eigengene Expression"
#                  ,"\nX-axis:"
#                  ,"\nKang: Stage 1, 2, 3, 4, 5, 6, 7, 8"
#                  ,"\nhNPs: 1, 4, 8, 12 weeks"
#                  ,"\nMiller: VZ, SZ, IZ, SP, CP, MZ, SG"
#                  ,"\n"))
# ggsave(paste0(outGraphPfx, "Line_Graphs.pdf"), width = 14, height = 12)
################################################################################