
# Damon Polioudakis
# 2014-04-13
# Plot modules defined in Kang dataset across hNP and Kang time courses
################################################################################

rm(list = ls())

require("WGCNA")
require(ggplot2)
require(reshape2)

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

# hNPs
HNPAnnot = read.csv("../orig.data/HNPData1.4.8.12/annot.csv", row.names = 1)
HNPExpr = read.csv("../orig.data/HNPData1.4.8.12/exprdata.csv", row.names = 1)
HNPMeta = read.csv("../orig.data/HNPData1.4.8.12/sampleinfo.csv", row.names = 1)


## Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text=element_text(size=18)))
theme_update(plot.title = element_text(size = 16))


## Variables

graphCodeTitle <- "Kang_ModuleEG_Expression.R"
outGraphPfx <- "../analysis/graphs/Kang_ModuleEG_Expression_"


## Output Directory
dir.create("../analysis/graphs", recursive = TRUE)
################################################################################

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

## Limit to ENTREZ IDs present in all three datasets
## Arbitrary for Module preservation analysis, but mandatory for follow up
## procedures

if (TRUE) {
  ENTREZset = intersect(intersect(KangAnnotnet2$ENTREZ_ID
                              , MillerAnnot$ENTREZ_ID), HNPAnnot$ENTREZ_GENE_ID)
  
  Index = match(ENTREZset, KangAnnotnet2$ENTREZ_ID)
  KangExpr = KangExpr[Index, ]
  KangAnnot = KangAnnotnet2[Index, ]
  
  Index = match(ENTREZset, MillerAnnot$ENTREZ_ID)
  MillerExpr = MillerExpr[Index, ]
  MillerAnnot = MillerAnnot[Index, ]
  
  Index = match(ENTREZset, HNPAnnot$ENTREZ_GENE_ID)
  HNPExpr = HNPExpr[Index, ]
  HNPAnnot = HNPAnnot[Index, ]
}

if ((nrow(HNPExpr) == nrow(KangExpr)) == (nrow(HNPExpr) == nrow(MillerExpr)))
  cat('there are', nrow(HNPExpr), 'genes in all three datasets\n')
# there are 16451 genes in all three datasets
################################################################################

### Plot ME expression from Kang modules in datasets

## Calculate MEs
netColors <- net2$colors[match(ENTREZset, KangAnnotnet2$ENTREZ_ID)]
kangMEs <- moduleEigengenes(t(KangExpr), netColors)$eigengenes
hNPMEs <- moduleEigengenes(t(HNPExpr), netColors)$eigengenes
millerMEs <- moduleEigengenes(t(MillerExpr), netColors)$eigengenes

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

# Calculate mean ME expression at each time point or brain region
Dev_Time_Point_Mean_MEexpr <- function (timePoint, MEs) {
  df <- cbind(data.frame(Time_Point = timePoint), MEs)
  tmpL <- split(df, df$Time_Point)
  # Remove Time_Point column
  tmpL <- lapply(tmpL, function(df) df[ ,-1])
  outM <- sapply(tmpL, function(timePoint) {
    apply(timePoint, 2, mean)
  })
  outM
}

# Kang
ggM <- Dev_Time_Point_Mean_MEexpr(kangStage, kangMEs)
ggDF <- data.frame(t(ggM))
ggDF$Time_Point <- row.names(ggDF)
ggKangDF <- melt(ggDF, id.vars = "Time_Point"
                , variable.name = "ME", value.name = "Expression")
ggKangDF$Data_Set <- "Kang"

# hNPs
ggM <- Dev_Time_Point_Mean_MEexpr(hNPwK, hNPMEs)
ggDF <- data.frame(t(ggM))
ggDF$Time_Point <- row.names(ggDF)
gghNPdF <- melt(ggDF, id.vars = "Time_Point"
                 , variable.name = "ME", value.name = "Expression")
gghNPdF$Data_Set <- "hNP"

# Miller
ggM <- Dev_Time_Point_Mean_MEexpr(millerZones, millerMEs)
ggDF <- data.frame(t(ggM))
ggDF$Time_Point <- row.names(ggDF)
# Convert zones to numeric for plotting on Kang time course
ggDF$Time_Point <- as.numeric(factor(ggDF$Time_Point
                        , levels = c("VZ", "SZ", "IZ", "SP", "CP", "MZ", "SG")))
ggMillerdF <- melt(ggDF, id.vars = "Time_Point"
                , variable.name = "ME", value.name = "Expression")
ggMillerdF$Data_Set <- "Miller"

# Combine for ggplot2 plotting
ggDF <- rbind(ggKangDF, gghNPdF, ggMillerdF)
ggDF$Time_Point <- as.numeric(ggDF$Time_Point)

ggplot(ggDF, aes(x = Time_Point, y = Expression, color = Data_Set)) +
  facet_wrap(~ME) +
  geom_line() +
  xlab("Stage / Week / Zone") +
  ylab("ME Expression") +
  theme(axis.text.x = element_blank()) +
  ggtitle(paste0(graphCodeTitle
                 ,"\nKang Module Eigengene Expression"
                 ,"\nX-axis:"
                 ,"\nKang: Stage 1, 2, 3, 4, 5, 6, 7, 8"
                 ,"\nhNPs: 1, 4, 8, 12 weeks"
                 ,"\nMiller: VZ, SZ, IZ, SP, CP, MZ, SG"
                 ,"\n"))
ggsave(paste0(outGraphPfx, "Line_Graphs.pdf"), width = 14, height = 12)
######################################################################################################################################################
# WIP HERE
######################







######
#Module preservation algo
####

# convert tables such as columns are genes  and colnames are ENTREZ_IDS

KangExprt = t(KangExpr)
MillerExprt = t(MillerExpr)
HNPExprt = t(HNPExpr)


colnames(KangExprt) = KangAnnot$ENTREZ_ID
colnames(MillerExprt) = MillerAnnot$ENTREZ_ID
colnames(HNPExprt) = HNPAnnot$ENTREZ_GENE_ID

# Set 1 Compare Module Preservation Kang to Miller 

setLabels = c("KangData", "MillerLCMDE")

multiExpr = list(KangData = list(data = KangExprt), MillerLCMDE = list(data = MillerExprt)) 

multiColor = list(KangData = net2$colors[match(ENTREZset, KangAnnotnet2$ENTREZ_ID)])

dir.create(paste(home, args[2], "Modulepreservation", sep = "/"), showWarnings = F)

setwd(paste(home, args[2], "Modulepreservation", sep = "/"))

if(F){
  system.time( {
    mp = modulePreservation(multiExpr, multiColor, 
                            referenceNetworks = 1, 
                            nPermutations = 200, 
                            randomSeed = 1, 
                            quickCor = 0, 
                            verbose = 3)
  })
  
  #### Save the results
  save(mp, file  = "ModulePreservation_KangvsMiller.RData");
}

# Set 2 Compare Module Preservation Kang to HNPData 


setLabels = c("KangData", "HNPData")
multiExpr = list(KangData = list(data = KangExprt), HNPData = list(data = HNPExprt)) 

if(F){
  system.time( {
    mpKH = modulePreservation(multiExpr, multiColor, 
                              referenceNetworks = 1, 
                              nPermutations = 200, 
                              randomSeed = 1, 
                              quickCor = 0, 
                              verbose = 3)
  }
  )
  
  
  #### Save the results
  save(mpKH, file  = "ModulePreservation_KangvsHNP.RData");
}

setwd(home)

# plot the files 
if(T){
  load("../Output/Modulepreservation/modulePreservation_KangvsMiller.RData")
  load("../Output/Modulepreservation/ModulePreservation_KangvsHNP.RData")
  load("../orig.data/InVivoData/net2_power16_cutHeigt0.15.RData")
}

setwd(paste(home, args[2], "Modulepreservation", sep = "/"))


require(WGCNA)



#############################
#
# Plot Preservation 
#
#########################

# which modules to plot (all but named)

modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];

plotMods = !(modColors %in% c("grey", "gold"));

# plot Kang versus Miller


ref = 1
test = 2
statsObs  = cbind(mp$quality$observed[[ref]][[test]][, -1], 
                mp$preservation$observed[[ref]][[test]][, -1])
statsZ  = cbind(mp$quality$Z[[ref]][[test]][, -1], 
              mp$preservation$Z[[ref]][[test]][, -1]);
statslogp  = cbind(mp$quality$log.pBonf[[ref]][[test]][, -1], 
                 mp$preservation$log.pBonf[[ref]][[test]][, -1]);

print(cbind(statsObs[, c("medianRank.pres", "medianRank.qual")], 
            signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2), 
            signif(statslogp[, c("log.p.Bonfsummary.pres", "log.p.Bonfsummary.qual")], 3)))


text = modColors[plotMods]

refplotData = plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], 
                             mp$preservation$Z[[ref]][[test]][, 2], 
                             - mp$preservation$log.pBonf[[ref]][[test]][, 2])

mains = c("Preservation Median rank \nKang vs Miller", 
          "Preservation Zsummary \nKang vs Miller", 
          "Preservation -log10 Bonfp.val \nKang vs Miller")

pdf(file = "KangvsMiller-modulePreservation-Zsummary-medianRank.pdf", wi = 15, h = 5)


par(mfrow = c(1, 3))
par(mar = c(4.5, 6, 2.5, 1))
for (p in 1:3)
{min = min(refplotData[plotMods, p], na.rm = TRUE);
max = max(refplotData[plotMods, p], na.rm = TRUE);
# Adjust ploting ranges appropriately
if (p == 2 | p == 3)
{
  if (min > -max/10) min = -max/10
  ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
}
else
  ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))

plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21, 
     main = mains[p], 
     cex = 2.4, 
     ylab = mains[p], xlab = "Module size", log = "x", 
     ylim = ylim, 
     xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main  = 1.4)

labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 0.7, offs = 0.08)

# For Zsummary, add threshold lines
if (p == 2)
{
  abline(h = 0)
  abline(h = 2, col = "blue", lty = 2)
  abline(h = 10, col = "darkgreen", lty = 2)
}

if (p == 3)
{
  abline(h = 0, col = 1)
  abline(h = -log10(0.05), col = "blue", lty = 2)
  abline(h = -log10(0.01), col = "darkgreen", lty = 2)
}
}
# If plotting into a file, close it
dev.off();



# plot Kang versus HNP


ref = 1
test = 2
statsObs  = cbind(mpKH$quality$observed[[ref]][[test]][, -1], 
                mpKH$preservation$observed[[ref]][[test]][, -1])
statsZ  = cbind(mpKH$quality$Z[[ref]][[test]][, -1], 
              mpKH$preservation$Z[[ref]][[test]][, -1]);
statslogp  = cbind(mpKH$quality$log.pBonf[[ref]][[test]][, -1], 
                 mpKH$preservation$log.pBonf[[ref]][[test]][, -1]);

print(cbind(statsObs[, c("medianRank.pres", "medianRank.qual")], 
            signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2), 
            signif(statslogp[, c("log.p.Bonfsummary.pres", "log.p.Bonfsummary.qual")], 3)))



plotData = cbind(mpKH$preservation$observed[[ref]][[test]][, 2], 
                 mpKH$preservation$Z[[ref]][[test]][, 2], 
                 - mpKH$preservation$log.pBonf[[ref]][[test]][, 2])


modColors = rownames(mpKH$preservation$observed[[ref]][[test]])
moduleSizes = mpKH$preservation$Z[[ref]][[test]][, 1];


text = modColors[plotMods]

mains = c("Preservation Median rank \nKang vs HNP", 
          "Preservation Zsummary \nKang vs HNP", 
          "Preservation -log10 Bonfp.val \nKang vs HNP")


pdf(file = "KangvsHNP-modulePreservation-Zsummary-medianRank.pdf", wi = 15, h = 5)


par(mfrow = c(1, 3))
par(mar = c(4.5, 6, 2.5, 1))
for (p in 1:3)
{min = min(refplotData[plotMods, p], na.rm = TRUE);
max = max(refplotData[plotMods, p], na.rm = TRUE);
# Adjust ploting ranges appropriately
if (p == 2 | p == 3)
{
  if (min > -max/10) min = -max/10
  ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
}
else
  ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))

plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21, 
     main = mains[p], 
     cex = 2.4, 
     ylab = mains[p], xlab = "Module size", log = "x", 
     ylim = ylim, 
     xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main  = 1.4)

labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 0.7, offs = 0.08)

# For Zsummary, add threshold lines
if (p == 2)
{
  abline(h = 0)
  abline(h = 2, col = "blue", lty = 2)
  abline(h = 10, col = "darkgreen", lty = 2)
}

if (p == 3)
{
  abline(h = 0, col = 1)
  abline(h = -log10(0.05), col = "blue", lty = 2)
  abline(h = -log10(0.01), col = "darkgreen", lty = 2)
}
}
# If plotting into a file, close it
dev.off();



#####################################################################
# nStatistics of Module Preservation
#
#
#
#####################################################################


Prescomp  = data.frame(row.names = modColors, 
                     KangvsMillerZ =  mp$preservation$Z[[ref]][[test]][, 2], 
                     KangvsMillerlogp =  -mp$preservation$log.pBonf[[ref]][[test]][, 2], 
                     KangvsHNPZ = mpKH$preservation$Z[[ref]][[test]][, 2], 
                     KangvsHNPlogp =  -mpKH$preservation$log.pBonf[[ref]][[test]][, 2])

Prescomp$Ratio = (Prescomp[, 1]/Prescomp[, 2])

write.table(Prescomp, file = "../Output/Modulepreservation/Preservationscores.csv", sep = ", ")


MOI = modColors[which(Prescomp$KangvsMillerlogp > -log10(0.01) &
                      Prescomp$KangvsHNPlogp < -log10(0.05))]

MOI = MOI[!MOI%in%c("gold", "grey")]

cat("modules not preserved in HNP dataset:\n", MOI, "\n")

consMOI = modColors[which(Prescomp$KangvsMillerlogp > -log10(0.05)&
                          Prescomp$KangvsHNPlogp > -log10(0.05))]

consMOI = consMOI[!consMOI%in%c("gold", "grey")]
cat("modules preserved in HNP dataset:\n", consMOI, "\n")

ucMOI = modColors[!modColors %in% c(MOI, consMOI)]

ucMOI = ucMOI[!ucMOI%in%c("gold", "grey")]

cat("modules not picked in any dataset:\n", ucMOI, "\n") 


MEs = net2$MEs


pdf(file = "ModMEs_notPreserved.pdf", wi = 20, he = 8)
par(mfrow = c(1, 2), mar = c(4.5, 5, 4, 10))

limit = max(abs(MEs[, paste0("ME", MOI)]))*2

plot(0, 0, 
     xlim = c(min(sampleKey$Stage), max(sampleKey$Stage)), 
     ylim = c(-limit, limit), 
     xlab =  "Stage", 
     ylab = "Module Eigengene, norm to min", 
     main = "non-preserved Modules", 
     type = "n", 
     cex = 1.5)
legend(x = 8.5, y = 0.2, legend = MOI, pch = 16, col = MOI, bty = "n", xpd = T)
for (mod in MOI){
  set = MEs[, paste0("ME", mod)]
  set2 = spline(set~sampleKey$Stage, n = 100)
  set = set2$y-set2$y[1]
  lines(set2$x, set , col = mod, lwd = 3)
}


limit = max(abs(MEs[, paste0("ME", MOI)]))*2

plot(0, 0, 
     xlim = c(min(sampleKey$Stage), max(sampleKey$Stage)), 
     ylim = c(-limit, limit), 
     xlab =  "Stage", 
     ylab = "Module Eigengene", 
     main = "non-preserved Modules", 
     type = "n", 
     cex = 1.5)
legend(x = 8.5, y = 0.2, legend = MOI, pch = 16, col = MOI, bty = "n", xpd = T)
for (mod in MOI){
  set = MEs[, paste0("ME", mod)]
  set2 = spline(set~sampleKey$Stage, n = 100)
  lines(set2$x, set2$y , col = mod, lwd = 3)
}

dev.off()

npMOIOgenes = KangAnnotnet2[net2$colors %in% MOI, ]

npMOIOgenes.BKGD = KangAnnotnet2[net2$colors %in% c(MOI, consMOI), ]

write.table(npMOIOgenes, 
            file = paste(home, args[2], "Genelists/MP_not_cons_mods.csv", sep = "/"), 
            sep = ", ", row.names = F)

write.table(npMOIOgenes.BKGD, 
            file = paste(home, args[2], "Genelists/MP_not_cons_modsBKGD.csv", sep = "/"), 
            sep = ", ", row.names = F)


#save.image("ModulePreservation.RData")

#########################################################################
# Identifiy genes that have a different module membership
########################################################################


# select genes of conserved module


kMEBinder =  list()

MEKang = moduleEigengenes(KangExprt, multiColor[[1]])$eigengenes
MEMiller = moduleEigengenes(MillerExprt, multiColor[[1]])$eigengenes
MEHNP = moduleEigengenes(HNPExprt, multiColor[[1]])$eigengenes


for (module in consMOI) {
  colindex = multiColor[[1]] == module
  Memod = paste0("ME", module)
  Kangset = KangExprt[, colindex]
  Millerset = MillerExprt[, colindex]
  HNPset = HNPExprt[, colindex]
  kMEBinder[[module]] = data.frame(
    kMEKang = apply(Kangset, 2, function(x){cor(x, MEKang[, Memod], method = "spearman")}), 
    kMEMiller = apply(Millerset, 2, function(x){cor(x, MEMiller[, Memod], method = "spearman")}), 
    kMEHNP = apply(HNPset, 2, function(x){cor(x, MEHNP[, Memod], method = "spearman")}))
}

require(scatterplot3d)

pdf("kMEplots.pdf", 20, 20)

par(mfrow = c(ceiling(sqrt(length(consMOI))), ceiling(sqrt(length(consMOI)))), mar = c(3, 3, 3, 3))

for (module in consMOI){
  
  a = scatterplot3d(kMEBinder[[module]], color = "gray", pch = 21, 
                    type = "h", lty.hplot = 2, bg = module, 
                    xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1), 
                    angle = 210, main = module, cex.axis = 1.2)
  
  a.coords <- a$xyz.convert(kMEBinder[[module]])
  labels = KangAnnot$SYMBOL[match(rownames(kMEBinder[[module]]), KangAnnot$ENTREZ_ID)]
  text(a.coords$x, a.coords$y, labels = labels, pos = 4, cex = 0.5)
}

dev.off()


pdf("kMEdetails.pdf", 15, 8)
par(mfrow = c(2, 4), mar = c(5, 5, 8, 5), oma = c(1, 1, 5, 3))


#plot plots with critical genes named and generate results list  

module = "salmon"

for (module in consMOI){
  
  set2use = kMEBinder[[module]]
  refMean = mean(set2use$kMEKang)
  refSD = sd(set2use$kMEKang)
  set2use$SYMBOL = KangAnnot$SYMBOL[match(rownames(set2use), KangAnnot$ENTREZ_ID)] 
  set2use$ZkMEKang = (set2use$kMEKang-refMean)/refSD
  set2use$ZkMEMiller = (set2use$kMEMiller-refMean)/refSD
  set2use$ZkMEHNP = (set2use$kMEHNP-refMean)/refSD
  set2use$pkMEKang = p.adjust(2*(1-pnorm(abs(set2use$ZkMEKang))), n = nrow(set2use))
  set2use$pkMEMiller = p.adjust(2*(1-pnorm(abs(set2use$ZkMEMiller))), n = nrow(set2use))
  set2use$pkMEHNP = p.adjust(2*(1-pnorm(abs(set2use$ZkMEHNP))), n = nrow(set2use))
  
  set2use$kMEcrit =  set2use$pkMEMiller>0.1&
    set2use$pkMEHNP<0.05
  
  set2use$kMEbckgrd = set2use$pkMEMiller>0.1 |
    set2use$pkMEHNP>0.1
  
  a = scatterplot3d(set2use[, 1:3], color = "gray", pch = 21, 
                    type = "h", lty.hplot = 2, bg = module, 
                    mar = c(5, 5, 8, 7), 
                    xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1), 
                    angle = 210, main = "Module Membership")
  a.coords <- a$xyz.convert(set2use[, 1:3])
  labels = set2use$SYMBOL
  labels[!set2use$kMEcrit] = ""
  
  text(a.coords$x, a.coords$y, labels = labels, pos = 4, cex = 0.5)
  text(0, 6, paste("kMEs", module, "module"), xpd = T, cex = 2, font = 2)
  
  for (col in 1:3)
    hist(set2use[, col], breaks = 10, main = colnames(set2use)[col], xlim = c(-1, 1), col = module, 
         xlab = "kME score")
  kMEBinder[[module]] = set2use 
}

dev.off()



DkME.results = data.frame()


for(module in consMOI){
  set2use = kMEBinder[[module]]
  set2use$ENTREZ = rownames(set2use)
  res = set2use[set2use$kMEcrit == T, ]
  res$Module = rep(module, nrow(res))
  DkME.results = rbind(DkME.results, res)
}


DkME.results = DkME.results[, c(4, 13, 14, 1:3, 5:10)]

cat("non conserved genes : ", nrow(DkME.results), "\n", 
    "out of" , sum(unlist(multiColor) %in% consMOI), "genes in conserved modules\n")

cat(paste(unique(DkME.results$Module), tapply(DkME.results$ENTREZ, DkME.results$Module, length), "\n"))


summarytab = data.frame(module = unique(DkME.results$Module))
summarytab$Genes = sapply(summarytab$module, function(x){sum(unlist(multiColor) == x)})    
summarytab$notPreGenes = as.vector(tapply(DkME.results$ENTREZ, DkME.results$Module, length))
summarytab$ratio = summarytab$notPreGenes/summarytab$Genes

write.table(summarytab, file = paste(home, args[2], "Modulepreservation/Summary.csv", sep = "/"), sep = ", ", row.names = F)

write.table(DkME.results, 
            file = paste(home, args[2], "Genelists", "DkME_targets.csv", sep = "/"), 
            sep = ", ", row.names = F)

DkME.bckgrd = data.frame()


for(module in consMOI){
  set2use = kMEBinder[[module]]
  set2use$ENTREZ = rownames(set2use)
  res = set2use[set2use$kMEbckgrd == T, ]
  res$Module = rep(module, nrow(res))
  DkME.bckgrd = rbind(DkME.bckgrd, res)
}

DkME.bckgrd = DkME.bckgrd[, c(4, 13, 12, 14, 1:3, 5:10)]

write.table(DkME.bckgrd, 
            file = paste(home, args[2], "Genelists", "DkME_BKGD.csv", sep = "/"), 
            sep = ", ", row.names = F)  

save.image("ModulePreservation.RData")




