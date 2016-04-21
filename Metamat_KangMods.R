rm(list=ls())
options(stringsAsFactors=FALSE)
setwd("/geschwindlabshares/Hi-C_Disease/Schizophrenia/geneenrichment")

library(biomaRt) ## First, construct a matrix containing all possible protein coding genes. We will use human gene symbols. Although other identifiers, such as ENSG IDs are preferred, many studies only report human gene symbols, so we will use this as the common identifier.
library(WGCNA) ## For plotting heatmaps

## First, we get all the genes from Gencode v19
#getinfo <- c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position","strand","band","gene_biotype")
#mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="feb2014.archive.ensembl.org") ## Using Gencode v19 annotations
#geneAnno1 <- getBM(attributes = getinfo,filters=c("chromosome_name"),values=c(seq(1,22,by=1),"X"),mart=mart)
#geneAnno1 <- geneAnno1[geneAnno1[,"gene_biotype"]=="protein_coding",]; geneAnno1 <- geneAnno1[!duplicated(geneAnno1[,"hgnc_symbol"]) & geneAnno1[,"hgnc_symbol"]!="",] ## 19163 with hgnc symbols after removing duplicates. We use gene symbols here because the RDNV data comes in gene symbol format. Otherwise, I prefer to use ENSG IDs.


## Get mouse orthology information
#getinfo <- c("ensembl_gene_id","mmusculus_homolog_ensembl_gene","mmusculus_homolog_orthology_type")
#m.gene.anno <- getBM(attributes = getinfo,filters="ensembl_gene_id",values=geneAnno1[,1],mart=mart)
#m.gene.1to1 <- unique(m.gene.anno[m.gene.anno[,"mmusculus_homolog_orthology_type"]=="ortholog_one2one","ensembl_gene_id"])
#m.gene.1to1 <- geneAnno1[match(m.gene.1to1,geneAnno1[,"ensembl_gene_id"]),"hgnc_symbol"]
#m.gene.1to1 <- unique(m.gene.1to1[m.gene.1to1!=""])

geneA = load("/geschwindlabshares/atacseq/scripts/geneAnnotation_biomart.rda");

## Construct a matrix to store the gene lists - this starts as a matrix with gene information and list membership information, and new lists are appended. If the new list contains the gene, it gets a 1, if not 0, if the gene is not even in the study (not in the background) we mark an NA
metaMat <- matrix(NA,ncol=1,nrow=nrow(geneAnno1)) ## We can add columns as necessary
rownames(metaMat) <- geneAnno1[,"hgnc_symbol"]
ENSGIDs <- geneAnno1[,"ensembl_gene_id"]
geneLength <- geneAnno1[,"end_position"]-geneAnno1[,"start_position"]
names(geneLength) <- geneAnno1[,"hgnc_symbol"] ## For use later in making the covariate matrix
metaX = geneAnno1[match(rownames(metaMat), geneAnno1$hgnc_symbol),"chromosome_name"] 

## a) Gene lists from Luis 
ORAgenelist <- read.csv("/geschwindlabshares/Hi-C/traitsenrichment/data/Lists/GO_ORAGeneLists.csv", fill=T, header=T)
for(i in 1:length(table(ORAgenelist$SortClass))){
  thislist <- ORAgenelist[ORAgenelist$SortClass==i,1]
  listname <- gsub(" ", "_",unique(ORAgenelist[ORAgenelist$SortClass==i,4]))
  colnames(metaMat)[dim(metaMat)[2]] <- listname
  matchlist <- match(rownames(metaMat), thislist)
  metaMat[!is.na(matchlist),listname] <- 1
  metaMat[is.na(matchlist),listname] <- 0 
  metaMat <- cbind(metaMat,rep(NA,nrow(metaMat)))
}
metaMat = metaMat[, -dim(metaMat)[2]] ##Remove NA in last column

## Cell types in fetal brain (Pollen et al S3)
load(file="/geschwindlabshares/atacseq/HumanEnhancers/Pollen-2015-SingleCell/TableS3-CellTypesGenes-ENSG/PollenS3GenesENSG.R");
for(i in 1:length(PollenS3GenesENSG)){
    thislist = PollenS3GenesENSG[[i]]
    gename = geneAnno1[geneAnno1$ensembl_gene_id %in% thislist, 2]
    gename = gename[gename!=""]
    matchlist = match(rownames(metaMat), gename)
    matchlist = ifelse(is.na(matchlist),0,1)
    metaMat = cbind(metaMat, matchlist)
    colnames(metaMat)[dim(metaMat)[2]] = names(PollenS3GenesENSG)[i]
}

## Cell types in fetal brain (Pollen et al S4)
load(file="/geschwindlabshares/atacseq/HumanEnhancers/Pollen-2015-SingleCell/HumanRegGenesbyCellType/Non-Unique/PollenCellTypesGenes.R");
for(i in 1:length(CellTypesGenes)){
    thislist = CellTypesGenes[[i]]
    gename = geneAnno1[geneAnno1$ensembl_gene_id %in% thislist, 2]
    gename = gename[gename!=""]
    matchlist = match(rownames(metaMat), gename)
    matchlist = ifelse(is.na(matchlist),0,1)
    metaMat = cbind(metaMat, matchlist)
    colnames(metaMat)[dim(metaMat)[2]] = names(CellTypesGenes)[i]
}

## LCM DE (Miller et al) Use this for all possible layer comparisons

#dir = "/geschwindlabshares/CoNTExT/2012-198/2012-198GeschwindHT-12/AllenLCM/RRHO";
#dir = "/geschwindlabshares/CoNTExT/2012-198/2012-198GeschwindHT-12/scripts/WebsiteCode/v2/LCMDE";
#names = c("VZ","SZ","IZ","SP","CPi","CPo","MZ","SG");

#for(i in 1:7){
 #   for(j in 2:8){
  #       if (i>=j) {
   #      next
    #}
     #    thislist = read.csv(paste(dir,"/LCM_DE_",names[i],"vs",names[j],".csv",sep=""));
      #   DEup = which(thislist$Pval.FDR<=0.05 & thislist$LogFC>0);
       #  DEdown = which(thislist$Pval.FDR<=0.05 & thislist$LogFC<0);
        # DEup = thislist[DEup,2];
         #DEdown = thislist[DEdown,2];
         
#         matchlist = match(rownames(metaMat), DEup);
 #        matchlist = ifelse(is.na(matchlist),0,1);
  #       metaMat = cbind(metaMat, matchlist);
   #      colnames(metaMat)[dim(metaMat)[2]] = paste("LCM-DE-Up-",names[i],"_",names[j],sep="");
#
 #        matchlist = match(rownames(metaMat), DEdown);
  #       matchlist = ifelse(is.na(matchlist),0,1);
   #      metaMat = cbind(metaMat, matchlist);
    #     colnames(metaMat)[dim(metaMat)[2]] = paste("LCM-DE-Down-",names[i],"_",names[j],sep="");   
     #}}

## LCM DE (Miller et al) Use for biologically-relevant layer transitions only

#dir = "/geschwindlabshares/CoNTExT/2012-198/2012-198GeschwindHT-12/AllenLCM/RRHO";
dir = "/geschwindlabshares/CoNTExT/2012-198/2012-198GeschwindHT-12/scripts/WebsiteCode/v2/LCMDE";
names = c("VZvsSZ", "SZvsIZ", "IZvsSP","IZvsCPi","IZvsCPo","SPvsCPi", "SPvsCPo");

for(i in 1:7){
         thislist = read.csv(paste(dir,"/LCM_DE_",names[i],".csv",sep=""));
         DEup = which(thislist$Pval.FDR<=0.05 & thislist$LogFC>0);
         DEdown = which(thislist$Pval.FDR<=0.05 & thislist$LogFC<0);
         DEup = thislist[DEup,2];
         DEdown = thislist[DEdown,2];
         
         matchlist = match(rownames(metaMat), DEup);
         matchlist = ifelse(is.na(matchlist),0,1);
         metaMat = cbind(metaMat, matchlist);
         colnames(metaMat)[dim(metaMat)[2]] = paste("LCM-DE-Up-",names[i],sep="");

         matchlist = match(rownames(metaMat), DEdown);
         matchlist = ifelse(is.na(matchlist),0,1);
         metaMat = cbind(metaMat, matchlist);
         colnames(metaMat)[dim(metaMat)[2]] = paste("LCM-DE-Down-",names[i],sep="");   
     }


## b) Gene lists from candidate studies and exome sequencing
metaMat = cbind(metaMat, rep(NA,nrow(metaMat))) #Re-added the NAs at the end of metaMat
thislist <- unique(read.csv("/geschwindlabshares/Hi-C/traitsenrichment/data/Lists/ConstrainedGenes.csv")[,"gene"])
listname <- "Constrained_Genes"
colnames(metaMat)[ncol(metaMat)] <- listname
matchlist <- match(rownames(metaMat),thislist)
metaMat[!is.na(matchlist),listname] <- 1
metaMat[is.na(matchlist),listname] <- 0

## Iossifov ASD
listMat <- read.csv("/geschwindlabshares/Hi-C/traitsenrichment/data/Lists/AllDeNovoByGene.csv") ## From the supplmental tables of Iossifov et al., 2014 - note that some gene symbols have the excel conversion error in this data, but it will get removed when we overlap things
listMat <- listMat[!duplicated(listMat[,1]),]
rownames(listMat) <- listMat[,"gene"]
exomeLength <- listMat[,"codingLenInTarget"] ## For later use in making the covariate matrix
names(exomeLength) <- rownames(listMat)
listMat <- listMat[,c(5:9,12,17:18,23:24,29)] ## Only keeping the desired columns
listMat[listMat>1] <- 1 ## If there is more than one count, just count it as 1 - we are using "yes" vs "no" criteria
metaMat <- cbind(metaMat,listMat[match(rownames(metaMat),rownames(listMat)),])
colnames(metaMat)[(ncol(metaMat)-5):ncol(metaMat)] <- paste("ASD_Iossifov",colnames(metaMat)[(ncol(metaMat)-5):ncol(metaMat)], sep="_")

iossmat = listMat[,c(6,7)]
iossmat$gene = rownames(iossmat)

## de Rubeis ASD
listMat <- read.table("/geschwindlabshares/Hi-C/traitsenrichment/data/Lists/deRubeis_mutationlist.txt", header=T, sep="\t")
rownames(listMat) = listMat$Gene
listMat = listMat[,3:7]
colnames(listMat) = paste0("Rubeis_", colnames(listMat))
metaMat <- cbind(metaMat,listMat[match(rownames(metaMat),rownames(listMat)),])
rubmat = listMat
rubmat$gene = rownames(rubmat)

## iossifov + de Rubeis
asdmat = merge(rubmat, iossmat, by=intersect("gene", "gene"))
asdmat$combined_dn_prob_LOF = asdmat$dnv_LGDs_prb + asdmat$Rubeis_dn.LoF
asdmat$combined_dn_sib_LOF = asdmat$dnv_LGDs_sib
asdmat$combined_prob_LOF = asdmat$Rubeis_dn.LoF + asdmat$Rubeis_case.LoF + asdmat$Rubeis_trans.LoF + asdmat$Rubeis_ntrans.LoF + asdmat$dnv_LGDs_prb
asdmat$combined_sib_LOF = asdmat$Rubeis_ctrl.LoF + asdmat$dnv_LGDs_sib
listMat = asdmat[,c(9:12)]
rownames(listMat) = asdmat$gene
colnames(listMat) = paste("ASD",colnames(listMat),sep="_")

metaMat <- cbind(metaMat,listMat[match(rownames(metaMat),rownames(listMat)),])

## de novo ID: NEJM + Lancet
listMat <- read.table("/geschwindlabshares/Hi-C/traitsenrichment/data/Lists/ID_denovo_deLigt_NEJM.txt", header=T, sep="\t", fill=T)
nejmmat = listMat[listMat$nature_of_mutation=="D",1]
listMat <- read.table("/geschwindlabshares/Hi-C/traitsenrichment/data/Lists/ID_denovo_Rauch_Lancet.txt", header=T, sep="\t", fill=T)
lancetmat = listMat[listMat$Type %in% c("frameshift", "nonsense", "splice"),1]
thislist = unique(c(lancetmat, nejmmat))
listname <- "DeNovoLGDsInID"

metaMat <- cbind(metaMat,rep(NA,nrow(metaMat)))
colnames(metaMat)[ncol(metaMat)] <- listname
matchlist <- match(rownames(metaMat), thislist)
metaMat[!is.na(matchlist), listname] <- 1
metaMat[is.na(matchlist), listname] <- 0 

## c) Laminae and cell types
listMat <- read.csv("/geschwindlabshares/Hi-C/traitsenrichment/data/Laminae/FetalLaminarAssignmentWithTtest.csv") ## From Miller et al. 2014 - after running DGE in this list, we threshold for expression above other layers in each layer, and use significant events as the gene list for that layer
listMat <- listMat[!duplicated(listMat[,"gene_symbol"]),]
rownames(listMat) <- listMat[,"gene_symbol"]
pMat <- listMat[,c(26:32)]
bMat <- listMat[,c(19:25)]
tmpMat <- matrix(p.adjust(as.numeric(data.matrix(pMat)),method="BH"),nrow=nrow(pMat),ncol=ncol(pMat))
rownames(tmpMat) <- rownames(listMat)
colnames(tmpMat) <- colnames(pMat)
listMat <- tmpMat
listMat[tmpMat<=0.05 & bMat>0] <- 1
listMat[!(tmpMat<=0.05 & bMat>0)] <- 0
metaMat <- cbind(metaMat,listMat[match(rownames(metaMat),rownames(listMat)),])

listMat <- read.csv("/geschwindlabshares/Hi-C/traitsenrichment/data/Laminae/LaminarAssignmentWithTtest.csv") ## From Bernard et al., 2012 - similar treatment as the fetal list above
listMat <- listMat[!duplicated(listMat[,"hgnc_symbol"]),]
rownames(listMat) <- listMat[,"hgnc_symbol"]
pMat <- listMat[,c(18:22)]
bMat <- listMat[,c(13:17)]
tmpMat <- matrix(p.adjust(as.numeric(data.matrix(pMat)),method="BH"),nrow=nrow(pMat),ncol=ncol(pMat))
rownames(tmpMat) <- rownames(listMat)
colnames(tmpMat) <- colnames(pMat)
listMat <- tmpMat
listMat[tmpMat<=0.05 & bMat>0] <- 1
listMat[!(tmpMat<=0.05 & bMat>0)] <- 0
metaMat <- cbind(metaMat,listMat[match(rownames(metaMat),rownames(listMat)),])

listMat <- read.csv("/geschwindlabshares/Hi-C/traitsenrichment/data/CellTypes/ZhangEtAlCellTypeAssignmentWithTtest.csv") ## From Zhang et al., 2014 - similar treatment as layer lists above, using the RNA-seq data from Zhang et al., 2014
listMat <- listMat[!duplicated(listMat[,"hgnc_symbol"]) & !is.na(listMat[,"hgnc_symbol"]) & listMat[,"hgnc_symbol"]!="",]
rownames(listMat) <- listMat[,"hgnc_symbol"]
pMat <- listMat[,c(37,33,36,35,34)]
bMat <- listMat[,c(21,17,20,19,18)]
tmpMat <- matrix(p.adjust(as.numeric(data.matrix(pMat)),method="BH"),nrow=nrow(pMat),ncol=ncol(pMat))
rownames(tmpMat) <- rownames(listMat)
colnames(tmpMat) <- colnames(pMat)
listMat <- tmpMat
listMat[tmpMat<=0.05 & bMat>0] <- 1
listMat[!(tmpMat<=0.05 & bMat>0)] <- 0
metaMat <- cbind(metaMat,listMat[match(rownames(metaMat),rownames(listMat)),])
colnames(listMat) <- paste("2fold",colnames(listMat),sep=".")
listMat[tmpMat<=0.05 & bMat>1] <- 1
listMat[!(tmpMat<=0.05 & bMat>1)] <- 0
metaMat <- cbind(metaMat,listMat[match(rownames(metaMat),rownames(listMat)),])

listMat <- read.csv("/geschwindlabshares/Hi-C/traitsenrichment/data/CellTypes/DoyleEtAlCorticalCellTypeAssignmentWithTtest.csv") ## From Doyle et al., 2009 - treated similarly as the data in Zhang et al. 2014
listMat <- listMat[!duplicated(listMat[,"hgnc_symbol"]) & !is.na(listMat[,"hgnc_symbol"]) & listMat[,"hgnc_symbol"]!="",]
rownames(listMat) <- listMat[,"hgnc_symbol"]
pMat <- listMat[,c(37:46)]
bMat <- listMat[,c(17:26)]
tmpMat <- matrix(p.adjust(as.numeric(data.matrix(pMat)),method="BH"),nrow=nrow(pMat),ncol=ncol(pMat))
rownames(tmpMat) <- rownames(listMat)
colnames(tmpMat) <- colnames(pMat)
listMat <- tmpMat
listMat[tmpMat<=0.05 & bMat>0] <- 1
listMat[!(tmpMat<=0.05 & bMat>0)] <- 0
metaMat <- cbind(metaMat,listMat[match(rownames(metaMat),rownames(listMat)),])
colnames(listMat) <- paste("2fold",colnames(listMat),sep=".")
listMat[tmpMat<=0.05 & bMat>1] <- 1
listMat[!(tmpMat<=0.05 & bMat>1)] <- 0
##listMat[apply(listMat,1,sum)>1,] <- 0
metaMat <- cbind(metaMat,listMat[match(rownames(metaMat),rownames(listMat)),])

## d) Transcriptional regulators and RBPs
listMat <- read.csv("/geschwindlabshares/Hi-C/traitsenrichment/data/Regulators/CHD8_binding.csv") ## CHD8 targets from Sugathan et al., 2014
colnames(listMat) <- listMat[1,]
listMat <- listMat[-c(1),]
listMat <- listMat[!duplicated(listMat[,"gene symbol"]),]
genenames <- as.character(listMat[,"gene symbol"])
listMat <- as.numeric(data.matrix(listMat[c(4)]))
names(listMat) <- genenames
metaMat <- cbind(metaMat,listMat[match(rownames(metaMat),names(listMat))])
colnames(metaMat)[ncol(metaMat)] <- "CHD8.target"

listMat <- read.table("/geschwindlabshares/Hi-C/traitsenrichment/data/Regulators/DarnellEtAl_GENCODEV10_FMRP_one2one_orthologs.txt",sep="\t") ## FMRP targets from Darnell et al., 2011
listMat <- listMat[!duplicated(listMat[,6]),6]
backMat <- m.gene.1to1
targets <- as.numeric(!is.na(match(backMat,listMat)))
names(targets) <- backMat
metaMat <- cbind(metaMat,targets[match(rownames(metaMat),names(targets))])
colnames(metaMat)[ncol(metaMat)] <- "FMRP.target"

listMat <- read.table("/geschwindlabshares/Hi-C/traitsenrichment/data/Regulators/CellReports_GENCODEV10_RBFOX1_one2one_orthologs.txt",sep="\t") ## Rbfox1 targets from Weyn-Vanhentenryck et al., 2014
listMat <- listMat[!duplicated(listMat[,15]),15]
backMat <- m.gene.1to1
targets <- as.numeric(!is.na(match(backMat,listMat)))
names(targets) <- backMat
metaMat <- cbind(metaMat,targets[match(rownames(metaMat),names(targets))])
colnames(metaMat)[ncol(metaMat)] <- "RBFOX1.target"

datLists <- read.table("/geschwindlabshares/Hi-C/traitsenrichment/data/Lists/CandidateGeneLists.txt",sep="\t",header=TRUE)

metaMat <- cbind(metaMat, rep(NA, nrow(metaMat)))
thislist <- datLists[, "FOXP2targetsPFSKSKNMC"]
listname <- "FOXP2target"
colnames(metaMat)[ncol(metaMat)] <- listname
matchlist <- match(rownames(metaMat), thislist)
metaMat[!is.na(matchlist), listname] <- 1
metaMat[is.na(matchlist), listname] <- 0

metaMat <- cbind(metaMat, rep(NA, nrow(metaMat)))
thislist <- datLists[, "FOXP2_UP_human"]
listname <- "FOXP2UpHuman"
colnames(metaMat)[ncol(metaMat)] <- listname
matchlist <- match(rownames(metaMat), thislist)
metaMat[!is.na(matchlist), listname] <- 1
metaMat[is.na(matchlist), listname] <- 0

metaMat <- cbind(metaMat, rep(NA, nrow(metaMat)))
thislist <- datLists[, "FOXP2_DOWN_human"]
listname <- "FOXP2DownHuman"
colnames(metaMat)[ncol(metaMat)] <- listname
matchlist <- match(rownames(metaMat), thislist)
metaMat[!is.na(matchlist), listname] <- 1
metaMat[is.na(matchlist), listname] <- 0 

## e) Network studies
listMat <- read.table("/geschwindlabshares/Hi-C/traitsenrichment/data/Networks/VoineaguEtAl_GENCODEV10_M12.txt",sep="\t") ## asdM12
listMat <- listMat[!duplicated(listMat[,1]),1]
backMat <- read.table("/geschwindlabshares/Hi-C/traitsenrichment/data/Networks/Voineagu_kmes_geneInfoCTX.csv",sep=",",header=TRUE) ## Background
backMat <- backMat[!duplicated(backMat[,2]),2]
targets <- as.numeric(!is.na(match(backMat,listMat)))
names(targets) <- backMat
metaMat <- cbind(metaMat,targets[match(rownames(metaMat),names(targets))])
colnames(metaMat)[ncol(metaMat)] <- "asdM12"

listMat <- read.table("/geschwindlabshares/Hi-C/traitsenrichment/data/Networks/VoineaguEtAl_GENCODEV10_M16.txt",sep="\t") ## asdM16
listMat <- listMat[!duplicated(listMat[,1]),1]
backMat <- read.table("/geschwindlabshares/Hi-C/traitsenrichment/data/Networks/Voineagu_kmes_geneInfoCTX.csv",sep=",",header=TRUE) ## Background
backMat <- backMat[!duplicated(backMat[,2]),2]
targets <- as.numeric(!is.na(match(backMat,listMat)))
names(targets) <- backMat
metaMat <- cbind(metaMat,targets[match(rownames(metaMat),names(targets))])
colnames(metaMat)[ncol(metaMat)] <- "asdM16"

## Autism co-expression and differential expression gene lists - from new data
datCpat <- read.csv("/geschwindlabshares/Hi-C/traitsenrichment/data/ASD_RNAseq/CPAT_FCvTCinASDandCTL.csv")
rownames(datCpat) <- datCpat[,2]
datCpat <- datCpat[,-c(1,2)]
FCvTCcpat <- as.numeric(datCpat[,"ASD.FDR"]>=0.05 & datCpat[,"CTL.FDR"]<0.05)
names(FCvTCcpat) <- datCpat[,"hgnc_symbol"]
metaMat <- cbind(metaMat,FCvTCcpat[match(rownames(metaMat),names(FCvTCcpat))])
colnames(metaMat)[ncol(metaMat)] <- "Attenuated.Patterning"

## Network and DGE data
datNet <- read.table("/geschwindlabshares/Hi-C/traitsenrichment/data/ASD_RNAseq/CTXgeneAssociationSummaries.txt",header=TRUE)
UpInASDCTX <- as.numeric(p.adjust(datNet[,"p.condition"],method="BH") < 0.05 & datNet[,"beta.condition"]>0)
DownInASDCTX <- as.numeric(p.adjust(datNet[,"p.condition"],method="BH") < 0.05 & datNet[,"beta.condition"]<0)
names(UpInASDCTX) <- names(DownInASDCTX) <- datNet[,"hgnc_symbol"]
metaMat <- cbind(metaMat,UpInASDCTX[match(rownames(metaMat),names(UpInASDCTX))])
colnames(metaMat)[ncol(metaMat)] <- "UpInASD.CTX"
metaMat <- cbind(metaMat,DownInASDCTX[match(rownames(metaMat),names(DownInASDCTX))])
colnames(metaMat)[ncol(metaMat)] <- "DownInASD.CTX"

modNames <- datNet[,"moduleLabels"]
names(modNames) <- datNet[,"hgnc_symbol"]
for (i in c(c(1:17))) { ## Using all modules - can use only ASD related modules
  thesegenes <- names(modNames)[modNames==i]
  thesegenes <- thesegenes[thesegenes!=""]
  targets <- as.numeric(!is.na(match(names(modNames),thesegenes)))
  names(targets) <- names(modNames)
  metaMat <- cbind(metaMat,targets[match(rownames(metaMat),names(targets))])
  colnames(metaMat)[ncol(metaMat)] <- paste("ASD.RNAseq.M",i,sep="")
}

## f) Neel's module 
listMat <- read.table("/geschwindlabshares/Hi-C/traitsenrichment/data/Networks/NeelModule.txt",header=FALSE)
backMat = listMat[,2]
for(i in c(1:6,8:18)){
  modulename = paste("M",i,sep="")
  modulemat = listMat[listMat[,3]==modulename, 2]
  modulemat = modulemat[!duplicated(modulemat) & !is.na(modulemat)]
  targets <- as.numeric(!is.na(match(backMat, modulemat)))
  names(targets) <- backMat
  metaMat <- cbind(metaMat,targets[match(rownames(metaMat),names(targets))])
  colnames(metaMat)[ncol(metaMat)] <- paste("Parikshak_", modulename, sep="")
}

## g) SCZ de novo mutations: Fromer et al., 2014 Nature 
# 1-1. de novo mutations: LOF vs missense
#denovo = read.table("/home/hwon/bio/Hi-C/disease/Schizophrenia/denovo/SCZ_denovo.txt", sep="\t", header=T) #Permission denied
denovo = read.table("/geschwindlabshares/atacseq/Enrichment_Lists/SCZ_denovo.txt", sep="\t", header=T)
sczall = denovo$Genes
silmatch = grep(pattern="silent", denovo$Gene.annotations)
nonsil = denovo[-silmatch,] # nonsilent 

sczlof = nonsil[nonsil$Gene.annotations %in% c("frameshift", "nonsense", "codon-deletion", "codon-insertion", "esplice"),5]
sczmatch = sczlof[match(rownames(metaMat),sczlof)]
sczlist = ifelse(is.na(sczmatch), 0, 1)
metaMat <- cbind(metaMat,sczlist)
colnames(metaMat)[ncol(metaMat)] <- "SCZ_Fromer_LGD"

sczmis = nonsil[nonsil$Gene.annotations=="missense",5]
sczmatch = sczmis[match(rownames(metaMat),sczmis)]
sczlist = ifelse(is.na(sczmatch), 0, 1)
metaMat <- cbind(metaMat,sczlist)
colnames(metaMat)[ncol(metaMat)] <- "SCZ_Fromer_mis"

# 1-3. SCZ DEG: microarray
#sczdeg = read.csv("/geschwindlabshares/Hi-C/traitsenrichment/data/SCZ_DEG/scz_meta-analysis_Chen_Maycox_Lanz_Iwa_Narayan_062615.csv") #Permission denied
sczdeg = read.csv("/geschwindlabshares/atacseq/Enrichment_Lists/scz_meta-analysis_Chen_Maycox_Lanz_Iwa_Narayan_062615.csv")
sczdeg$FDR = p.adjust(sczdeg$p, method="BH")

deg0.05 = sczdeg[sczdeg$FDR<0.05, 1]
deg0.01 = sczdeg[sczdeg$FDR<0.01, 1]

dr5 = sczdeg[sczdeg$FDR<0.05 & sczdeg$B<0, 1]
ur5 = sczdeg[sczdeg$FDR<0.05 & sczdeg$B>0, 1]

dr1 = sczdeg[sczdeg$FDR<0.01 & sczdeg$B<0, 1]
ur1 = sczdeg[sczdeg$FDR<0.01 & sczdeg$B>0, 1]

drfc5 = sczdeg[sczdeg$FDR<0.05 & sczdeg$B>0.2, 1]
urfc5 = sczdeg[sczdeg$FDR<0.05 & sczdeg$B<(-0.2),1]

geneA = geneAnno1[geneAnno1$ensembl_gene_id %in% sczdeg$X, 2]

gene = geneAnno1[geneAnno1$ensembl_gene_id %in% deg0.05, 2]
gene = unique(gene)
gene = gene[gene!=""]
sczmatch = gene[match(rownames(metaMat),gene)]
sczlist = ifelse(is.na(sczmatch), 0, 1)
sczlist[!(rownames(metaMat) %in% geneA)] <- NA
metaMat <- cbind(metaMat,sczlist)
colnames(metaMat)[ncol(metaMat)] <- "SCZ_microarray_DEG_0.05"

gene = geneAnno1[geneAnno1$ensembl_gene_id %in% dr5, 2]
gene = unique(gene)
gene = gene[gene!=""]
sczmatch = gene[match(rownames(metaMat),gene)]
sczlist = ifelse(is.na(sczmatch), 0, 1)
sczlist[!(rownames(metaMat) %in% geneA)] <- NA
metaMat <- cbind(metaMat,sczlist)
colnames(metaMat)[ncol(metaMat)] <- "SCZ_microarray_downreg_0.05"

gene = geneAnno1[geneAnno1$ensembl_gene_id %in% ur5, 2]
gene = unique(gene)
gene = gene[gene!=""]
sczmatch = gene[match(rownames(metaMat),gene)]
sczlist = ifelse(is.na(sczmatch), 0, 1)
sczlist[!(rownames(metaMat) %in% geneA)] <- NA
metaMat <- cbind(metaMat,sczlist)
colnames(metaMat)[ncol(metaMat)] <- "SCZ_microarray_upreg_0.05"

gene = geneAnno1[geneAnno1$ensembl_gene_id %in% deg0.01, 2]
gene = unique(gene)
gene = gene[gene!=""]
sczmatch = gene[match(rownames(metaMat),gene)]
sczlist = ifelse(is.na(sczmatch), 0, 1)
sczlist[!(rownames(metaMat) %in% geneA)] <- NA
metaMat <- cbind(metaMat,sczlist)
colnames(metaMat)[ncol(metaMat)] <- "SCZ_microarray_DEG_0.01"

gene = geneAnno1[geneAnno1$ensembl_gene_id %in% dr1, 2]
gene = unique(gene)
gene = gene[gene!=""]
sczmatch = gene[match(rownames(metaMat),gene)]
sczlist = ifelse(is.na(sczmatch), 0, 1)
sczlist[!(rownames(metaMat) %in% geneA)] <- NA
metaMat <- cbind(metaMat,sczlist)
colnames(metaMat)[ncol(metaMat)] <- "SCZ_microarray_downreg_0.01"

gene = geneAnno1[geneAnno1$ensembl_gene_id %in% ur1, 2]
gene = unique(gene)
gene = gene[gene!=""]
sczmatch = gene[match(rownames(metaMat),gene)]
sczlist = ifelse(is.na(sczmatch), 0, 1)
sczlist[!(rownames(metaMat) %in% geneA)] <- NA
metaMat <- cbind(metaMat,sczlist)
colnames(metaMat)[ncol(metaMat)] <- "SCZ_microarray_upreg_0.01"

# 1-4. SCZ DEG: RNA-seq
#sczdeg = read.csv("/geschwindlabshares/Hi-C_Disease/Schizophrenia/scz_DEG/commonMind/DGE_SCZ_commondMind.cqn.balanced.regressFirst.012716.csv") #Permission denied
sczdeg = read.csv("/geschwindlabshares/atacseq/Enrichment_Lists/DGE_SCZ_commondMind.cqn.balanced.regressFirst.012716.csv")
deg0.001 = sczdeg[sczdeg$adj.p.value<0.001, 1]
deg0.01 = sczdeg[sczdeg$adj.p.value<0.01, 1]
deg0.05 = sczdeg[sczdeg$adj.p.value<0.05, 1]

dr0.001 = sczdeg[sczdeg$adj.p.value<0.001 & sczdeg$log2FC<0, 1]
ur0.001 = sczdeg[sczdeg$adj.p.value<0.001 & sczdeg$log2FC>0, 1]

dr0.01 = sczdeg[sczdeg$adj.p.value<0.01 & sczdeg$log2FC<0, 1]
ur0.01 = sczdeg[sczdeg$adj.p.value<0.01 & sczdeg$log2FC>0, 1]

dr0.05 = sczdeg[sczdeg$adj.p.value<0.05 & sczdeg$log2FC<0, 1]
ur0.05 = sczdeg[sczdeg$adj.p.value<0.05 & sczdeg$log2FC>0, 1]

geneA = geneAnno1[geneAnno1$ensembl_gene_id %in% sczdeg$X, 2]

gene = geneAnno1[geneAnno1$ensembl_gene_id %in% deg0.001, 2]
gene = unique(gene)
gene = gene[gene!=""]
sczmatch = gene[match(rownames(metaMat),gene)]
sczlist = ifelse(is.na(sczmatch), 0, 1)
sczlist[!(rownames(metaMat) %in% geneA)] <- NA
metaMat <- cbind(metaMat,sczlist)
colnames(metaMat)[ncol(metaMat)] <- "SCZ_DEG_RNAseq_0.001"

gene = geneAnno1[geneAnno1$ensembl_gene_id %in% dr0.001, 2]
gene = unique(gene)
gene = gene[gene!=""]
sczmatch = gene[match(rownames(metaMat),gene)]
sczlist = ifelse(is.na(sczmatch), 0, 1)
sczlist[!(rownames(metaMat) %in% geneA)] <- NA
metaMat <- cbind(metaMat,sczlist)
colnames(metaMat)[ncol(metaMat)] <- "SCZ_downreg_RNAseq_0.001"

gene = geneAnno1[geneAnno1$ensembl_gene_id %in% ur0.001, 2]
gene = unique(gene)
gene = gene[gene!=""]
sczmatch = gene[match(rownames(metaMat),gene)]
sczlist = ifelse(is.na(sczmatch), 0, 1)
sczlist[!(rownames(metaMat) %in% geneA)] <- NA
metaMat <- cbind(metaMat,sczlist)
colnames(metaMat)[ncol(metaMat)] <- "SCZ_upreg_RNAseq_0.001"

gene = geneAnno1[geneAnno1$ensembl_gene_id %in% deg0.01, 2]
gene = unique(gene)
gene = gene[gene!=""]
sczmatch = gene[match(rownames(metaMat),gene)]
sczlist = ifelse(is.na(sczmatch), 0, 1)
sczlist[!(rownames(metaMat) %in% geneA)] <- NA
metaMat <- cbind(metaMat,sczlist)
colnames(metaMat)[ncol(metaMat)] <- "SCZ_DEG_RNAseq_0.01"

gene = geneAnno1[geneAnno1$ensembl_gene_id %in% dr0.01, 2]
gene = unique(gene)
gene = gene[gene!=""]
sczmatch = gene[match(rownames(metaMat),gene)]
sczlist = ifelse(is.na(sczmatch), 0, 1)
sczlist[!(rownames(metaMat) %in% geneA)] <- NA
metaMat <- cbind(metaMat,sczlist)
colnames(metaMat)[ncol(metaMat)] <- "SCZ_downreg_RNAseq_0.01"

gene = geneAnno1[geneAnno1$ensembl_gene_id %in% ur0.01, 2]
gene = unique(gene)
gene = gene[gene!=""]
sczmatch = gene[match(rownames(metaMat),gene)]
sczlist = ifelse(is.na(sczmatch), 0, 1)
sczlist[!(rownames(metaMat) %in% geneA)] <- NA
metaMat <- cbind(metaMat,sczlist)
colnames(metaMat)[ncol(metaMat)] <- "SCZ_upreg_RNAseq_0.01"

gene = geneAnno1[geneAnno1$ensembl_gene_id %in% deg0.05, 2]
gene = unique(gene)
gene = gene[gene!=""]
sczmatch = gene[match(rownames(metaMat),gene)]
sczlist = ifelse(is.na(sczmatch), 0, 1)
sczlist[!(rownames(metaMat) %in% geneA)] <- NA
metaMat <- cbind(metaMat,sczlist)
colnames(metaMat)[ncol(metaMat)] <- "SCZ_DEG_RNAseq_0.05"

gene = geneAnno1[geneAnno1$ensembl_gene_id %in% dr0.05, 2]
gene = unique(gene)
gene = gene[gene!=""]
sczmatch = gene[match(rownames(metaMat),gene)]
sczlist = ifelse(is.na(sczmatch), 0, 1)
sczlist[!(rownames(metaMat) %in% geneA)] <- NA
metaMat <- cbind(metaMat,sczlist)
colnames(metaMat)[ncol(metaMat)] <- "SCZ_downreg_RNAseq_0.05"

gene = geneAnno1[geneAnno1$ensembl_gene_id %in% ur0.05, 2]
gene = unique(gene)
gene = gene[gene!=""]
sczmatch = gene[match(rownames(metaMat),gene)]
sczlist = ifelse(is.na(sczmatch), 0, 1)
sczlist[!(rownames(metaMat) %in% geneA)] <- NA
metaMat <- cbind(metaMat,sczlist)
colnames(metaMat)[ncol(metaMat)] <- "SCZ_upreg_RNAseq_0.05"


############
#Insert ASD gene lists here
###########

#################
#ASDWGS-Regulatory
#################

dir = "/geschwindlabshares/atacseq/AutismWGS";
names = c("ASD.sig.all","ASD.sig.filteredTSS");

for(i in 1:length(names)){
         thislist = read.csv(paste(dir,"/ATAC.Enh.GR.",names[i],".csv",sep=""));
         DEup = which(thislist$padj<0.05 & thislist$log2FoldChange>0);
         DEdown = which(thislist$padj<0.05 & thislist$log2FoldChange<0);
         DEup = thislist[DEup,6];
         DEdown = thislist[DEdown,6];
         All= thislist[,6];
         AllDESig = c(DEup,DEdown);

         gename = geneAnno1[geneAnno1$ensembl_gene_id %in% All, 2];
         gename = gename[gename!=""];
         matchlist = match(rownames(metaMat), gename);
         matchlist = ifelse(is.na(matchlist),0,1);
         metaMat = cbind(metaMat, matchlist);
         colnames(metaMat)[dim(metaMat)[2]] = paste("ASDWGS-Regulatory-All",names[i],sep="");

         gename = geneAnno1[geneAnno1$ensembl_gene_id %in% AllDESig, 2];
         gename = gename[gename!=""];
         matchlist = match(rownames(metaMat), gename);
         matchlist = ifelse(is.na(matchlist),0,1);
         metaMat = cbind(metaMat, matchlist);
         colnames(metaMat)[dim(metaMat)[2]] = paste("ASDWGS-Regulatory-DESig",names[i],sep="");
         
         gename = geneAnno1[geneAnno1$ensembl_gene_id %in% DEup, 2];
         gename = gename[gename!=""];
         matchlist = match(rownames(metaMat), gename);
         matchlist = ifelse(is.na(matchlist),0,1);
         metaMat = cbind(metaMat, matchlist);
         colnames(metaMat)[dim(metaMat)[2]] = paste("ASDWGS-Regulatory-VZ-",names[i],sep="");

         gename = geneAnno1[geneAnno1$ensembl_gene_id %in% DEdown, 2];
         gename = gename[gename!=""];
         matchlist = match(rownames(metaMat), gename);
         matchlist = ifelse(is.na(matchlist),0,1);
         metaMat = cbind(metaMat, matchlist);
         colnames(metaMat)[dim(metaMat)[2]] = paste("ASDWGS-Regulatory-CP-",names[i],sep="");
     }


#################
#ASDWGS-All
#################

dir = "/geschwindlabshares/atacseq/AutismWGS";
names = c("ASD_all.sig.all","ASD_all.sig.filteredTSS","SIB_all.sig.all","SIB_all.sig.filteredTSS");

for(i in 1:length(names)){
         thislist = read.csv(paste(dir,"/ATAC.Enh.GR.",names[i],".noncoding.csv",sep=""));
         DEup = which(thislist$padj<0.05 & thislist$log2FoldChange>0);
         DEdown = which(thislist$padj<0.05 & thislist$log2FoldChange<0);
         DEup = thislist[DEup,6];
         DEdown = thislist[DEdown,6];
         All = thislist[,6];
         AllDESig = c(DEup,DEdown);

         gename = geneAnno1[geneAnno1$ensembl_gene_id %in% All, 2];
         gename = gename[gename!=""];
         matchlist = match(rownames(metaMat), gename);
         matchlist = ifelse(is.na(matchlist),0,1);
         metaMat = cbind(metaMat, matchlist);
         colnames(metaMat)[dim(metaMat)[2]] = paste("ASDWGS-All-All",names[i],sep="");

         gename = geneAnno1[geneAnno1$ensembl_gene_id %in% AllDESig, 2];
         gename = gename[gename!=""];
         matchlist = match(rownames(metaMat), gename);
         matchlist = ifelse(is.na(matchlist),0,1);
         metaMat = cbind(metaMat, matchlist);
         colnames(metaMat)[dim(metaMat)[2]] = paste("ASDWGS-All-DESig",names[i],sep="");

         gename = geneAnno1[geneAnno1$ensembl_gene_id %in% DEup, 2];
         gename = gename[gename!=""];
         matchlist = match(rownames(metaMat), gename);
         matchlist = ifelse(is.na(matchlist),0,1);
         metaMat = cbind(metaMat, matchlist);
         colnames(metaMat)[dim(metaMat)[2]] = paste("ASDWGS-All-VZ-",names[i],sep="");

         gename = geneAnno1[geneAnno1$ensembl_gene_id %in% DEdown, 2];
         gename = gename[gename!=""];
         matchlist = match(rownames(metaMat), gename);
         matchlist = ifelse(is.na(matchlist),0,1);
         metaMat = cbind(metaMat, matchlist);
         colnames(metaMat)[dim(metaMat)[2]] = paste("ASDWGS-All-CP-",names[i],sep="");   
     }

############
#ATAC-Hi-C Evolutionary Target Genes
############
#load(file="/geschwindlabshares/atacseq/HumanEnhancers/ATAC-Seq-HiC-Overlap/ATACSeqPCGenes2.R")
#for(i in 1:length(ATACSeqPCGenes)){
    #thislist = ATACSeqPCGenes[[i]]
    #gename = geneAnno1[geneAnno1$ensembl_gene_id %in% thislist, 2]
    #gename = gename[gename!=""]
    #matchlist = match(rownames(metaMat), gename)
    #matchlist = ifelse(is.na(matchlist),0,1)
    #metaMat = cbind(metaMat, matchlist)
    #colnames(metaMat)[dim(metaMat)[2]] = names(ATACSeqPCGenes)[i]
#}

############
#ATAC-Hi-C Evolutionary Target Genes-Lenient
############
inputdir = "/geschwindlabshares/atacseq/HumanEnhancers/ATAC-Seq-HiC-Overlap-Lenient/GOELITE/backup/input/";
files = dir("/geschwindlabshares/atacseq/HumanEnhancers/ATAC-Seq-HiC-Overlap-Lenient/GOELITE/backup/input")
for(i in 1:length(files)){
    thislist = read.csv(paste(inputdir,files[[i]],sep=""),sep="\t");
    thislist = thislist[,1];
    gename = geneAnno1[geneAnno1$ensembl_gene_id %in% thislist, 2]
    gename = gename[gename!=""]
    matchlist = match(rownames(metaMat), gename)
    matchlist = ifelse(is.na(matchlist),0,1)
    metaMat = cbind(metaMat, matchlist)
    colnames(metaMat)[dim(metaMat)[2]] = files[i]
}

##########
#Other Target Lists of interest
#########


setwd("/geschwindlabshares/atacseq/HumanEnhancers/ATAC-Seq-HiC-Overlap-Lenient/Enrichment")
save(metaMat, file="ATAC-HiC-Evo_metaMat.rda")

## remove PCDH from metaMat  #This was done by HJ because these genes are close together in the genome and HiC resolution target them all
#PCDHgrep = grep("PCDH", rownames(metaMat))
#PCDHrep = c(which(rownames(metaMat)=="PCDHA1"), which(rownames(metaMat)=="PCDH1"), which(rownames(metaMat)=="PCDHB1"), which(rownames(metaMat)=="PCDHGA1"))
#PCDHgrep = setdiff(PCDHgrep, PCDHrep)
#metaMat = metaMat[-PCDHgrep, ]

#geneLength = geneLength[-PCDHgrep]

##########
#Remove not-necessary columns from metMat
##########

all=c(1:233);
#layers = c(83:107,130:141);
#neuron_activity_column = c(19:23,1:17,47:56,61:64,60,78:82,105:121)
#binding_column = c(18,147:150)
#gene_family_column = c(65:77,83:88)
#disease_column = c(24,29,89:103,# autism
                  # 27,104, # ID
                  # 40,192:208, # schizophrenia
                  # 44:46,58:59)# Downsyndrome
#network_column = c(30:31,175:191)
target_column = c(234:243);
important_column = c(all,target_column)

metaMat = metaMat[,important_column]

### 2) Run the overlaps using gene length as a covariate
## These enrichments will be performed in a logistic regression framework

setwd("/geschwindlabshares/atacseq/HumanEnhancers/ATAC-Seq-HiC-Overlap-Lenient/Enrichment")

intGenes <- intersect(names(geneLength),names(exomeLength))
cor(geneLength[match(intGenes,names(geneLength))],
    exomeLength[match(intGenes,names(exomeLength))],
    method="spearman") ## Only 0.51 correlation (lower with Pearson) between coding and total gene length - this is interesting and good to keep in mind

exomeLength <- exomeLength[match(names(geneLength),names(exomeLength))]
dnvcols <- unique(c(grep("LGD",colnames(metaMat)), grep("denovo", colnames(metaMat)), grep("dnv", colnames(metaMat)), 
                    grep("LOF", colnames(metaMat)), grep("mis", colnames(metaMat)), grep("disrupting", colnames(metaMat)),
                    grep("sibling", colnames(metaMat)), grep("proband", colnames(metaMat))))
#genecols <- unique(c(grep("_gene", colnames(metaMat)), grep("gene", colnames(metaMat)), 
                    # grep("MAGMA", colnames(metaMat)), grep("clst", colnames(metaMat)), grep("Hi-C", colnames(metaMat))))

metaMat[metaMat>0] <- 1
table(as.numeric(data.matrix(metaMat))) ## Confirm all values are 1s and 0s

enrichOR <- enrichP <- matrix(NA,nrow=ncol(metaMat),ncol=ncol(metaMat))
colnames(enrichOR) <- rownames(enrichP) <- rownames(enrichOR) <- rownames(enrichP) <- colnames(metaMat)


for (i in 234:243) { ## Query Lists; Loop through columns of the data
  for (j in 1:233) { ## Biology/Interpretation Lists
    if (!is.na(match(j,dnvcols))) { ## if using de novo gene sets, use exome length
      thiscovar <- exomeLength
      print(paste(colnames(metaMat)[i],"vs",colnames(metaMat)[j]))
      dat1 <- as.numeric(metaMat[,i])
      dat2 <- as.numeric(metaMat[,j])
      glm.out <- glm(dat1~dat2+thiscovar,family=binomial)
      keepdat <- !is.na(dat1)&!is.na(dat2)
      
    #} else if (!is.na(match(j, genecols)) | !is.na(match(i, genecols)) & j<120) {
     # thiscovar <- geneLength ## For other gene sets, use gene length
     # print(paste(colnames(metaMat)[i],"vs",colnames(metaMat)[j]))
     # dat1 <- as.numeric(metaMat[,i])
     # dat2 <- as.numeric(metaMat[,j])
     # glm.out <- glm(dat1~dat2+thiscovar,family=binomial)
     # keepdat <- !is.na(dat1)&!is.na(dat2)
      
    } else{
      print(paste(colnames(metaMat)[i],"vs",colnames(metaMat)[j]))
      dat1 <- as.numeric(metaMat[,i])
      dat2 <- as.numeric(metaMat[,j])
      glm.out <- glm(dat1~dat2,family=binomial)
      keepdat <- !is.na(dat1)&!is.na(dat2)
    }
    ## Compute glm that asks how well dat2 predicts dat1 when controlling for thiscovar
    
    enrichP[i,j] <- summary(glm.out)$coefficients[2,4]
    enrichOR[i,j] <- summary(glm.out)$coefficients[2,1]
  }
}

## Save the output in case you don't want to run all of the above again
save(enrichP,enrichOR,metaMat, file="AllEnrichments-HumanEvo-noncoding.Rdata")

Pmat <- enrichP
diag(Pmat) <- 0
Pmat[upper.tri(Pmat)] <- t(Pmat)[upper.tri(Pmat)] ## Fill in the other half, keep a symmetric matrix

ORmat <- enrichOR
diag(ORmat) <- Inf
ORmat[upper.tri(Pmat)] <- t(ORmat)[upper.tri(ORmat)] ## Fill in the other half, keep a symmetric matrix

#moduleSizes <- apply(metaMat,2,sum,na.rm=TRUE)

### 3) Take a subset of the relationships for plotting
ASD = c(114:129,172:216); 
SCZ = c(217:233);
Layers = c(19:23,83:107,130:171);

#keepcols <- c(234:243) ###################################### Subset Target Lists here
keepcols2 = c(1:233);        #Interpretation Lists
#keepcols1 = c(238:239,246:247,254:255)  #Target Lists
keepcols1 = c(234:243);

subORmat <- ORmat[keepcols1,keepcols2]
subPmat <- Pmat[keepcols1,keepcols2]

dispMat <- log2(exp(subORmat)) ## Tranform to log2 odds ratio
rmval <- subPmat > 0.05 ## Display only values with p > 0.05
dispMat[rmval] <- 0
rownames(dispMat) <- rownames(subORmat)

#Use this with single query list vs all bio lists

#load("/geschwindlabshares/atacseq/HumanEnhancers/ATAC-Seq-HiC-Overlap/Enrichment/AllEnrichments-ATAC-HiC-Evo.Rdata")
#setwd("/geschwindlabshares/atacseq/HumanEnhancers/ATAC-Seq-HiC-Overlap/Enrichment")

#pdf("AllEnrichments-ATAC-HiC-Evo.pdf",height=5,width=20)
#FDRmat = subPmat; FDRfilt <- FDRmat<0.05

#textMat <- paste(signif(dispMat,2),"\n(",signif(FDRmat,2),")",sep="")
#textMat[!FDRfilt|subORmat==Inf] <- ""

#dispMat = matrix(dispMat, ncol=length(dispMat), nrow=1) #This is the color gradient for OR
#textMat = matrix(textMat, ncol=length(textMat), nrow=1) #This is the actual OR and p-values
#colnames(dispMat)=colnames(textMat)=colnames(metaMat)[-dim(metaMat)[2]]
#rownames(dispMat)=rownames(textMat)="ASDWGS"

#labeledHeatmap(Matrix=dispMat,
 #              textMatrix=textMat,
  #             xLabels=colnames(dispMat),
   #            yLabels=rownames(dispMat),
    #           colors=blueWhiteRed(100),
     #          cex.lab.x=0.4,
      #         cex.lab.y=0.4,
       #        zlim=c(-2,2),
        #       cex.text=0.2,
         #      setStdMargins=FALSE)
#dev.off()

#Use this with multiple query lists vs all bio lists

#load("/geschwindlabshares/atacseq/HumanEnhancers/ATAC-Seq-HiC-Overlap/Enrichment/AllEnrichments-ASDWGS-noncoding.Rdata")
setwd("/geschwindlabshares/atacseq/HumanEnhancers/ATAC-Seq-HiC-Overlap-Lenient/Enrichment")

pdf("AllEnrichments-HumanEvo.pdf",height=5,width=20)
FDRmat = subPmat; FDRfilt <- FDRmat<0.05

textMat <- matrix(paste(signif(dispMat,2),"\n(",signif(FDRmat,2),")",sep=""),nrow=nrow(dispMat),ncol=ncol(dispMat))
textMat[!FDRfilt|subORmat==Inf] <- ""

labeledHeatmap(Matrix=dispMat,
               textMatrix=textMat,
               xLabels=colnames(dispMat),
               yLabels=rownames(dispMat),
               colors=blueWhiteRed(100),
               cex.lab.x=0.4,
               cex.lab.y=0.4,
               zlim=c(-2,2),
               cex.text=0.4,
               setStdMargins=TRUE)
dev.off();
