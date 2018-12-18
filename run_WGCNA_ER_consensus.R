##################Consensus Network WGCNA####################
# Consensus performed between Downs syndrome microarray data in cortex regions (CBC, DFC, FC, IPC, ITC, MFC, OFC, S1C, STC, V1C) from GSE59630,
# Alzheimer's RNA seq data in temporal cortex from Mayo Clinic, and Alzheimer's microarray data in frontal cortex from Zhang 
setwd("~/WGCNA")
getwd()
library(WGCNA)
library(flashClust)
library(gplots)
library(cluster) 
library(igraph); #for part 8
library(RColorBrewer); #for part 8
options(stringsAsFactors=FALSE)
#Only when not on server, comment out:
enableWGCNAThreads()
#=====================================================================================
#
#  Part 1: Data Input, cleaning, and pre-processing
#
#=====================================================================================
# Data from GSE59630 (downs syndrome microarray) prepared with GEO2R - everything until data cleanup generated online


# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Tue Oct 30 13:30:48 EDT 2018
################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO
gset <- getGEO("GSE59630", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL5175", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("00000000000000000000000000000000000000000000000000",
               "00000000111111111111111111111111111111111111111111",
               "1111111111111111")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
#exprs(gset) is expression matrix
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","GB_LIST","SPOT_ID","RANGE_GB","RANGE_STRAND","RANGE_START"))
write.table(tT, file="DE_GEO.txt", row.names=F, sep="\t")

##########################Clean up data######################################
targets.ALL = pData(gset);
datExpr = as.data.frame(exprs(gset));

#remove hippocampus regions - 6 total so 110 patients left
targets.ALL.select  = targets.ALL[!(targets.ALL$region == "HIP"),]
#remove extra columns
targets.ALL.select = as.data.frame(targets.ALL.select)
targets.ALL.select2 = targets.ALL.select[ , which(names(targets.ALL.select) %in% c("geo_accession","disease status:ch1", 
                                                              "age:ch1", "region:ch1", "Sex:ch1"))]
#remove hippocampus regions from expression
targets.ALL.HIP  = targets.ALL[(targets.ALL$region == "HIP"),]
GSM = targets.ALL.HIP$geo_accession;
datExpr.select = datExpr[ , -which(names(datExpr) %in% c(GSM))]

#gene IDs are affymetrix, convert to ensembl:
annot=read.csv("/home/erebboah/WGCNA/ensembl_GEO.csv")
expr = as.integer(rownames(datExpr.select));
ensembl = annot$From;
ind=intersect(ensembl, expr)
#16597 x 110
datExpr.select.subset=datExpr.select[na.omit(match(ind,rownames(datExpr.select))),]
annot1= annot[na.omit(match(ind,annot$From)),]
#16077 x 110
datExpr_ENSGID=collapseRows(datExpr.select.subset, annot1$To, as.character(annot1$From))$datETcollapsed  # changing to ENSG_IDs using collapseRows
save(datExpr_ENSGID, targets.ALL.select2, file="datExpr_GSE59630.rda")

#########################################################################
#AD data input (temporal cortex)
load("/home/vivek/AMP_AD/Mayo/Analysis/Step04_WGCNA/Discovery_Set/AD/rWGCNA_Mayo_ForPreservation.rda") #Alzheimer's data
#make expression and target variables for AD
datExpr.Ref.AD_1 = datExpr.Ref
targets.Ref.AD_1 = targets.Ref

#AD data input #2 (frontal cortex)
load("/home/vivek/AD/Zhang/Zhang_Regressed.rda")
datExpr.Ref.AD_2 = t(normExpr.reg);
targets.Ref.AD_2 = targets; 

#DS data input
load("datExpr_GSE59630.rda") 
#Make sure genes in columns and samples in rows, make variables for expression and target variables for DS
datExpr.Ref.DS=t(datExpr_ENSGID);
targets.Ref.DS=targets.ALL.select2;

#Use only intersecting genes - 11810 tot intersect. DS has 110 samples, AD_1 has 104, AD_2 has 465
gnS=intersect(colnames(datExpr.Ref.DS),intersect(colnames(datExpr.Ref.AD_1),colnames(datExpr.Ref.AD_2)))
datExpr.Ref.DS=datExpr.Ref.DS[,match(gnS,colnames(datExpr.Ref.DS))]
datExpr.Ref.AD_1 =datExpr.Ref.AD_1 [,match(gnS,colnames(datExpr.Ref.AD_1))]
datExpr.Ref.AD_2 =datExpr.Ref.AD_2 [,match(gnS,colnames(datExpr.Ref.AD_2))]

#Make multiExpr list of 3 data frames
nSets=3;
setLabels=c("Downs syndrome samples", "Alzheimer's samples 1", "Alzheimer's samples 2")
#Form multi-set expression data
multiExpr = vector(mode="list", length=nSets)

#Remove rownames and column names from datExpr dfs, re-assign in multiExpr 
genes = colnames(datExpr.Ref.DS);
samples = rownames(datExpr.Ref.DS);
rownames(datExpr.Ref.DS)=NULL;
colnames(datExpr.Ref.DS)=NULL;
multiExpr[[1]] = list(data = datExpr.Ref.DS);
colnames(multiExpr[[1]]$data) = genes;
rownames(multiExpr[[1]]$data) = samples;

genes = colnames(datExpr.Ref.AD_1);
samples = rownames(datExpr.Ref.AD_1);
rownames(datExpr.Ref.AD_1)=NULL;
colnames(datExpr.Ref.AD_1)=NULL;
multiExpr[[2]] = list(data = datExpr.Ref.AD_1);
names(multiExpr[[2]]$data) = genes;
rownames(multiExpr[[2]]$data) = samples;

genes = colnames(datExpr.Ref.AD_2);
samples = rownames(datExpr.Ref.AD_2);
rownames(datExpr.Ref.AD_2)=NULL;
colnames(datExpr.Ref.AD_2)=NULL;
multiExpr[[3]] = list(data = datExpr.Ref.AD_2);
colnames(multiExpr[[3]]$data) = genes;
rownames(multiExpr[[3]]$data) = samples;

exprSize = checkSets(multiExpr)
exprSize #Check data is in correct format; checkSets function checks whether given sets have the correct format and retrieves dimensions.

#$nGenes
#[1] 11810
#$nSamples
#[1] 110 104 465

#Make multiMeta of traits for all sets - list of lists
multiMeta=list(aging = list(data=targets.Ref.DS), AD1 = list(data=targets.Ref.AD_1), AD2 = list(data=targets.Ref.AD_2))

#Check that all genes and samples have sufficiently low numbers of missing values
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK

#If NOT all ok - 
if (!gsg$allOK)
{
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], 
                                              collapse = ", ")))
  for (set in 1:exprSize$nSets)
  {
    if (sum(!gsg$goodSamples[[set]]))
      printFlush(paste("In set", setLabels[set], "removing samples",
                       paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
  }
  # Update exprSize
  exprSize = checkSets(multiExpr)
}

#Save everything
nGenes=exprSize$nGenes;
nSamples=exprSize$nSamples;
save(multiExpr, multiMeta, nGenes, nSamples, setLabels, 
     targets.Ref.AD_1, targets.Ref.AD_2, targets.Ref.DS, 
     datExpr.Ref.AD_1, datExpr.Ref.AD_2, datExpr.Ref.DS, file = "Consensus_dataInput_112418.rda");
#=====================================================================================
#
#  Part 2: Choose soft thresholding power
#
#=====================================================================================
# Load the data saved in the first part
load(file = "Consensus_dataInput_112418.rda");

# Get the number of sets in the multiExpr structure.
nSets = checkSets(multiExpr)$nSets

# Choose a set of soft-thresholding powers
powers = c(seq(1,10,by=1), seq(12,30, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 5, networkType="signed", corFnc="bicor")[[2]]);

# Plot the results
pdf("1_Power_DSAD_112418.pdf", height=10, width=18)
colors = c("blue", "red", "black")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity");

# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}
dev.off()

#=====================================================================================
#
#  Part 3: Network Construction
#
#=====================================================================================
#Set soft thresholding power to number based on plots
softpower=16;
# Auto network
load(file = "Consensus_dataInput_112418.rda");
net=blockwiseConsensusModules(multiExpr, blocks = NULL,
                              maxBlockSize = 30000, ## This should be set to a smaller size if the user has limited RAM
                              randomSeed = 12345,
                              corType = "pearson", ## no use for bicor
                              power = softpower,
                              consensusQuantile = 0.2,
                              networkType = "signed",
                              TOMType = "unsigned",
                              TOMDenom = "min",
                              scaleTOMs = TRUE, scaleQuantile = 0,
                              sampleForScaling = TRUE, sampleForScalingFactor = 1000,
                              useDiskCache = TRUE, chunkSize = NULL,
                              deepSplit = 2,
                              detectCutHeight = 0.99999999, minModuleSize = 20,
                              mergeCutHeight = 0.2,
                              saveConsensusTOMs = TRUE,
                              consensusTOMFilePattern = "ConsensusTOM-block.%b.rda")
save(list=ls(),file="consensus_112418.rda")

consMEs = net$multiMEs;
moduleColors = net$colors;
table(moduleColors) 

load("ConsensusTOM-block.1.rda") # consensus TOM
consTree= hclust(1-consTomDS,method="average");
save(list=ls(),file="consensus_112418.rda")

# Relate modules with traits for each dataset
##################################AD 1################################################
load("Consensus_dataInput_112418.rda");

Diagnosis=factor(as.character(targets.Ref.AD_1$Diagnosis))
Diagnosis<-relevel(Diagnosis,'Control')
Diagnosis=as.numeric(Diagnosis)
Age=as.numeric(targets.Ref.AD_1$AgeAtDeath)
RIN=as.numeric(targets.Ref.AD_1$RIN)
PC1_Sequencing=as.numeric(targets.Ref.AD_1$Seq.PC1)
PC2_Sequencing=as.numeric(targets.Ref.AD_1$Seq.PC2)
Brain.Bank=as.numeric(factor(targets.Ref.AD_1$Source))
Gender=as.numeric(factor(targets.Ref.AD_1$Gender))

geneSigsAD1=matrix(NA,nrow=7,ncol=ncol(datExpr.Ref.AD_1)) 
for(i in 1:ncol(geneSigsAD1)) {
  exprvec=as.numeric(datExpr.Ref.AD_1[,i])
  ager=bicor(Age,exprvec,use="pairwise.complete.obs")
  sexr=bicor(exprvec, Gender,use="pairwise.complete.obs")
  sourcer=bicor(exprvec, Brain.Bank,use="pairwise.complete.obs")
  conditionr=bicor(exprvec, Diagnosis,use="pairwise.complete.obs")
  rinr=bicor(RIN,exprvec,use="pairwise.complete.obs")
  pc1r=bicor(PC1_Sequencing,exprvec,use="pairwise.complete.obs")
  pc2r=bicor(PC2_Sequencing,exprvec,use="pairwise.complete.obs")
  geneSigsAD1[,i]=c(ager, sexr, sourcer,conditionr,rinr,pc1r,pc2r)
  cat('Done for gene...',i,'\n')
}

geneSigsAD1[1,] =numbers2colors(as.numeric(geneSigsAD1[1,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigsAD1[2,] =numbers2colors(as.numeric(geneSigsAD1[2,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigsAD1[3,] =numbers2colors(as.numeric(geneSigsAD1[3,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigsAD1[4,] =numbers2colors(as.numeric(geneSigsAD1[4,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigsAD1[5,] =numbers2colors(as.numeric(geneSigsAD1[5,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigsAD1[6,] =numbers2colors(as.numeric(geneSigsAD1[6,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigsAD1[7,] =numbers2colors(as.numeric(geneSigsAD1[7,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))

rownames(geneSigsAD1)=c("Age","Gender","Brain.Bank","Diagnosis","RIN","PC1_Sequencing","PC2_Sequencing")

##################################AD 2################################################
Diagnosis=factor(as.character(targets.Ref.AD_2$Disease))
Diagnosis<-relevel(Diagnosis,'Control')
Diagnosis=as.numeric(Diagnosis)
Age=as.numeric(targets.Ref.AD_2$Age)
GEO=as.numeric(factor(targets.Ref.AD_2$geo_accession))
Gender=as.numeric(factor(targets.Ref.AD_2$Sex))

geneSigsAD2=matrix(NA,nrow=4,ncol=ncol(datExpr.Ref.AD_2)) 
for(i in 1:ncol(geneSigsAD2)) {
  exprvec=as.numeric(datExpr.Ref.AD_2[,i])
  ager=bicor(Age,exprvec,use="pairwise.complete.obs")
  sexr=bicor(exprvec, Gender,use="pairwise.complete.obs")
  conditionr=bicor(exprvec, Diagnosis,use="pairwise.complete.obs")
  geor=bicor(GEO,exprvec,use="pairwise.complete.obs")
  geneSigsAD2[,i]=c(ager, sexr,conditionr,geor)
  cat('Done for gene...',i,'\n')
}

geneSigsAD2[1,] =numbers2colors(as.numeric(geneSigsAD2[1,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigsAD2[2,] =numbers2colors(as.numeric(geneSigsAD2[2,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigsAD2[3,] =numbers2colors(as.numeric(geneSigsAD2[3,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigsAD2[4,] =numbers2colors(as.numeric(geneSigsAD2[4,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))

rownames(geneSigsAD2)=c("Age","Gender","Diagnosis", "GEO")

###############################DS########################################
targets.Ref.DS$Age = c("0", "0", "0", "0", "0.333", "0.333", "0.333", "0.5", "0.5", "0.5",
                       "0.5", "0.5", "0.5", "0.833", "0.833", "0.833", "0.833", "0.833", "0.833", "0.833", 
                       "0.833", "0.833", "0.833", "1", "1", "2", "2", "3", "3", "3",
                       "3", "3", "8", "8", "8", "15", "15", "15", "15", "18", 
                       "18", "18", "22", "22", "22", "30", "30", "30", "30", "30",
                       "30", "42", "42", "42", "42", "0", "0", "0", "0", "0.083", 
                       "0.083", "0.083", "0.5", "0.5", "0.5", "0.5", "0.5", "0.5", "0.75", "0.75",
                       "0.75", "0.75", "0.75", "0.75", "0.75", "0.75", "0.75", "0.75", "1.17", "1.17", 
                       "3", "3", "3", "3", "3", "2", "2", "10", "10", "10",
                       "13", "13", "13", "13", "19", "19", "19", "22", "22", "22",
                       "39", "39", "39", "39", "39", "39", "40", "40", "40", "40")

Age=as.numeric(targets.Ref.DS$Age)
Gender=as.numeric(factor(targets.Ref.DS$Sex))
Diagnosis = as.numeric(factor(targets.Ref.DS$disease))
Region = as.numeric(factor(targets.Ref.DS$region))

geneSigsDS=matrix(NA,nrow=4,ncol=ncol(datExpr.Ref.DS)) 
for(i in 1:ncol(geneSigsAD2)) {
  exprvec=as.numeric(datExpr.Ref.DS[,i])
  ager=bicor(Age,exprvec,use="pairwise.complete.obs")
  sexr=bicor(exprvec, Gender,use="pairwise.complete.obs")
  conditionr=bicor(exprvec, Diagnosis,use="pairwise.complete.obs")
  regionr=bicor(Region,exprvec,use="pairwise.complete.obs")
  
  geneSigsDS[,i]=c(ager, sexr,conditionr,regionr)
  cat('Done for gene...',i,'\n')
}

geneSigsDS[1,] =numbers2colors(as.numeric(geneSigsDS[1,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigsDS[2,] =numbers2colors(as.numeric(geneSigsDS[2,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigsDS[3,] =numbers2colors(as.numeric(geneSigsDS[3,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
geneSigsDS[4,] =numbers2colors(as.numeric(geneSigsDS[4,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))

rownames(geneSigsDS)=c("Age","Gender","Diagnosis", "Region")

# Calculate modules for each set of parameters
load("consensus_112418.rda")
mColorh <- mLabelh <- colorLabels <- NULL
for (minModSize in c(40,100,160)) {
  for (dthresh in c(0.1,0.2,0.25)) {
    for (ds in c(2,4)) {
      print("Trying parameters:")
      print(c(minModSize,dthresh,ds))
      tree = cutreeHybrid(dendro = consTree, pamStage=FALSE,
                          minClusterSize = minModSize, cutHeight = 0.99999999,
                          deepSplit = ds, distM = as.matrix(1-consTomDS))
      
      merged <- mergeCloseModules(exprData = multiExpr,colors = tree$labels,
                                  cutHeight = dthresh)
      mColorh <- cbind(mColorh,labels2colors(merged$colors))
      mLabelh <- c(mLabelh,paste("DS=",ds," mms=\n",minModSize," dcor=",dthresh))
    }
  }
}

# Plotting modules for each set of params and traits 
mColorh1=cbind(mColorh,geneSigsAD1[1,],geneSigsAD1[2,],geneSigsAD1[3,],geneSigsAD1[4,],geneSigsAD1[5,],geneSigsAD1[6,],geneSigsAD1[7,],
               geneSigsAD2[1,], geneSigsAD2[2,], geneSigsAD2[3,], geneSigsAD2[4,],
               geneSigsDS[1,], geneSigsDS[2,], geneSigsDS[3,], geneSigsDS[4,])
rownames_geneSigs = c(rownames(geneSigsAD1), rownames(geneSigsAD2), rownames(geneSigsDS))
mLabelh1=c(mLabelh,rownames_geneSigs)

pdf("ConsensusTOM_MultiDendro_DSAD_112418.pdf",height=25,width=20)
plotDendroAndColors(consTree, mColorh1, groupLabels = mLabelh1,addGuide=TRUE,dendroLabels=FALSE,main= paste("Signed consensus network with power = 16"));
dev.off()

##############Choose parameters for final dendrogram#########################
load('/home/erebboah/WGCNA/ConsensusTOM-block.1.rda') #results of network analysis

mms=160
ds=4
dthresh=0.1

tree = cutreeHybrid(dendro = consTree, pamStage=FALSE,
                    minClusterSize = mms, cutHeight = 0.99999999,
                    deepSplit = ds, distM = as.matrix(1-consTomDS))

merged <- mergeCloseModules(exprData = multiExpr,colors = tree$labels,
                            cutHeight = dthresh)

# Eigengenes of the new merged modules:
MEs_DS = merged$newMEs[[1]]$data;
MEs_AD1 = merged$newMEs[[2]]$data;
MEs_AD2 = merged$newMEs[[3]]$data;

# Module color associated with each gene
moduleColor.cons <- labels2colors(merged$colors)

mColorh <- cbind(labels2colors(merged$colors))
mLabelh <- c("Merged Colors")

mColorh1=cbind(mColorh,geneSigsAD1[1,],geneSigsAD1[2,],geneSigsAD1[3,],geneSigsAD1[4,],geneSigsAD1[5,],geneSigsAD1[6,],geneSigsAD1[7,],
               geneSigsAD2[1,], geneSigsAD2[2,], geneSigsAD2[3,], geneSigsAD2[4,],
               geneSigsDS[1,], geneSigsDS[2,], geneSigsDS[3,], geneSigsDS[4,])
rownames_geneSigs = c(rownames(geneSigsAD1), rownames(geneSigsAD2), rownames(geneSigsDS))
mLabelh1=c(mLabelh,rownames_geneSigs)

pdf("ConsensusTOM_FinalDendro_DSAD_112718.pdf",height=10,width=16)
plotDendroAndColors(consTree, mColorh1, groupLabels = mLabelh1,addGuide=TRUE,dendroLabels=FALSE,main= paste("Signed bicor network with power = 16, mms=",mms,"ds=",ds,"dthresh=",dthresh,"cquant=0.2"));
dev.off()

#Just look at diagnosis traits
mColorh1=cbind(mColorh,geneSigsAD1[4,], geneSigsAD2[3,], geneSigsDS[3,])
rownames_geneSigs = c("AD1 Diagnosis", "AD2 Diagnosis", "DS Diagnosis")
mLabelh1=c(mLabelh,rownames_geneSigs)

pdf("ConsensusTOM_FinalDendro_DSAD_diag_112718.pdf",height=10,width=16)
plotDendroAndColors(consTree, mColorh1, groupLabels = mLabelh1,addGuide=TRUE,dendroLabels=FALSE,main= paste("Signed bicor network with power = 16, mms=",mms,"ds=",ds,"dthresh=",dthresh,"cquant=0.2"));
dev.off()

# Convert numerical lables to colors for labeling of modules 
MEColors=labels2colors(as.numeric(substring(names(MEs_DS), 3)));
MEColorNames = paste("ME", MEColors, sep="");

colnames(MEs_DS)=c(MEColorNames)
colnames(MEs_AD1)=c(MEColorNames)
colnames(MEs_AD2)=c(MEColorNames)

#Calculate kMEs and p value, Z score of kMEs for each gene (connectivity of each gene to a module eigengene) across datasets
consKME1=consensusKME(multiExpr=multiExpr, moduleColor.cons,
                      multiEigengenes = NULL,
                      consensusQuantile = 0.2,
                      signed = TRUE)
#Consensus kMEs:
consensus.KMEs=consKME1[,regexpr('consensus.kME',names(consKME1))>0]

save(consensus.KMEs, consKME1, multiExpr, MEColorNames, moduleColor.cons, consTree, 
     MEs_DS, MEs_AD1, MEs_AD2, file="Consensus_MEs_112418.rda")

#=====================================================================================
#
#  Part 4: Get annotated gene lists and their associated kMEs
#
#=====================================================================================
load("Consensus_MEs_112418.rda")
ensembl=read.csv("/home/vivek/FTD_Seeley/Analysis_Nov2017/ENSG85_Human.csv.gz") # Convert Ensembl gene ID to gene names
consensus.KMEs$Ensembl.Gene.ID=paste(rownames(consensus.KMEs))

merged=merge(consensus.KMEs,ensembl,by.x="Ensembl.Gene.ID",by.y="Ensembl.Gene.ID",all.x=T)
ind=match(consensus.KMEs$Ensembl.Gene.ID,merged$Ensembl.Gene.ID)
merged1=merged[ind,]
consensus.KMEs.annot=merged1

geneInfo.cons=as.data.frame(cbind(consensus.KMEs.annot$Ensembl.Gene.ID,consensus.KMEs.annot$Associated.Gene.Name,
                                  moduleColor.cons,consensus.KMEs))
geneInfo.cons=geneInfo.cons[,-ncol(geneInfo.cons)] # check if last column is Ensembl gene id

colnames(geneInfo.cons)[1]= "Ensembl.Gene.ID"
colnames(geneInfo.cons)[2]= "GeneSymbol"
colnames(geneInfo.cons)[3]= "Initially.Assigned.Module.Color"

write.csv(geneInfo.cons,'geneInfo.cons.DSAD_112418.csv') #Final annotated geneInfo file is input to the rest of the analysis
save(list=ls(),file='geneInfo.cons.112418.rda')

#=====================================================================================
#
#  Part 5: Module-trait relationships
#
#=====================================================================================
load("geneInfo.cons.112418.rda")

##################################AD 1################################################
nSamples = nrow(datExpr.Ref.AD_1);
nGenes = ncol(datExpr.Ref.AD_1);

Diagnosis=factor(as.character(targets.Ref.AD_1$Diagnosis))
Diagnosis<-relevel(Diagnosis,'Control')
Diagnosis=as.numeric(Diagnosis)
Age=as.numeric(targets.Ref.AD_1$AgeAtDeath)
RIN=as.numeric(targets.Ref.AD_1$RIN)
PC1_Sequencing=as.numeric(targets.Ref.AD_1$Seq.PC1)
PC2_Sequencing=as.numeric(targets.Ref.AD_1$Seq.PC2)
Brain.Bank=as.numeric(factor(targets.Ref.AD_1$Source))
Gender=as.numeric(factor(targets.Ref.AD_1$Gender))

factors1_AD1=cbind(Diagnosis, Age, Gender, Brain.Bank, RIN, PC1_Sequencing, PC2_Sequencing)

PCvalues<-MEs_AD1[,-ncol(MEs_AD1)] #exclude grey
dim(PCvalues)

moduleTraitCor_AD1 = cor(PCvalues, factors1_AD1, use = "p");
moduleTraitPvalue_AD1 = corPvalueStudent(moduleTraitCor_AD1, nSamples);
colnames(moduleTraitPvalue_AD1) = paste("p.value.", colnames(moduleTraitCor_AD1), sep="");

## Use the text function with the FDR filter in labeledHeatmap to add asterisks, e.g.
txtMat <-  signif(moduleTraitPvalue_AD1,2)
txtMat[txtMat>=0.05] <- ""
txtMat[txtMat <0.05&txtMat >0.01] <- "*"
txtMat[txtMat <0.01&txtMat >0.005] <- "**"
#txtMat[txtMat <0.005&txtMat >0] <- "***"

txtMat1 <- signif( moduleTraitCor_AD1,2)
#we only want to look at pearson correlations in certain range
txtMat1[txtMat1> -0.3&txtMat1<0.2] <- "" 
textMatrix1 = paste( txtMat1, '\n', '(',txtMat ,')', sep = '');
textMatrix1= matrix(textMatrix1,ncol=ncol( moduleTraitPvalue_AD1),nrow=nrow(moduleTraitPvalue_AD1))

#Plot heatmap
pdf(paste('NetworkPlot_AD1_consensus_112418.pdf'),width=16,height=30)
par( mar = c(8, 12, 3, 3) );
labeledHeatmap( Matrix = moduleTraitCor_AD1,
                xLabels = colnames(factors1_AD1),
                yLabels = rownames(moduleTraitPvalue_AD1),
                ySymbols = rownames(moduleTraitPvalue_AD1),
                colorLabels = FALSE,
                colors = blueWhiteRed(50),
                textMatrix = textMatrix1,
                setStdMargins = FALSE,
                cex.text = 1.5,
                zlim = c(-1, 1),
                cex.lab.x = 1.2,
                main = paste("Module-trait relationships")
);

#Plot eigengene heatmap
par(cex = 1.0)
plotEigengeneNetworks(MEs_AD1, "Eigengene Network", marHeatmap = c(3,4,2,2), 
                      marDendro = c(0,4,1,2),cex.adjacency = 0.3,plotDendrograms = TRUE, 
                      xLabelsAngle = 90,heatmapColors=blueWhiteRed(100)[51:100])

#Plot boxplots, scatterplots
toplot=t(MEs_AD1)
cols=substring(colnames(MEs_AD1),3,20)
par(mfrow=c(4,4))
par(mar=c(5,6,4,2))
for (i in 1:nrow(toplot)) {
  boxplot(toplot[i,]~factor(as.vector(as.factor(targets.Ref.AD_1$Diagnosis)),c('Control','AD')),col=cols[i],ylab="ME",main=rownames(toplot)[i],xlab=NULL,las=2)
  verboseScatterplot(x=as.numeric(targets.Ref.AD_1$Age),y=toplot[i,],xlab="Age",ylab="ME",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=cols[i],pch=19)
  boxplot(toplot[i,]~factor(targets.Ref.AD_1$Gender),col=cols[i],ylab="ME",main=rownames(toplot)[i],xlab=NULL,las=2)
}

dev.off()

##################################AD 2################################################
nSamples = nrow(datExpr.Ref.AD_2);
nGenes = ncol(datExpr.Ref.AD_2);

Diagnosis=factor(as.character(targets.Ref.AD_2$Disease))
Diagnosis<-relevel(Diagnosis,'Control')
Diagnosis=as.numeric(Diagnosis)
Age=as.numeric(targets.Ref.AD_2$Age)
GEO=as.numeric(factor(targets.Ref.AD_2$geo_accession))
Gender=as.numeric(factor(targets.Ref.AD_2$Sex))

factors1_AD2=cbind(Diagnosis, Age, Gender, GEO)

PCvalues<-MEs_AD2[,-ncol(MEs_AD2)] #exclude grey
dim(PCvalues)

moduleTraitCor_AD2 = cor(PCvalues, factors1_AD2, use = "p");
moduleTraitPvalue_AD2 = corPvalueStudent(moduleTraitCor_AD2, nSamples);
colnames(moduleTraitPvalue_AD2) = paste("p.value.", colnames(moduleTraitCor_AD2), sep="");

## Use the text function with the FDR filter in labeledHeatmap to add asterisks, e.g.
txtMat <-  signif(moduleTraitPvalue_AD2,2)
txtMat[txtMat>=0.05] <- ""
txtMat[txtMat <0.05&txtMat >0.01] <- "*"
txtMat[txtMat <0.01&txtMat >0.005] <- "**"
#txtMat[txtMat <0.005&txtMat >0] <- "***"

txtMat1 <- signif( moduleTraitCor_AD2,2)
#we only want to look at pearson correlations in certain range
txtMat1[txtMat1> -0.3&txtMat1<0.2] <- "" 
textMatrix1 = paste( txtMat1, '\n', '(',txtMat ,')', sep = '');
textMatrix1= matrix(textMatrix1,ncol=ncol( moduleTraitPvalue_AD2),nrow=nrow(moduleTraitPvalue_AD2))

#Plot heatmap
pdf(paste('NetworkPlot_AD2_consensus_112418.pdf'),width=16,height=30)
par( mar = c(8, 12, 3, 3) );
labeledHeatmap( Matrix = moduleTraitCor_AD2,
                xLabels = colnames(factors1_AD2),
                yLabels = rownames(moduleTraitPvalue_AD2),
                ySymbols = rownames(moduleTraitPvalue_AD2),
                colorLabels = FALSE,
                colors = blueWhiteRed(50),
                textMatrix = textMatrix1,
                setStdMargins = FALSE,
                cex.text = 1.5,
                zlim = c(-1, 1),
                cex.lab.x = 1.2,
                main = paste("Module-trait relationships")
);

#Plot eigengene heatmap
par(cex = 1.0)
plotEigengeneNetworks(MEs_AD2, "Eigengene Network", marHeatmap = c(3,4,2,2), 
                      marDendro = c(0,4,1,2),cex.adjacency = 0.3,plotDendrograms = TRUE, 
                      xLabelsAngle = 90,heatmapColors=blueWhiteRed(100)[51:100])

#Plot boxplots, scatterplots
toplot=t(MEs_AD2)
cols=substring(colnames(MEs_AD2),3,20)
par(mfrow=c(4,4))
par(mar=c(5,6,4,2))
for (i in 1:nrow(toplot)) {
  boxplot(toplot[i,]~factor(as.vector(as.factor(targets.Ref.AD_2$Disease)),c('Control','AD')),col=cols[i],ylab="ME",main=rownames(toplot)[i],xlab=NULL,las=2)
  verboseScatterplot(x=as.numeric(targets.Ref.AD_2$Age),y=toplot[i,],xlab="Age",ylab="ME",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=cols[i],pch=19)
  boxplot(toplot[i,]~factor(targets.Ref.AD_2$Sex),col=cols[i],ylab="ME",main=rownames(toplot)[i],xlab=NULL,las=2)
}

dev.off()

###############################DS########################################
nSamples = nrow(datExpr.Ref.DS);
nGenes = ncol(datExpr.Ref.DS);

targets.Ref.DS$Age = c("0", "0", "0", "0", "0.333", "0.333", "0.333", "0.5", "0.5", "0.5",
                       "0.5", "0.5", "0.5", "0.833", "0.833", "0.833", "0.833", "0.833", "0.833", "0.833", 
                       "0.833", "0.833", "0.833", "1", "1", "2", "2", "3", "3", "3",
                       "3", "3", "8", "8", "8", "15", "15", "15", "15", "18", 
                       "18", "18", "22", "22", "22", "30", "30", "30", "30", "30",
                       "30", "42", "42", "42", "42", "0", "0", "0", "0", "0.083", 
                       "0.083", "0.083", "0.5", "0.5", "0.5", "0.5", "0.5", "0.5", "0.75", "0.75",
                       "0.75", "0.75", "0.75", "0.75", "0.75", "0.75", "0.75", "0.75", "1.17", "1.17", 
                       "3", "3", "3", "3", "3", "2", "2", "10", "10", "10",
                       "13", "13", "13", "13", "19", "19", "19", "22", "22", "22",
                       "39", "39", "39", "39", "39", "39", "40", "40", "40", "40")

Age=as.numeric(targets.Ref.DS$Age)
Gender=as.numeric(factor(targets.Ref.DS$Sex))
Diagnosis = as.numeric(factor(targets.Ref.DS$disease))
Region = as.numeric(factor(targets.Ref.DS$region))

factors1_DS=cbind(Diagnosis, Age, Gender, Region)

PCvalues<-MEs_DS[,-ncol(MEs_DS)] #exclude grey
dim(PCvalues) #number of modules

moduleTraitCor_DS= cor(PCvalues, factors1_DS, use = "p");
moduleTraitPvalue_DS = corPvalueStudent(moduleTraitCor_DS, nSamples);
colnames(moduleTraitPvalue_DS) = paste("p.value.", colnames(moduleTraitCor_DS), sep="");

## Use the text function with the FDR filter in labeledHeatmap to add asterisks, e.g.
txtMat <-  signif(moduleTraitPvalue_DS,2)
txtMat[txtMat>=0.05] <- ""
txtMat[txtMat <0.05&txtMat >0.01] <- "*"
txtMat[txtMat <0.01&txtMat >0.005] <- "**"
#txtMat[txtMat <0.005&txtMat >0] <- "***"

txtMat1 <- signif( moduleTraitCor_DS,2)
#we only want to look at pearson correlations in certain range
txtMat1[txtMat1> -0.3&txtMat1<0.2] <- "" 

textMatrix1 = paste( txtMat1, '\n', '(',txtMat ,')', sep = '');
textMatrix1= matrix(textMatrix1,ncol=ncol( moduleTraitPvalue_DS),nrow=nrow(moduleTraitPvalue_DS))

pdf(paste('NetworkPlot_DS_consensus_112418.pdf'),width=16,height=30)
par( mar = c(8, 12, 3, 3) );
labeledHeatmap( Matrix = moduleTraitCor_DS,
                xLabels = colnames(factors1_DS),
                yLabels = rownames(moduleTraitPvalue_DS),
                ySymbols = rownames(moduleTraitPvalue_DS),
                colorLabels = FALSE,
                colors = blueWhiteRed(50),
                textMatrix = textMatrix1,
                setStdMargins = FALSE,
                cex.text = 1.5,
                zlim = c(-1, 1),
                cex.lab.x = 1.2,
                main = paste("Module-trait relationships")
);

#Plot eigengene heatmap
par(cex = 1.0)
plotEigengeneNetworks(MEs_DS, "Eigengene Network", marHeatmap = c(3,4,2,2), marDendro = c(0,4,1,2),cex.adjacency = 0.3,plotDendrograms = TRUE, xLabelsAngle = 90,heatmapColors=blueWhiteRed(100)[51:100])

#Plot boxplots, scatterplots
toplot=t(MEs_DS)
cols=substring(colnames(MEs_DS),3,20)
par(mfrow=c(4,4))
par(mar=c(5,6,4,2))
for (i in 1:nrow(toplot)) {
  boxplot(toplot[i,]~factor(as.vector(as.factor(targets.Ref.DS$disease)),c('CTL','DS')),col=cols[i],ylab="ME",main=rownames(toplot)[i],xlab=NULL,las=2)
  verboseScatterplot(x=as.numeric(targets.Ref.DS$Age),y=toplot[i,],xlab="Age",ylab="ME",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=cols[i],pch=19)
  boxplot(toplot[i,]~factor(targets.Ref.DS$Sex),col=cols[i],ylab="ME",main=rownames(toplot)[i],xlab=NULL,las=2)
}

dev.off()


moduleTraitPval=cbind(moduleTraitPvalue_DS, moduleTraitPvalue_AD1, moduleTraitPvalue_AD2)
moduleTraitCor=cbind(moduleTraitCor_DS, moduleTraitCor_AD1, moduleTraitCor_AD2)
write.csv(moduleTraitPval,'moduleTraitPvalue_DSAD_112418.csv')
write.csv(moduleTraitCor,'moduleTraitCor_DSAD_112418.csv')

#=====================================================================================
#
#  Part 6: GO analysis 
#
#=====================================================================================
# Makes bar plots of top enriched GO terms for each module

dir.create("./geneInfo")
dir.create("./geneInfo/background/")
dir.create("./geneInfo/input/")
dir.create("./geneInfo/output/")

geneInfo.cons$SystemCode =rep("En",length=nrow(geneInfo.cons))
background=geneInfo.cons[,"Ensembl.Gene.ID"]
background=as.data.frame(background)

## Output files for GO elite
background <- cbind(background,rep("En",length=length(background)))
colnames(background) <- c("Source Identifier","SystemCode")
write.table(background,"./geneInfo/background/denominator.txt",row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

uniquemodcolors=unique(moduleColor.cons)
uniquemodcolors=uniquemodcolors[uniquemodcolors!='grey']

# i = Number of modules
for(i in 1:length(uniquemodcolors)){
  thismod= uniquemodcolors[i]
  ind=which(colnames(geneInfo.cons)==paste("consensus.kME",thismod,sep=""))
  thisInfo=geneInfo.cons[geneInfo.cons$Initially.Assigned.Module.Color==thismod, c(1, 19, ind)] ##18=Ensembl.ID, 21="SystemCode",ind=kME value
  colnames(thisInfo) <- c("Source Identifier","SystemCode","kME")
  write.table(thisInfo,file=paste("./geneInfo/input/",thismod,"_Module.txt",sep=""),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
}

# Run GO elite as nohupped shell script:
codedir <- "/home/vivek/bin/GO-Elite_v.1.2.5-Py"
pathname <- "/home/erebboah/WGCNA/geneInfo"
nperm=10000
system(paste("nohup python ",codedir,"/GO_Elite.py --species Hs --mod Ensembl --permutations ",
             nperm,"  --method \"z-score\" --zscore 1.96 --pval 0.01 --num 5 --input ",pathname,
             "/input --denom ",pathname,"/background --output ",pathname,"/output &",sep=""))

# Plotting the GO Output
pathname <- "/home/erebboah/WGCNA/geneInfo/output/GO-Elite_results/CompleteResults"

uniquemodcolors=uniquemodcolors[-c(2, 14)] # For some reason sometimes modules are not run correctly, therefore won't be able to be plotted so they are excluded

pdf("GOElite_plot_Modules_112718.pdf",height=8,width=12)
for(i in 1:length(uniquemodcolors)){
  thismod= uniquemodcolors[i]
  tmp=read.csv(file=paste(pathname,"/ORA_pruned/",thismod,"_Module-GO_z-score_elite.txt",sep=""),sep="\t")
  tmp=subset(tmp,Ontology.Type!='cellular_component')
  tmp=tmp[,c(2,9)] ## Select GO-terms and Z-score
  tmp=tmp[order(tmp$Z.Score,decreasing=T),] #
  if (nrow(tmp)<10){
    tmp1=tmp ## Take top 10 Z-score
    tmp1 = tmp1[order(tmp1$Z.Score),] ##Re-arrange by increasing Z-score
    par(mar=c(5,40,5,2))
    barplot(tmp1$Z.Score,horiz=T,col="blue",names.arg= tmp1$Ontology.Name,cex.names=1.2,las=1,main=paste("Gene Ontology Plot of",thismod,"Module"),xlab="Z-Score")
    abline(v=2,col="red")
  } else {
    tmp1=tmp[c(1:10),] ## Take top 10 Z-score
    tmp1 = tmp1[order(tmp1$Z.Score),] ##Re-arrange by increasing Z-score
    par(mar=c(5,40,5,2))
    barplot(tmp1$Z.Score,horiz=T,col="blue",names.arg= tmp1$Ontology.Name,cex.names=1.2,las=1,main=paste("Gene Ontology Plot of",thismod,"Module"),xlab="Z-Score")
    abline(v=2,col="red")
    }
  
  cat('Done ...',thismod,'\n')
}

dev.off()

#=====================================================================================
#
#  Part 7: Cell type enrichment 
#
#=====================================================================================
geneInfo=read.csv('geneInfo.cons.DSAD_112418.csv')

datKME <- geneInfo[,c("Ensembl.Gene.ID","Initially.Assigned.Module.Color")] ## Get a list of genes to test for enrichment, e.g. genes with modules defined
testbackground <- as.character(geneInfo$Ensembl.Gene.ID) # background list
datKME=subset(datKME,Initially.Assigned.Module.Color!="grey")
namestestlist <- names(table(datKME[,2])) ## module
multiTest <- vector(mode = "list", length = length(namestestlist))
names(multiTest) <- namestestlist


for (i in 1:length(multiTest))
{
  multiTest[[i]] <- datKME[datKME[,2]==namestestlist[i],1]
}

datCells <- read.csv("/home/vivek//bin/ZhangEtAlCellTypeList_humanENSG.csv") ## From Zhang et al., 2014 - ## CSV file with reference genes in 1st column, annotated list with Category in 2nd

## Set up reference lists
namesreflist <- names(table(datCells[,2])) ## category or module color
multiRef <- vector(mode = "list", length = length(namesreflist))
names(multiRef) <- namesreflist
for (i in 1:length(multiRef))
{
  multiRef[[i]] <- datCells[datCells[,2]==namesreflist[i],3]
}

refbackground<- testbackground

source('/home/vivek//AD/Zhang/ORA.R')

ORA.OR = matrix(NA,nrow=length(multiTest),ncol=length(multiRef));
colnames(ORA.OR) = names(multiRef);
rownames(ORA.OR) = names(multiTest);
ORA.P = matrix(NA,nrow=length(multiTest),ncol=length(multiRef));
colnames(ORA.P) = names(multiRef);
rownames(ORA.P) = names(multiTest);

for (i in 1:length(multiRef)) {
  for (j in 1:length(multiTest)) {
    result = ORA(multiTest[[j]],multiRef[[i]],testbackground,refbackground);
    ORA.OR[j,i] = result[1];
    ORA.P[j,i] = result[2];
  }
}

ORA.OR<-apply(ORA.OR,2,as.numeric)
dim(ORA.OR)<-dim(ORA.P)

##FDR correct

FDRmat.Array <- matrix(p.adjust( ORA.P,method="fdr"),nrow=nrow( ORA.P),ncol=ncol( ORA.P))
rownames(  FDRmat.Array)=rownames(ORA.P)
colnames(  FDRmat.Array)=colnames(ORA.P)

ORA.P=matrix(as.numeric(ORA.P),nrow=nrow( ORA.P),ncol=ncol( ORA.P))
ORA.OR=matrix(as.numeric(ORA.OR),nrow=nrow( ORA.OR),ncol=ncol( ORA.OR))
rownames(ORA.P) <- rownames(ORA.OR) <- rownames(  FDRmat.Array)
colnames(ORA.P) <- colnames(ORA.OR) <- colnames(  FDRmat.Array)

dispMat <- ORA.OR ## You can change this to be just log2(Bmat) if you want the color to reflect the odds ratios
#Use the text function with the FDR filter in labeledHeatmap to add asterisks
txtMat <-  ORA.OR
txtMat[FDRmat.Array >0.05] <- ""
txtMat[FDRmat.Array <0.05&FDRmat.Array >0.01] <- "*"
txtMat[FDRmat.Array <0.01&FDRmat.Array >0.005] <- "**"
txtMat[FDRmat.Array <0.005] <- "***"

txtMat1 <- signif( ORA.OR,2)
txtMat1[txtMat1<2] <- ""

textMatrix1 = paste( txtMat1, '\n', txtMat , sep = '');
textMatrix1= matrix(textMatrix1,ncol=ncol( ORA.P),nrow=nrow( ORA.P))

# Make heatmap of modules and cell types: neurons, microglia, myelinating oligodendrocytes, astrocytes, and endothelial cells
pdf("CellTypeEnrich_WGCNAMods_110618.pdf", width=6,height=10)
labeledHeatmap(Matrix=dispMat,
               yLabels=rownames(dispMat),
               yColorLabels=TRUE,
               xLabels= colnames(dispMat),
               colors=blueWhiteRed(40),
               textMatrix = textMatrix1,
               cex.lab.x=1.0,
               zlim=c(-0.1,3),
               main="Cell-type enrichment Heatmap")
dev.off()


#=====================================================================================
#
#  Part 8: TOM network plot
#
#=====================================================================================
# Will make node and edge plots with hub genes in the center surrounded by all other genes in each module

load("geneInfo.cons.112418.rda")
load("consensus_112418.rda")

#Get the top connected genes in the module
uniquemodcolors = unique(moduleColor.cons);
uniquemodcolors <- uniquemodcolors[!uniquemodcolors %in% "grey"]
TOM.matrix = as.matrix(consTomDS);

pdf("ModuleNetworks.pdf",height=9,width=10);
for (mod in uniquemodcolors)  {
  numgenesingraph = 100;
  numconnections2keep = 1500;
  cat('module:',mod,'\n');
  geneInfo.cons=geneInfo.cons[geneInfo.cons$GeneSymbol!="NA",]
  colind = which(colnames(geneInfo.cons)==paste("consensus.kME",mod, sep=""));
  rowind = which(geneInfo.cons[,3]==mod);
  cat(' ',length(rowind),'probes in module\n');
  submatrix = geneInfo.cons[rowind,];
  orderind = order(submatrix[,colind],decreasing=TRUE);
  if (length(rowind) < numgenesingraph) {
    numgenesingraph = length(rowind);
    numconnections2keep = numgenesingraph * (numgenesingraph - 1);
  }
  cat('Making network graphs, using top',numgenesingraph,'probes and',numconnections2keep,'connections of TOM\n');
  submatrix = submatrix[orderind[1:numgenesingraph],];
  #Identify the columns in the TOM that correspond to these hub probes
  matchind = match(submatrix$Ensembl.Gene.ID,colnames(multiExpr));
  reducedTOM = TOM.matrix[matchind,matchind];
  
  orderind = order(reducedTOM,decreasing=TRUE);
  connections2keep = orderind[1:numconnections2keep];
  reducedTOM = matrix(0,nrow(reducedTOM),ncol(reducedTOM));
  reducedTOM[connections2keep] = 1;
  
  g0 <- graph.adjacency(as.matrix(reducedTOM[1:10,1:10]),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutMata <- layout.circle(g0)
  
  g0 <- graph.adjacency(as.matrix(reducedTOM[11:50,11:50]),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutMatb <- layout.circle(g0)
  
  g0 <- graph.adjacency(as.matrix(reducedTOM[51:ncol(reducedTOM),51:ncol(reducedTOM)]),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutMatc <- layout.circle(g0)
  g1 <- graph.adjacency(as.matrix(reducedTOM),mode="undirected",weighted=TRUE,diag=FALSE)
  layoutMat <- rbind(layoutMata*0.25,layoutMatb*0.8, layoutMatc)

  plot(g1,edge.color="grey",vertex.color=mod,vertex.label=as.character(submatrix$GeneSymbol),vertex.label.cex=0.7,vertex.label.dist=0.45,vertex.label.degree=-pi/4,vertex.label.color="black",layout= layoutMat,vertex.size=submatrix[,colind]^2*8,main=paste(mod,"module"))
  
}
dev.off();
