rm(list = ls())
setwd("J:/V/R")
# CODE FROM HERE 
METHOD-1 up to 790 lines: #https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html 
install.packages("BiocManager")
BiocManager::install("WGCNA")
BiocManager::install("GO.db")
install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival", "BiocManager"))
BiocManager::install(c("GO.db", "preprocessCore", "impute"));
library(GO.db)
library(WGCNA)
library(DESeq2)
library(pheatmap)
library(tidyr)
library(ggplot2)
library(RCy3)
library(GO.db)
#BiocManager::install("org.Mm.eg.db")
library("org.Mm.eg.db")
1. Data input and cleaning
#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in the female liver data set
femData = read.csv("KIRC1.csv", header = T, row.names = 1);
# Take a quick look at what is in the data set:
dim(femData);
names(femData);
datExpr0=log2(femData+1)
datExpr0=t(datExpr0)
class(datExpr0)
datExpr0 = as.data.frame(datExpr0)
#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================
#datExpr0 = as.data.frame(t(femData[, -c(1:8)]));
#names(datExpr0) = femData$substanceBXH;
#rownames(datExpr0) = names(femData)[-c(1:8)];
#=====================================================================================
#
#  Code chunk 3 :Check Data for excessive missingness
#You should get TRUE for this dataset given the parameters above. TRUE means that all clusters have passed the cut. This means that when you limited your clusters to the most common ones above that you didn't leave any in that had too many zeros. It is still possible that you were too choosy though. If you got FALSE for your data then you have to follow some other steps that I don't go over here.
#=====================================================================================
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
#=====================================================================================
#
#  Code chunk 5
# Eucledean diatance: Graph determine how many clusters are within the dendrogram and cut the dendrogram at a certain tree height to separate the data into different groups.
#=====================================================================================
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================
# Plot a line to show the cut
abline(h = 100, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 100, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#The variable datExpr now contains the expression data ready for network analysis.
#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================
# Must be numeric data
traitData = read.csv("cesc miR.csv", row.names = 1);
dim(traitData)
names(traitData)
# remove columns that hold information we do not need.
#allTraits = traitData[, -c(31, 16)];
#allTraits = allTraits[, c(2, 11:36) ];
#dim(allTraits)
#names(allTraits)
# Form a data frame analogous to expression data that will hold the clinical traits.
femaleSamples = rownames(datExpr);
traitRows = match(femaleSamples, row.names(traitData));
# traitRows extraction from all traits data of matched data from expression profile 
datTraits = traitData[traitRows,];
#datTraits[,1]
#rownames(datTraits)<-datTraits[,1] # we have to replace row names 
#datTraits[,1]
#View(datTraits)
#datTraits<-datTraits[,-1]
View(datTraits)

collectGarbage();
#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
# white means a low value, red means a high value, and gray means missing entry.
#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================
save(datExpr, datTraits, file = "FemaleLiver-01-dataInput_kirc100genes.RData")

2. Network construction and module detection
this is common for a, b, c 
0 Preliminaries: setting up the R session
2.a Automatic network construction and module detection
2.a.1 Choosing the soft-thresholding power: analysis of network topology
2.b Step-by-step network construction and module detection
2.b.1 Choosing the soft-thresholding power: analysis of network topology
2.c Automatic block-wise network construction and module detection
2.c.1 Choosing the soft-thresholding power: analysis of network topology
#=====================================================================================
#
#  Code chunk 9 Preliminaries: setting up the R session
#
#=====================================================================================
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "FemaleLiver-01-dataInput_kirc100genes.RData");
#The variable lnames contains the names of loaded variables.
lnames
#=====================================================================================
#
#  Code chunk 2 Automatic construction of the gene network and identification of modules
#
#=====================================================================================
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
#We choose the power 6, which is the lowest power for which the scale-free topology fit index curve flattens out upon reaching a high value (in this case, roughly 0.90).

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

2.a.2 One-step network construction and module detection
#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================
power=sft$powerEstimate 
cor <- WGCNA::cor
net = blockwiseModules(datExpr, power = power,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM", 
                       verbose = 3)

cor<- stats::cor
#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "FemaleLiver-02-networkConstruction-auto_kirc_1000genes.RData")

2.b.2 Co-expression similarity and adjacency
#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================
softPower = 5;
adjacency = adjacency(datExpr, power = softPower);

2.b.3 Topological Overlap Matrix (TOM)
#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================
2.b.4 Clustering using TOM
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================
2.b.5 Merging of modules whose expression profiles are very similar
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================
sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()
#=====================================================================================
#
#  Code chunk 11
#
#=====================================================================================
# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "FemaleLiver-02-networkConstruction-stepByStep.RData")


2.c.2 Block-wise network construction and module detection
#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================
bwnet = blockwiseModules(datExpr, maxBlockSize = 2000,
                         power = 6, TOMType = "unsigned", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "femaleMouseTOM-blockwise_kirc",
                         verbose = 3)
#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================
# Load the results of single-block analysis
load(file = "FemaleLiver-02-networkConstruction-auto_kirc.RData");
# Relabel blockwise modules
bwLabels = matchLabels(bwnet$colors, moduleLabels);
# Convert labels to colors for plotting
bwModuleColors = labels2colors(bwLabels)
#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================
# open a graphics window
sizeGrWindow(6,6)
# Plot the dendrogram and the module colors underneath for block 1
plotDendroAndColors(bwnet$dendrograms[[1]], bwModuleColors[bwnet$blockGenes[[1]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 1", 
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for block 2
plotDendroAndColors(bwnet$dendrograms[[2]], bwModuleColors[bwnet$blockGenes[[2]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 2", 
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================
2.c.3 Comparing the single block and block-wise network analysis
sizeGrWindow(12,9)
plotDendroAndColors(geneTree,
                    cbind(moduleColors, bwModuleColors),
                    c("Single block", "2 blocks"),
                    main = "Single block gene dendrogram and module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================
singleBlockMEs = moduleEigengenes(datExpr, moduleColors)$eigengenes;
blockwiseMEs = moduleEigengenes(datExpr, bwModuleColors)$eigengenes;
#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================
#Next we match the single-block and block-wise eigengenes by name and calculate the correlations of the corresponding eigengenes:
single2blockwise = match(names(singleBlockMEs), names(blockwiseMEs))
signif(diag(cor(blockwiseMEs[, single2blockwise], singleBlockMEs)), 3)

3.Relating modules to external clinical traits and identifying important genes
#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "FemaleLiver-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "FemaleLiver-02-networkConstruction-auto.RData");
lnames
#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================
3 Relating modules to external clinical traits
3.a Quantifying module-trait associations
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
#Each row corresponds to a module eigengene, column to a trait. 
#Each cell contains the corresponding correlation and p-value. The table is color-coded by correlation according to the color legend.
#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================
3.b Gene relationship to trait and important modules: Gene Significance and Module
Membership
# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$weight_g);
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");
#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================
3.c Intramodular analysis: identifying genes with high GS and MM
module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================
3.d Summary output of network analysis results
names(datExpr)
#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================
names(datExpr)[moduleColors=="brown"]
#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================
annot = read.csv(file = "GeneAnnotation.csv");
dim(annot)
names(annot)
probes = names(datExpr)
probes2annot = match(probes, annot$substanceBXH)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.

#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================
# Create the starting data frame
geneInfo0 = data.frame(substanceBXH = probes,
                       geneSymbol = annot$gene_symbol[probes2annot],
                       LocusLinkID = annot$LocusLinkID[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight));
geneInfo = geneInfo0[geneOrder, ]
#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================
write.csv(geneInfo, file = "geneInfo.csv")

4. Interfacing network analysis with other data such as functional annotation and gene ontology 
4.a Output gene lists for use with online software and services#
# Read in the probe annotation
annot = read.csv(file = "GeneAnnotation.csv");
# Match probes in the data set to the probe IDs in the annotation file 
probes = names(datExpr)
probes2annot = match(probes, annot$substanceBXH)
# Get the corresponding Locuis Link IDs
allLLIDs = annot$LocusLinkID[probes2annot];
# $ Choose interesting modules
intModules = c("brown", "red", "salmon")
for (module in intModules)
{
  # Select module probes
  modGenes = (moduleColors==module)
  # Get their entrez ID codes
  modLLIDs = allLLIDs[modGenes];
  # Write them into a file
  fileName = paste("LocusLinkIDs-", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE)
}
# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("LocusLinkIDs-all.txt", sep="");
write.table(as.data.frame(allLLIDs), file = fileName,
            row.names = FALSE, col.names = FALSE)


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

4.b Enrichment analysis directly within R
GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, organism = "mouse", nBestP = 10);
#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================
tab = GOenr$bestPTerms[[4]]$enrichment
#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================
names(tab)
#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================
write.table(tab, file = "GOEnrichmentTable.csv", sep = ",", quote = TRUE, row.names = FALSE)
#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================
keepCols = c(1, 2, 5, 6, 7, 12, 13);
screenTab = tab[, keepCols];
# Round the numeric columns to 2 decimal places:
numCols = c(3, 4);
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)
# Truncate the the term name to at most 40 characters
screenTab[, 7] = substring(screenTab[, 7], 1, 40)
# Shorten the column names:
colnames(screenTab) = c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name");
rownames(screenTab) = NULL;
# Set the width of R's output. The reader should play with this number to obtain satisfactory output.
options(width=95)
# Finally, display the enrichment table:
screenTab

5. Network visualization using WGCNA functions
5.a Visualizing the gene network
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "FemaleLiver-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "FemaleLiver-02-networkConstruction-auto.RData");
lnames
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================
5.b Visualizing the network of eigengenes
nSelect = 400
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================
# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate weight from the clinical traits
weight = as.data.frame(datTraits$weight_g);
names(weight) = "weight"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, weight))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================
# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)

6. Exporting a gene network to external visualization software
#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================
6.a Exporting to VisANT
# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr, power = 6);
# Read in the annotation file
annot = read.csv(file = "GeneAnnotation.csv");
# Select module
module = "brown";
# Select module probes
probes = names(datExpr)
inModule = (moduleColors==module);
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol) )
#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================
nTop = 30;
IMConn = softConnectivity(datExpr[, modProbes]);
top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(modTOM[top, top],
                            file = paste("VisANTInput-", module, "-top30.txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol) )
#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================
6.b Exporting to Cytoscape
# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 6);
# Read in the annotation file
annot = read.csv(file = "GeneAnnotation.csv");
# Select modules
modules = c("brown", "red");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])


Method-2:
library(WGCNA)
library(DESeq2)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#===============================================================================
#
#  Download the files
#
#==============================================================================
# Step 1: Download files from https://github.com/PengYuMaize/Yu2021NaturePlants
# Rename the sample ID in the "gene_counts_table" and "sample_info" files, so that they look alike
# Rename for the sample ID start with "787_" and "78371A_"
# or you may just download my modified "gene_counts_table_WGCNA_LC.txt" for smooth tutorial experience
#===============================================================================
#
#  Read the gene counts table and plot the sample tree
#
#===============================================================================
# Read the gene counts table 
data0=read.csv("KIRC1.csv",header=T,row.names=1)
datExpr=t(log2(data0+1))
datExpr=as.data.frame(datExpr)
# Normalization with log2(FPKM+1)
sample_metadata = read.csv(file = "clinical m0 and m1.csv")
b<-toupper(sample_metadata$s) 
#dataExpr_deseq <- DESeqDataSetFromMatrix(countData = data0[,-181],colData = sample_metadata,design = ~ Zone)
#mcols(dataExpr_deseq)$basepairs = data0$geneLengt1
#fpkm_matrix = fpkm(dataExpr_deseq)
#datExpr = t(log2(fpkm_matrix+1))

head(datExpr[1:5,1:5]) # samples in row, genes in column

match(sample_metadata$, colnames(data0))
datExpr <- datExpr[,1:5000]

gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK
#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr[gsg$goodSamples, gsg$goodGenes]
}
# Calculate sample distance and cluster the samples
sampleTree = hclust(dist(datExpr), method = "average");

Plot-1 sample tree
# plot sample tree
#pdf(file = "1-n-sampleClustering.pdf", width = 40, height = 20);
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
dev.off()

# Plot a line to show the cut
abline(h = 110, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 110, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr1 = datExpr[keepSamples, ]
nGenes = ncol(datExpr1)
nSamples = nrow(datExpr1)
#===============================================================================
#
#  Choose soft threshold parameter
#
#===============================================================================

# Choose a set of soft threshold parameters
powers = c(c(1:10), seq(from = 10, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5) 

Plot-2 Scale and mean connectivity
# Scale-free topology fit index as a function of the soft-thresholding power
#pdf(file = "2-n-sft.pdf", width = 9, height = 5);
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red") 
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#===============================================================================
#
#  Turn data expression into topological overlap matrix
#
#===============================================================================

# Turn data expression into topological overlap matrix
power=sft$powerEstimate #4

# Option 1: automatic
cor <- WGCNA::cor
net = blockwiseModules(datExpr, power = power,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = FALSE,
                       verbose = 3)
cor<- stats::cor
# unsigned -> nodes with positive & negative correlation are treated equally 
# signed -> nodes with negative correlation are considered *unconnected*, treated as zero

sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
Plot-3 Dendogram in pdf
pdf(file = "4-module_tree_blockwise.pdf", width = 8, height = 6);
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

################################################################################
################################################################################
# Option 2a: step-by-step
power = power
adjacency = adjacency(datExpr, power = power)
TOM = TOMsimilarity(adjacency); # Turn adjacency into topological overlap
dissTOM = 1-TOM


# Option 2b: 
#TOM = TOMsimilarityFromExpr(datExpr, power = power)
#dissTOM = 1-TOM 
#dim(dissTOM)


#===============================================================================
#
#  Construct modules (proceed with the genetree from option 2b)
#
#===============================================================================
Plot-4: dendogram  in 
# Plot gene tree
geneTree = hclust(as.dist(dissTOM), method = "average");
#pdf(file = "3-gene_cluster.pdf", width = 12, height = 9);
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()

# Module identification using dynamic tree cut
# We like large modules, so we set the minimum module size relatively high:
# minModuleSize = 30;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, 
                            pamRespectsDendro = FALSE,minClusterSize = 30);
table(dynamicMods)
length(table(dynamicMods)) 
# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
pdf(file = "4-module_tree.pdf", width = 8, height = 6);
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE,
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
dev.off()

#===============================================================================
#
#  Merge modules
#
#===============================================================================
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
Plot-5: in clustering 
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# Merge close modules
MEDissThres=0.40
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3) 
mergedColors = merge$colors  
mergedMEs = merge$newMEs  
# Plot merged module tree
#pdf(file = "5-merged_Module_Tree.pdf", width = 12, height = 9)  
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)  
dev.off()
#write.table(merge$oldMEs,file="oldMEs.txt");
#write.table(merge$newMEs,file="newMEs.txt");

#===============================================================================
#
#  Export of networks to external software
#
#===============================================================================

# Export the gene list of old modules 
for (i in 1:length(merge$oldMEs)){
  modules = c(substring(names(merge$oldMEs)[i], 3));
  genes = colnames(datExpr)
  inModule = is.finite(match(dynamicColors,modules))
  modGenes = genes[inModule]
  modTOM=TOM[inModule,inModule]
  dimnames(modTOM)=list(modGenes,modGenes)
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("output_for_cytoscape/orign_CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("output_for_cytoscape/orign_CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE, threshold = -1, nodeNames = modGenes, nodeAttr = dynamicColors[inModule]);
}
# Export the gene list of new modules 
for (i in 1:length(merge$newMEs)){
  modules = c(substring(names(merge$newMEs)[i], 3));
  genes = colnames(datExpr)
  inModule = is.finite(match(dynamicColors,modules))
  modGenes = genes[inModule]
  modTOM=TOM[inModule,inModule]
  dimnames(modTOM)=list(modGenes,modGenes)
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("output_for_cytoscape/merge_CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("output_for_cytoscape/merge_CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE, threshold = -1, nodeNames = modGenes, nodeAttr = dynamicColors[inModule]);
}

#===============================================================================
#
#  PART 1: Correlate module eigen-genes and samples (or other discrete data)
#
#===============================================================================
# Heatmap of old module eigen-genes and samples
#pdf(file="oldMEs.pdf",heigh=80,width=20)
library("pheatmap")
rownames(merge$oldMEs)=names(data0[,-181])
pheatmap(merge$oldMEs,cluster_col=T,cluster_row=T,show_rownames=F,show_colnames=T,fontsize=6)
dev.off()


# Heatmap of new module eigen-genes and sample trait (e.g. Zone)
col_ann <- sample_metadata[,c(1,3)]
rownames(col_ann) <- col_ann[,1]
col_ann <- data.frame(col_ann)
col_ann$Zone <- as.factor(col_ann$Zone)
col_ann <- col_ann[order(col_ann$Zone),]
col_ann$sample_ID <- NULL
head(col_ann)
ann_color <- list("col_ann" = c("Z1" = "yellow",
                                "Z2" = "red",
                                "Z3" = "green"))

data <- data.frame(merge$newMEs)
data <- data[order(match(rownames(data), rownames(col_ann))),]
dim(merge$newMEs)

#pdf(file="newMEs.pdf",heigh=60,width=20)
rownames(merge$newMEs)=names(data0[,-181])
pheatmap(data,cluster_col=T,cluster_row=F,show_rownames=F,
         show_colnames=T,fontsize=6,
         annotation_row = col_ann, annotation_colors = ann_color)
dev.off()

#=====================================================================================
#
#  PART 2: Correlation between gene modules and microbial traits (continuous data)
#
#=====================================================================================

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, mergedColors)$eigengenes
MEs = orderMEs(MEs0)

# Read microbial data as traits
bac_traits = read.table("b_order_234.txt", header = T, sep = "\t")
rownames(bac_traits) = bac_traits[, 1]
bac_traits = bac_traits[, -1]
# sample names should be consistent in eigen genes and traits !!!!
bac_traits = bac_traits[match(rownames(MEs), rownames(bac_traits)), ]
table(rownames(MEs) == rownames(bac_traits))

# Calculate pearson correlation coefficients between module eigen-genes and traits
moduleTraitCor = cor(MEs, bac_traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
#write.table(moduleTraitCor,file="moduleTrait_correlation.txt");
#write.table(moduleTraitPvalue,file="moduleTrait_pValue.txt");


#  Plot heatmap of module-traits relationship
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf("5-module-traits-bacteria-order.pdf", width = 80, height = 15)
par(mar = c(15, 12, 5, 5));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(bac_traits),
               yLabels = colnames(MEs),
               ySymbols = colnames(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#  Plot heatmap of module-traits relationship
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor[1:25,1:25], 2), "\n(",
                    signif(moduleTraitPvalue[1:25,1:25], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[1:25,1:25])
pdf("5-module-traits-bacteria-order1.pdf", width = 20, height = 10)
par(mar = c(15, 12, 5, 5));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor[1:25,1:25],
               xLabels = colnames(bac_traits[1:25,1:25]),
               yLabels = colnames(MEs[1:25,1:25]),
               ySymbols = colnames(MEs[1:25,1:25]),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#=====================================================================================
#
#   Intramodular analysis: identifying genes with high geneModuleMembership & geneTraitSignficance
#
#=====================================================================================


# Define variable Verru containing the Verrucomicrobiales column of bac_traits
Verru = as.data.frame(bac_traits$Verrucomicrobiales);
names(Verru) = "Verrucomicrobiales"

# names (colors) of the modules
modNames = substring(names(MEs), 3)

MET = orderMEs(cbind(MEs, Verru))

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, Verru, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Verru), sep="");
names(GSPvalue) = paste("p.GS.", names(Verru), sep="");

# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(10,10,1,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)


module = "lightgreen"
# Rename to moduleColors
moduleColors = mergedColors
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Verrucomicrobiales",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "blue")


## Draw bubble plot for particular module
colsum_bac_traits <- colSums(bac_traits)
colsum_bac_traits <- data.frame(colsum_bac_traits)
colsum_bac_traits$b_order <- rownames(colsum_bac_traits)
library(tidyr)
moduleTraitCor_long <- data.frame(moduleTraitCor)
moduleTraitCor_long$module <- rownames(moduleTraitCor)
moduleTraitCor_long <- moduleTraitCor_long[,c(235,1:234)]
moduleTraitCor_long <- gather(moduleTraitCor_long, b_order, PCC, Pseudomonadales:Others, factor_key = TRUE)

moduleTraitPvalue_long <- data.frame(moduleTraitPvalue)
moduleTraitPvalue_long$module <- rownames(moduleTraitPvalue)
moduleTraitPvalue_long <- moduleTraitPvalue_long[,c(235,1:234)]
moduleTraitPvalue_long <- gather(moduleTraitPvalue_long, b_order, pval, Pseudomonadales:Others, factor_key = TRUE)

moduleTrait_long <- merge(moduleTraitCor_long, moduleTraitPvalue_long, by = c("module","b_order"))

bubble_Data <- merge(moduleTrait_long, colsum_bac_traits, by = "b_order")
#just want module = "lightgreen"
bubble_Data_lightgreen <- bubble_Data[which(bubble_Data$module == "MElightgreen"),]

library(ggplot2)
ggplot(bubble_Data_lightgreen, aes(x= colsum_bac_traits, y= PCC, size = colsum_bac_traits,
                                   color = PCC, label = b_order)) +
  geom_text(hjust = 1, size=3) +
  geom_point(alpha=1) + ylab("Module-taxon correlation") + xlab("Relative abundance (sum)") +
  theme_bw()


############# Summary ###################################

head(datExpr)[1:5,1:5] # transcriptome data

head(sample_metadata)[1:5,] # metadata (sample info)
head(bac_traits)[1:5,1:5] # external trait






#=====================================================================================
#
#   Cytoscape
#
#=====================================================================================


#if(!"RCy3" %in% installed.packages()){
 install.packages("BiocManager")
 BiocManager::install("RCy3")
#}

# https://cytoscape.org/cytoscape-automation/for-scripters/R/notebooks/
library(RCy3)

cytoscapePing () # make sure cytoscape is open
cytoscapeVersionInfo ()

###### for yellow module of the merged data (newMEs) #################################
edge <- read.delim("output_for_cytoscape/merge_CytoscapeInput-edges-lightgreen.txt")
colnames(edge)
colnames(edge) <- c("source", "target","weight","direction","fromAltName","toAltName")

node <- read.delim("output_for_cytoscape/merge_CytoscapeInput-nodes-lightgreen.txt")
colnames(node)  
colnames(node) <- c("id","altName","node_attributes") 

createNetworkFromDataFrames(node,edge[1:50,], title="my first network", collection="DataFrame Example")

################ customise the network visualization ##################################
# use other pre-set visual style
setVisualStyle('Marquee')

# set up my own style
style.name = "myStyle"
defaults <- list(NODE_SHAPE="diamond",
                 NODE_SIZE=30,
                 EDGE_TRANSPARENCY=120,
                 NODE_LABEL_POSITION="W,E,c,0.00,0.00")
nodeLabels <- mapVisualProperty('node label','id','p')
nodeFills <- mapVisualProperty('node fill color','node_attributes','d',c("A","B"), c("#FF9900","#66AAAA"))
arrowShapes <- mapVisualProperty('Edge Target Arrow Shape','interaction','d',c("activates","inhibits","interacts"),c("Arrow","T","None"))
edgeWidth <- mapVisualProperty('edge width','weight','p')

createVisualStyle(style.name, defaults, list(nodeLabels,nodeFills,arrowShapes,edgeWidth))
setVisualStyle(style.name)