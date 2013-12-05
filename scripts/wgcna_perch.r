# Options
rm(list=ls())
library(WGCNA)
options(stringsAsFactors = FALSE)



# Global variables
expression.file = "OUTPUT_lowess_data_for_wgcna.csv"
design.file = "files_to_fill/design.txt"
prepared.data = "OUTPUT_wgcna_prepared_data.RData"
wgcna.modules = "OUTPUT_wgcna_modules.RData"
tom.basename = "OUTPUT_wgcna_tom"

# Reading data
# tutorial file 1
data = read.csv(expression.file, sep="\t")
dim(data)
names(data)

# Formatting data
data.expression.temp = as.data.frame(t(data[, -1]))
names(data.expression.temp) = data$Gene
rownames(data.expression.temp) = names(data)[-1]

## Check for excessive missing values
gsg = goodSamplesGenes(data.expression.temp, verbose = 3)
gsg$allOK

# Remove gene and samples with too many missing values
#if (!gsg$allOK)
#{
#    # Optionally, print the gene and sample names that were removed:
#    if (sum(!gsg$goodGenes)>0)
#    printFlush(paste("Removing genes:", paste(names(data.expression.temp)[!gsg$goodGenes], collapse = ", ")))
#    if (sum(!gsg$goodSamples)>0)
#    printFlush(paste("Removing samples:", paste(rownames(data.expression.temp)[!gsg$goodSamples], collapse = ", ")))
#    # Remove the offending genes and samples from the data:
#    data.expression.temp = data.expression.temp[gsg$goodSamples, gsg$goodGenes]
#}

# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
jpeg("sample_clustering_and_outliers.jpg", width=1200, height=1000, quality=100)
    ## Cluster samples to detect outliers
    sampleTree = flashClust(dist(data.expression.temp), method = "average")
    par(cex = 0.6)
    par(mar = c(0,4,2,0))
    plot(sampleTree, main="Sample clustering", sub="", xlab="", cex.lab=1.5, cex.axis=1.5, cex.main=2)

    # Choose height threshold to remove potential outliers
    # Plot a line to show the cut
    cut.threshold = 25
    abline(h = cut.threshold, col = "red")
dev.off()

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = cut.threshold, minSize = 10)
table(clust)

# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = data.expression.temp[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
nGenes
nSamples

## Loading trait data
traitData = read.csv(design.file, sep="\t")

# WARNING: Change this for your project
traitData = traitData[,c("Fish_ID", "Temperature", "Ni", "Cd")]

# Eliminate data from removed samples
allTraits = traitData[keepSamples, ]
dim(allTraits)
names(allTraits)

# Form a data.frame analogous to expression data that will hold the traits.
perchSamples = rownames(datExpr)
traitRows = match(perchSamples, allTraits$Ind_name)
datTraits = allTraits
rownames(datTraits) = allTraits$Ind_name
datTraits = datTraits[, c("Temperature", "Ni", "Cd")]
collectGarbage()

## Visualize how traits to the sample dendogram
# Re-cluster samples
sampleTree2 = flashClust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE)
#traitColors = numbers2colors(traitRows, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
sizeGrWindow(10.25, 5.75)
plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap")

# Saving the data
save(datExpr, datTraits, file = prepared.data)



### START OF AUTOMATIC NETWORK CONSTRUCTION ###
# tutorial file 2

## Preliminaries
library(WGCNA)
options(stringsAsFactors = FALSE)

# Optional multi-threading within WGCNA. This helps speed up # certain calculations
# Any error here may be ignored but you may want to update WGCNA if you see one.
allowWGCNAThreads()

# Load the data saved in the first part
lnames = load(file = prepared.data)

#The variable lnames contains the names of loaded variables.
lnames

## Construction of the gene network and identification of modules
# Choose a set of soft-thresholding powers
powers = c(c(1:20))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(10.25, 5.75)
par(mfrow = c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red")

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")





# CAUTION Be sure to explore the parameters of blockwiseModules! There are more!
power = 7
minSize = 15
reasignThreshold = 0.2
deepSplit = 2
mergeHeight= 0.10
net = blockwiseModules(datExpr, power = power, minModuleSize = minSize, reassignThreshold = reasignThreshold, deepSplit=deepSplit, mergeCutHeight = mergeHeight, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = tom.basename, verbose = 2)

## Plot the result
# open a graphics window
#sizeGrWindow(10.25, 5.75)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
jpeg("module_dentogram.jpg", width=800, height=400, quality=100)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
table(net$colors)

geneTree = net$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree, file = wgcna.modules)


### START OF MODULE TO TRAIT RELATIONSHIP ###
# tutorial file 3

## Preliminaries
library(WGCNA)
options(stringsAsFactors = FALSE)
lnames = load(file = prepared.data)
lnames

## Relating modules to external traits
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

## Module-trait association
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Plot modules
jpeg("module-trait_relationships.jpg", height=800, width=400, quality=100)
    #sizeGrWindow(5.178, 11.79)
    # Will display correlations and their p-values
    textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
    dim(textMatrix) = dim(moduleTraitCor)
    par(mar = c(6, 8.5, 3, 3))
    # Display the correlation values within a heatmap plot
    labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits), yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE, colors = greenWhiteRed(50), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 1, zlim = c(-1,1), main = paste("Module-trait relationships"))
dev.off()



pvalues_all_traits = apply(moduleTraitPvalue, 1, min)
sign_modules = names(pvalues_all_traits)[pvalues_all_traits <= 0.05]
num_sign_modules = length(sign_modules)

sign_modules
pvalues_all_traits
cat("\n  Number of modules", length(MEs), "\n")
cat("  Number of significant modules", num_sign_modules, "\n")





#sign_modules=c("green")

pvalues = data.frame(moduleTraitPvalue)
genes = data.frame(1:1008)

for (trait in c("Temperature", "Ni", "Cd")){
    sign_modules = row.names(pvalues)[pvalues[,trait] <= 0.05]
    for (m in sign_modules){
        module = gsub("ME", "", m)
        cat(module, " module for trait ", trait, "\n", sep="")
        ## Gene-trait association
        # Define variable weight containing the weight column of datTrait
        weight = as.data.frame(datTraits[,trait]) #||| Trait |||#
        names(weight) = trait
        # names (colors) of the modules
        modNames = substring(names(MEs), 3)
        geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
        MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
        names(geneModuleMembership) = paste("MM", modNames, sep="")
        names(MMPvalue) = paste("p.MM", modNames, sep="")
        geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"))
        GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
        names(geneTraitSignificance) = paste("GS.", names(weight), sep="")
        names(GSPvalue) = paste("p.GS.", names(weight), sep="")
        ## Intramodular analysis
        column = match(module, modNames)
        moduleGenes = moduleColors==module
        MG = data.frame(moduleGenes)
        names(MG) = module
        genes = data.frame(genes, MG)
        jpeg(paste(trait, "_", module, "gene_significance.jpg", sep=""), width=800, height=800, quality=100)
        par(mfrow = c(1,1))
        verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]), abs(geneTraitSignificance[moduleGenes, 1]), xlab = paste("Module Membership in", module, "module"), ylab = paste("Gene significance for ", trait, sep=""), main = paste("Module membership vs. gene significance\n"), pch=21, cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black", bg=module)
        dev.off()
        
        MM = data.frame(abs(geneModuleMembership[, column]))
        names(MM) = paste("MM", module, sep=".")
        genes = data.frame(genes, MM)
        
        # Outliers
        jpeg(paste(trait, "_", module, "eigengene.jpg", sep=""), width=800, height=800)
        par(mfrow=c(2,1), mar=c(0.3, 5.5, 5, 2))
        which.module=module
        plotMat(t(scale(datExpr[, mergedColors==which.module ]) ),nrgcols=30,rlabels=T, clabels=T,rcols=which.module,main=which.module,cex.main=2)
        par(mar=c(5, 4.2, 0, 0.7))
        barplot(MEs0[,paste("ME",which.module,sep="")], col=which.module, main="", cex.main=2,ylab="eigengene expression",xlab="array sample")
        dev.off()
    }
    
    
    modules.temp = gsub("ME", "MM", sign_modules)
    modules_trait = gsub("MM", paste(trait, "_MM", sep=""), modules.temp)
    temp_modules = geneModuleMembership[, modules.temp]
    names(temp_modules) = modules_trait
    sign1 = geneTraitSignificance
    sign2 = abs(geneTraitSignificance)
    
    genes = data.frame(genes, geneTraitSignificance, abs(geneTraitSignificance))
}



names(genes) = gsub(".1", "_abs", names(genes))
names(genes)[1] = "Gene_number"
head(genes, 2)

write.table(genes, "OUTPUT_gene_significance.csv", row.names=F, quote=F, sep="\t")






gs_Temperature = as.numeric(cor(allTraits[,2], datExpr, use="p"))
gs_Temperature = abs(gs_Temperature)
# Next module significance is defined as average gene significance with kruskal-wallis test.
ModuleSignificance = tapply(gs_Temperature, mergedColors, mean, na.rm=T)
#cairo_pdf(filename = "(14)gene_significance_bl.pdf")
jpeg("Temperature_module_significance.jpg", width=800, height=800)
plotModuleSignificance(gs_Temperature,mergedColors,ylim=c(0,0.8), main = "Gene significance accross module (bl)", las=3)
dev.off()






gs_Ni = as.numeric(cor(allTraits[,3], datExpr, use="p"))
gs_Ni = abs(gs_Ni)
# Next module significance is defined as average gene significance with kruskal-wallis test.
ModuleSignificance = tapply(gs_Ni, mergedColors, mean, na.rm=T)
#cairo_pdf(filename = "(14)gene_significance_bl.pdf")
jpeg("Ni_module_significance.jpg", width=800, height=800)
plotModuleSignificance(gs_Ni,mergedColors,ylim=c(0,0.8), main = "Gene significance accross module (bl)", las=3)
dev.off()




gs_Cd = as.numeric(cor(allTraits[,4], datExpr, use="p"))
gs_Cd = abs(gs_Cd)
# Next module significance is defined as average gene significance with kruskal-wallis test.
ModuleSignificance = tapply(gs_Cd, mergedColors, mean, na.rm=T)
#cairo_pdf(filename = "(14)gene_significance_bl.pdf")
jpeg("Cd_module_significance.jpg", width=800, height=800)
plotModuleSignificance(gs_Cd,mergedColors,ylim=c(0,0.8), main = "Gene significance accross module (bl)", las=3)
dev.off()






# Calculate connectivity
GeneConnectivity=softConnectivity(datExpr, power=power, verbose = 3)
GeneConnectivity = GeneConnectivity/max(GeneConnectivity)
#sort(GeneConnectivity)

# Calculate gene significance for phenotype
Temperature_sign = standardScreeningNumericTrait(datExpr, datTraits[,1], alternative = "two.sided")
Ni_sign = standardScreeningNumericTrait(datExpr, datTraits[,2], alternative = "two.sided")
Cd_sign = standardScreeningNumericTrait(datExpr, datTraits[,3], alternative = "two.sided")


gene_choice = data.frame(GeneConnectivity, Temperature_sign$pvalueStudent, Ni_sign$pvalueStudent, Cd_sign$pvalueStudent)

plot(gene_choice$GeneConnectivity, gene_choice$Temperature_sign.pvalueStudent)
plot(gene_choice$GeneConnectivity, gene_choice$Ni_sign.pvalueStudent)
plot(gene_choice$GeneConnectivity, gene_choice$Cd_sign.pvalueStudent)





### Gene-trait association
## Define variable weight containing the weight column of datTrait
## traits_modules = list(c("Ni", "lightcyan"), c("Cd", "red"), c("Cd", "green"), c("Cd", "pink"))
#trait = "Temperature"
#module = "black"
#weight = as.data.frame(datTraits[,trait])
#names(weight) = trait

### Significant modules
## **red = Cd ~14, p=4.7e-08
## *pink = Cd ~4, p=0.03
## green = Cd ~3, p=0.42
## *lightcyan = Ni ~8 genes, p=0.066

## names (colors) of the modules
#modNames = substring(names(MEs), 3)
#geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
#MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
#names(geneModuleMembership) = paste("MM", modNames, sep="")
#names(MMPvalue) = paste("p.MM", modNames, sep="")
#geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"))
#GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
#names(geneTraitSignificance) = paste("GS.", names(weight), sep="")
#names(GSPvalue) = paste("p.GS.", names(weight), sep="")

### Intramodular analysis
#column = match(module, modNames)
#moduleGenes = moduleColors==module
#sizeGrWindow(10.25, 5.75)
#verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
#     abs(geneTraitSignificance[moduleGenes, 1]),
#     xlab = paste("Module Membership in", module, "module"),
#     ylab = paste("Gene significance for", names(weight)),
#     main = paste("Module membership vs. gene significance\n"),
#     pch=21, cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black", bg=module)








### Summary output
##names(datExpr)
##names(datExpr)[moduleColors=="darkmagenta"]

#annot = names(datExpr)#read.csv(file = "GeneAnnotation.csv")
#dim(annot)
#names(annot)
#probes = annot#names(datExpr)
#probes2annot = annot#match(probes, annot$substanceBXH)
## The following is the number or probes without annotation:
#sum(is.na(probes2annot))
## Should return 0.

## Create the starting data frame
####OLD### geneInfo0 = data.frame(substanceBXH = probes, geneSymbol = annot$gene_symbol[probes2annot], LocusLinkID = annot$LocusLinkID[probes2annot], moduleColor = moduleColors, geneTraitSignificance, GSPvalue)
#geneInfo0 = data.frame(substanceBXH = probes, geneSymbol = annot, LocusLinkID = annot, moduleColor = moduleColors, geneTraitSignificance, GSPvalue)
## Order modules by their significance for weight
#modOrder = order(-abs(cor(MEs, weight, use = "p")))
## Add module membership information in the chosen order
#for (mod in 1:ncol(geneModuleMembership))
#{
#    oldNames = names(geneInfo0)
#    geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],     MMPvalue[, modOrder[mod]])
#    names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),     paste("p.MM.", modNames[modOrder[mod]], sep=""))
#}
## Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
#geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight))
#geneInfo = geneInfo0[geneOrder, ]

#write.csv(geneInfo, file = "perch-mercury-geneInfo.csv")

##### END OF MODULE TO TRAIT RELATIONSHIP ###







#### START OF GO INTERFACING ###
## Tutorial file 4
#### END OF GO INTERFACING ###








#### START OF NETWORK VISUALISATION ###
## tutorial file 5

### Preliminaries
## Load the WGCNA package
#library(WGCNA)
## The following setting is important, do not omit.
#options(stringsAsFactors = FALSE)
## Load the expression and trait data saved in the first part
#lnames = load(file = prepared.data)
##The variable lnames contains the names of loaded variables.
#lnames
## Load network data saved in the second part.
#lnames = load(file = wgcna.modules)
#lnames
#nGenes = ncol(datExpr)
#nSamples = nrow(datExpr)

## Calculate topological overlap anew: this could be done more efficiently by saving the TOM
## calculated during module detection, but let us do it again here.
#dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 7)
## Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
#plotTOM = dissTOM^7
## Set diagonal to NA for a nicer plot
#diag(plotTOM) = NA
## Call the plot function
##sizeGrWindow(9,9)
##TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

### Restrict the number of plotted genes
#nSelect = 200
## For reproducibility, we set the random seed
##set.seed(10)
#select = sample(nGenes, size = nSelect)
#selectTOM = dissTOM[select, select]
## There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
#selectTree = flashClust(as.dist(selectTOM), method = "average")
#selectColors = moduleColors[select]
## Open a graphical window
#sizeGrWindow(10, 10)
## Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
## the color palette; setting the diagonal to NA also improves the clarity of the plot
#plotDiss = selectTOM^7
#diag(plotDiss) = NA
#TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")

### Visualize the network eigengenes
## Recalculate module eigengenes
#MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
## Isolate weight from the clinical traits
#weight = as.data.frame(datTraits$log_Inorghg)
#names(weight) = "log_Inorghg"
## Add the weight to existing module eigengenes
#MET = orderMEs(cbind(MEs, weight))
## Plot the relationships among the eigengenes and the trait
##sizeGrWindow(10.25, 5.75)
#par(cex = 0.9)
#plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle = 90)

## Plot the dendrogram
#sizeGrWindow(10.25, 5.75)
#par(cex = 1.0)
#plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0), plotHeatmaps = FALSE)
## Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
#par(cex = 1.0)
#plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)

#### END OF NETWORK VISUALISATION ###

#### START OF NETWORK EXPORT ###
## Tutorial file 6
#### END OF NETWORK EXPORT ###





































