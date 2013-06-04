# Microarray Analysis program - Version 3.3
#
# Eric Normandeau
# 2011-07-26

#rm(list = ls())
#memory.size(4000)
library(maanova)


### Global Variables

input.data = "OUTPUT_1c_data.txt"
design.file = "design.txt"
result.file = "OUTPUT_lakes_spring_2011_results.txt"
fold.change.file = "OUTPUT_anova_lakes_spring_2011_fold-changes.txt"
r.session.file = "OUTPUT_anova_lakes_spring_2011.RSession"

nb.perm = 2
my.model = ~Array + Dye + Region
my.random = ~Array
term.tested = "Region"

### Importing the data (in tab delimited file)

data = read.madata(input.data,
    arrayType='twoColor', log.trans=T, designfile=design.file, spotflag=F,
    probeid=1, metarow=2, metacol=3, row=4, col=5, intensity=6)

lowess.data = transform.madata (data, method="rlowess", f=0.5, degree=2,
    iter=3, draw="off")

print(dim(lowess.data$data))

write.table(lowess.data$data, "OUTPUT_lowess_data.txt", sep="\t", col.names=F, row.names=F)

##################################
# Creating and fitting ANOVA model

#anova.mix1 <- fitmaanova(lowess.data, formula = my.model,
#                         random = my.random)

######################################
# ANOVA testing and permutations
# May take a ***VERY*** long time to run!!!

### Statistical test (F tests, permutations)

#test.mix.Strain <- matest(lowess.data, anova.mix1, term=term.tested, n.perm=nb.perm,
#    shuffle.method="sample", pval.pool=T, verbose=T)


#### Adjusting p-values with FDR (False Discovery Rate)
#x <- test.mix.Strain
#test.FDR <- adjPval(x, method = "jsFDR")


################
## Volcano Plots for anova
##library(maanova)

#### Using tabulated p-values
#idx.mix.tab <- volcano(test.FDR, threshold=c(0.05,0.05), method=c("unadj","unadj"),
#    title="Tabulated Values", highlight.flag=F, onScreen=F)

#### Using permutation pvalues
#idx.mix.perm <- volcano(test.FDR, threshold=c(0.05,0.05), method=c("nominal","nominal"),
#    title="Permutated Values", highlight.flag=F, onScreen=F)

#### FDR adjusted permutation pvalues
#idx.mix.fdr <- volcano(test.FDR, threshold=c(0.05,0.05), method=c("fdrperm", "fdrperm"),
#    title="FDR corrected Permutated Values", highlight.flag=F, onScreen=F)

#### List of significant genes according to :
## Tabulated, Permutated and FDR corrected p-values
#summary(idx.mix.tab)
#summary(idx.mix.perm)
#summary(idx.mix.fdr)

#hist(test.FDR$Fs$Ptab,seq(0,1,0.05))

#write.table(test.FDR$Fs, result.file, col.names=T, 
#    row.names=F)

#save.image(r.session.file)

##library(maanova)
#fold.changes = calVolcanoXval(test.FDR)
#write.table(fold.changes, fold.change.file, row.names=F, col.names=F, quote=F)


