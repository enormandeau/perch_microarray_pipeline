###########################################################################################################
# Code adapted by Eric Normandeau from "Rcode_microarray_analysis_vIV.txt", originaly by Sebastien Renault
#
# What this code does:
#
# - Use gene expression data from a R-maanova analysis to:
#
#   - Prepare data for Heat-Map (remove unwanted effects)
#
#   Tested with R version 2.8.1 and maanova package version 1.13.3
#
# 2012 01 12
###########################################################################################################


###########################
# Import needed objects

load("OUTPUT_anova.RSession")

# Define used object names here:
design.file = "design.txt"
lowess.data = lowess.data
anova.model = anova.mix1
output.file = "OUTPUT_fitted_data_scaled.txt"
results = test.FDR

# Import design file
design = read.table(design.file, header = T)


################################
# Prepare data for Heat-Map


# Import design file
design = read.table(design.file, header = T)

# Fitted data (yhat) for significant genes
data.fitted = lowess.data$data

# Remove array effect

array.levels = levels(as.factor(design$Array))
nb.arrays = length(array.levels)

for(i in 1:nb.arrays)
{
    data.fitted[ , design$Array == array.levels[i]] = data.fitted[ , design$Array == array.levels[i]] - anova.model$Array[,i]
}


# Remove dye effect

dye.levels = levels(as.factor(design$Dye))
nb.dyes = length(dye.levels)

for(i in 1:nb.dyes)
{
    data.fitted[ , design$Dye == dye.levels[i]] = data.fitted[ , design$Dye == dye.levels[i]] - anova.model$Dye[,i]
}

# Remove other technical effects in the same way


#data.fitted = scale(data.fitted)
#boxplot(data.fitted)

write.table(data.fitted, output.file, sep="\t", col.names=F, row.names=F)


