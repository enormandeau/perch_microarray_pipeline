# perch_microarray_pipeline
# 4 extract fitted data

# Import needed objects
load("OUTPUT_anova.RSession")

# Define used object names here:
design.file = "files_to_fill/design.txt"
lowess.data = lowess.data
anova.model = anova.mix1
output.file = "OUTPUT_fitted_data.txt"

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

write.table(data.fitted, output.file, sep="\t", col.names=F, row.names=F)


