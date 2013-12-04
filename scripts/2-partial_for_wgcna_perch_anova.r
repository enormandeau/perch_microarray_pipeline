# perch_microarray_pipeline
# 2 partial - First half of 2-perch_anova.r to use lowess normalization

# Importing libraries
library(maanova)

# Global Variables
input.data = "OUTPUT_1c_data.txt"
design.file = "files_to_fill/design.txt"
result.file = "OUTPUT_lowess_data_for_wgcna.txt"

# Importing the data (in tab delimited file)
data = read.madata(input.data,
    arrayType='twoColor', log.trans=T, designfile=design.file, spotflag=F,
    probeid=1, metarow=2, metacol=3, row=4, col=5, intensity=6)

lowess.data = transform.madata (data, method="rlowess", f=0.5, degree=2,
    iter=3, draw="off")

data.wgcna = apply(lowess.data$data, 2, scale)

write.table(data.wgcna, result.file, sep="\t", col.names=F, row.names=F)

