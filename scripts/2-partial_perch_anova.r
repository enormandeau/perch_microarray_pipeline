# perch_microarray_pipeline
# 2 partial - First half of 2-perch_anova.r to use lowess normalization

# Importing libraries
library(maanova)

# Global Variables
input.data = "OUTPUT_1c_data.txt"
design.file = "design.txt"
result.file = "OUTPUT_lowess_data.txt"

# Importing the data (in tab delimited file)
data = read.madata(input.data,
    arrayType='twoColor', log.trans=T, designfile=design.file, spotflag=F,
    probeid=1, metarow=2, metacol=3, row=4, col=5, intensity=6)

lowess.data = transform.madata (data, method="rlowess", f=0.5, degree=2,
    iter=3, draw="off")

write.table(lowess.data$data, result.file, sep="\t", col.names=F, row.names=F)

