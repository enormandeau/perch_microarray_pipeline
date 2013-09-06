# perch_microarray_pipeline
# 1b Importing in maanova to impute missing data

### Global variables
input.file = "OUTPUT_1a_data.txt"
output.file = "OUTPUT_1b_data.txt"
spot.id.file = "spot_ids.txt"
design.file = "files_to_fill/design.txt"

# Importing libraries
library(maanova)

# Read data into maanova to impute missing data
data = read.madata(input.file,
    arrayType='twoColor', log.trans=F, designfile=design.file, spotflag=T,
    probeid=1, metarow=2, metacol=3, row=4, col=5, intensity=6)

# Write output
write.table(data$data, output.file, sep="\t", col.names=F,
    row.names=F, quote=F)

