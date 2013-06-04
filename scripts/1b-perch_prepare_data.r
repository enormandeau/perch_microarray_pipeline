# Microarray Analysis program - Version 3.3
# Preparing data for maanova
#
# Eric Normandeau
# 2011-07-26


### Global variables

input.file = "OUTPUT_1a_data.txt"
output.file = "OUTPUT_1b_data.txt"
spot.id.file = "spot_ids.txt"
design.file = "design.txt"

### Load maanova package
library(maanova)

data = read.madata(input.file,
    arrayType='twoColor', log.trans=F, designfile=design.file, spotflag=T,
    probeid=1, metarow=2, metacol=3, row=4, col=5, intensity=6)

write.table(data$data, output.file, sep="\t", col.names=F,
    row.names=F, quote=F)

