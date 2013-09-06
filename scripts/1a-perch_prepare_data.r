# perch_microarray_pipeline
# 1a Removing the background signal

### Global variables
input.file = "raw_data.txt"
output.file = "OUTPUT_1a_data.txt"
spot.id.file = "spot_ids.txt"
header.file = "header.txt"
filenames.file = "files_to_fill/filenames.txt"

# Calculate number of slides
filenames = read.table(filenames.file)
NUMBER.SLIDES = nrow(filenames)

# Importing libraries
library(maanova)

# Read data
data = read.table(input.file, sep="\t", header=T, fill=T)
dim(data)


# Subtracting backgrounds
data.sub=NULL
for (i in 1:NUMBER.SLIDES) {
    data.sub=cbind(data.sub, data[,5*(i-1) + 1] - data[,5*(i-1) + 2],
        data[,5*(i-1) + 3] - data[,5*(i-1) + 4], data[,5*(i-1) + 5])
}
dim(data.sub)

# Trick because one section is not needed in this script
data.crit = data.sub 

# Add 5 columns at beginning of file with :
# probeid, metarow, metacol, row, col
gene.location = read.table(spot.id.file, header=T, sep="\t")
data.crit=cbind(gene.location, data.crit)
print(dim(data.crit))

header = names(read.table(header.file, header=T))

names(data.crit) = header[1:ncol(data.crit)]

# Write output
write.table(data.crit, output.file, sep="\t", col.names=T,
    row.names=F, quote=F)

