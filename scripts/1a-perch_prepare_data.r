# Microarray Analysis program - Version 3.3
# Preparing data for maanova
#
# Eric Normandeau
# 2011-07-26


### Global variables
input.file = "raw_data.txt"
output.file = "OUTPUT_1a_data.txt"
spot.id.file = "spot_ids.txt"
header.file = "header.txt"
filenames.file = "filenames.txt"

# Calculate number of slides
filenames = read.table(filenames.file)
NUMBER.SLIDES = nrow(filenames)

### Load maanova package
library(maanova)


#########################################
# Open the ScanArray files
# Put all the expression values, 
# background and flags in the same matrix

data = read.table(input.file, sep="\t", header=T, fill=T)
dim(data)


#########################
# Subtracting backgrounds

data.sub=NULL
for (i in 1:NUMBER.SLIDES) {
    data.sub=cbind(data.sub, data[,5*(i-1) + 1] - data[,5*(i-1) + 2],
        data[,5*(i-1) + 3] - data[,5*(i-1) + 4], data[,5*(i-1) + 5])
}
dim(data.sub)


#############################
# Replace flagged genes by NA

data.crit = data.sub # Trick because one section is not needed in this script

#for (i in 1:(NUMBER.SLIDES)) {
#    data.crit[data.crit[,3*(i-1)+3] == 0, 3*(i-1) + 1] = NA
#    data.crit[data.crit[,3*(i-1)+3] == 0, 3*(i-1) + 2] = NA
#}
#print(dim(data.crit))


###########################################
# Add 5 columns at beginning of file with :
# probeid, metarow, metacol, row, col

gene.location = read.table(spot.id.file, header=T, sep="\t")
data.crit=cbind(gene.location, data.crit)
print(dim(data.crit))

header = names(read.table(header.file, header=T))

names(data.crit) = header[1:ncol(data.crit)]

write.table(data.crit, output.file, sep="\t", col.names=T,
    row.names=F, quote=F)

cat(output.file, "\n")

