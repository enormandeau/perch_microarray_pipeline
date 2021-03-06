# perch_microarray_pipeline
# 1c Averaging the 6 spot copies for each probe

### Global variables
input.file = "OUTPUT_1b_data.txt"
output.file = "OUTPUT_1c_data.txt"
design.file = "files_to_fill/design.txt"
header.file = "header.txt"
spot.id.file = "spot_ids.txt"

# Importing libraries
library(maanova)

# Open inputed data
data = read.table(input.file, sep="\t", header=F, fill=T)
dim(data)

# Add 5 columns at beginning of file with :
# probeid, metarow, metacol, row, col
gene.location = read.table(spot.id.file, header=T, sep="\t")
data.crit=cbind(gene.location, data)  
print(dim(data.crit))

# Average the 3 sub-arrays
sub1 = data.crit[1:2016, ]
sub2 = data.crit[2017:4032, ]
sub3 = data.crit[4033:6048, ]

data.mean = (sub1 + sub2 + sub3) / 3
data.mean[, 1:5] = sub1[, 1:5]

dim(data.mean)
head(data.mean, 12)

# Average the 2 replicates for each spot
rep1 = data.mean[data.mean[,1] %% 2 != 0, ]
rep2 = data.mean[data.mean[,1] %% 2 == 0, ]

data.rep = (rep1 + rep2) / 2
data.rep[, 1:5] = rep1[, 1:5]

dim(data.rep)
head(data.rep, 12)


# Visualizing quality of first 4 arrays
par(mfrow=c(4,2))

plot(log(data.rep[,6], 2),  log(data.rep[,7], 2))
plot(log(data.rep[,6], 2) + log(data.rep[,7], 2),
     log(data.rep[,6], 2) - log(data.rep[,7], 2), ylim=c(-2,2))

plot(log(data.rep[,9], 2),  log(data.rep[,10], 2))
plot(log(data.rep[,9], 2) + log(data.rep[,10], 2),
     log(data.rep[,9], 2) - log(data.rep[,10], 2), ylim=c(-2,2))

plot(log(data.rep[,12], 2),  log(data.rep[,13], 2))
plot(log(data.rep[,12], 2) + log(data.rep[,13], 2),
     log(data.rep[,12], 2) - log(data.rep[,13], 2), ylim=c(-2,2))

plot(log(data.rep[,15], 2),  log(data.rep[,16], 2))
plot(log(data.rep[,15], 2) + log(data.rep[,16], 2),
     log(data.rep[,15], 2) - log(data.rep[,16], 2), ylim=c(-2,2))


# Output file for Maanova
header = names(read.table(header.file, header=T))
names(data.rep) = header[1:ncol(data.rep)]

write.table(data.rep, output.file, sep="\t", col.names=T,
    row.names=F, quote=F)

