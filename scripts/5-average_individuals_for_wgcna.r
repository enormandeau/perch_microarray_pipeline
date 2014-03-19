rm(list=ls())

data = read.table("OUTPUT_lowess_data_for_averages.txt", header=T)
sample.names = read.table("sample_names.txt", stringsAsFactors=F)
sample.names = as.vector(t(sample.names))

averaged = data.frame(temp=seq(1008))

for (name in unique(sample.names)) {
    d = data.frame(data[,sample.names == name])
    if (ncol(d) > 1) {
        ok = data.frame(apply(d, 1, mean))
    } else {
        ok = d
    }
    names(ok) = name
    averaged = data.frame(averaged, ok)
}

averaged = averaged[, -1]

write.table(averaged, "OUTPUT_lowess_data_for_wgcna_averaged.txt", sep="\t", quote=F, row.names=F)

