data = read.table("OUTPUT_gene_significance.csv", header=T, stringsAsFactors=F)

module.names = names(data)
wanted = module.names[substr(module.names, 1, 1) %in% letters]
gs.names = module.names[substr(module.names, 1, 3) == "GS." &
    substr(module.names, nchar(module.names) - 3, nchar(module.names)) == "_abs"]
gs.cols = which(module.names %in% gs.names)
gene.numbers = data[,1]
results = data.frame(gene.numbers)

for (i in wanted) {
    col.name = which(module.names == i)
    col.mm = col.name + 1
    col.gs = gs.cols[gs.cols > col.mm][1]
    #cat("Wanted: ", i, "\n    ", col.name, " ", col.mm, " ", col.gs, "\n")
    good.genes = data[data[,col.name] == "TRUE" &
                      data[,col.mm] >= 0.6 &
                      data[,col.gs] >= 0.2, 1]
    flag = rep(0, 1008)
    flag[good.genes] = 1
    results = data.frame(results, flag)
    names(results)[length(names(results))] = paste0(names(data)[col.name], "_", gsub("_abs", "", gsub("GS.", "", names(data)[col.gs])))
}

write.table(results, "OUTPUT_gene_significance_flags.csv", quote=F, row.names=F, sep="\t")

