
options(scipen=999)

# list all the files in the trees directory
x_files <- list.files("windows", pattern="RAxML_bipartition*", full.names=T)

# find the chromosome, start, and end for each tree
x_names <- list.files("windows", pattern="RAxML_bipartition*")
x_chrom <- sapply(strsplit(sapply(strsplit(x_names, "RAxML_bipartitions."), "[[", 2), "__"), "[[", 1)
x_start <- sapply(strsplit(x_names, "__"), "[[", 2)
x_end <- sapply(strsplit(sapply(strsplit(x_names, "__"), "[[", 3), ".tre"), "[[", 1)

output <- data.frame(chrom=as.character(x_chrom), start=as.numeric(x_start), end=as.numeric(x_end))

# write tree info
write.table(output, file="zosterops_100kbp_tree_info.txt", sep="\t", quote=F, row.names=F, col.names=F)

# trees into one file
tree_list <- list()
for(a in 1:length(x_files)) {
	tree_list[[a]] <- scan(x_files[a], what="character")
}
# make sure only one tree per file (in case a tree was inferred twice in the array jobs)
for(a in 1:length(tree_list)) {
	if(length(tree_list[[a]]) > 1) {
	tree_list[[a]] <- tree_list[[a]][1]
	}
}
tree_list <- unlist(tree_list)
write(tree_list, file="zosterops_100kbp.trees", ncolumns=1)
