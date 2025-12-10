cat *stats >> ../genotyping.stats

# in R:
R
options(scipen=999)

n_inds <- 47
genome_size <- read.table("GCA_965231275.1_bZosLat1.hap1.1_genomic.fna.fai")
genome_size <- sum(genome_size[,2])
x <- scan("genotyping.stats", what="character", sep="\n")

taxon <- paste0(sapply(strsplit(x[seq(1:n_inds) * 7 - 6], "_"), "[[", 1), "_", sapply(strsplit(x[seq(1:n_inds) * 7 - 6], "_"), "[[", 2))
locality <- sapply(strsplit(x[seq(1:n_inds) * 7 - 6], "_"), "[[", 3)
id <- sapply(strsplit(x[seq(1:n_inds) * 7 - 6], "_"), "[[", 4)
depth <- as.numeric(sapply(strsplit(x[seq(1:n_inds) * 7 - 4], "=  "), "[[", 2)) / genome_size
prop_dupes <- as.numeric(x[seq(1:n_inds) * 7 - 2])
sites <- as.numeric(x[seq(1:n_inds) * 7 - 0])

output <- data.frame(taxon=as.character(taxon), locality=as.character(locality), id=as.character(id), depth= as.numeric(depth), prop_dupes=as.numeric(prop_dupes), sites= as.numeric(sites))

write.table(output, file="zosterops_genotyping_stats.txt", sep="\t", row.names=F, quote=F)
