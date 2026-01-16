options(scipen=999)

# input args <- input file, popmap
args <- commandArgs(trailingOnly = TRUE)

# base file name for parsing chromosome and window bounds
filename_simple <- args[1]
filename_simple <- strsplit(filename_simple, "windows/")[[1]][2]
filename_simple <- strsplit(filename_simple, ".recode")[[1]][1]
chromosome_info <- c(strsplit(filename_simple, "__")[[1]][1], 
                    as.numeric(strsplit(strsplit(filename_simple, "__")[[1]][2], "__")[[1]][1]),
                    as.numeric(strsplit(filename_simple, "__")[[1]][3]))


# determine base outname (edit as needed)
outname_base <- "filename_simple"

# define minimum number of sites to keep a file 
min_sites <- 10000

############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
# do not modify below this point
############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################


####################################
# read in a small vcf (don't use 
# for large vcf files)
####################################
read_vcf <- function(input_file) {
  header <- readLines(input_file)
  header <- header[grep('^#C', header)]
  header <- strsplit(header, "\t")[[1]]
  vcf <- read.table(input_file, header=F)
  colnames(vcf) <- header
  return(vcf)
}



####################################
# make MrBayes Nexus from VCF
# heterozygous sites choose one allele
# at random
####################################
create_nexus_from_vcf <- function(xxx, outname, num_sites_needed) {
	# keep going only if enough sites
	if(num_sites_needed <= nrow(xxx)) {
		
		# subset the genotypes from the allele info
		allele_info <- xxx[,c(4,5)]
		genotypes <- xxx[,10:ncol(xxx)]
		# define names of individuals
		output_names_nexus <- colnames(genotypes)
		
		# define heterozygous ambiguities
    		het <- rep("N", nrow(allele_info))
    		het[allele_info[,1] == "A" & allele_info[,2] == "C"] <- "M"
     	het[allele_info[,1] == "A" & allele_info[,2] == "G"] <- "R"
      	het[allele_info[,1] == "A" & allele_info[,2] == "T"] <- "W"
       	het[allele_info[,1] == "C" & allele_info[,2] == "A"] <- "M"
       	het[allele_info[,1] == "C" & allele_info[,2] == "G"] <- "S"
       	het[allele_info[,1] == "C" & allele_info[,2] == "T"] <- "Y"
       	het[allele_info[,1] == "G" & allele_info[,2] == "A"] <- "R"
       	het[allele_info[,1] == "G" & allele_info[,2] == "C"] <- "S"
       	het[allele_info[,1] == "G" & allele_info[,2] == "T"] <- "K"
       	het[allele_info[,1] == "T" & allele_info[,2] == "A"] <- "W"
       	het[allele_info[,1] == "T" & allele_info[,2] == "C"] <- "Y"
       	het[allele_info[,1] == "T" & allele_info[,2] == "G"] <- "K"
       	
       	# convert all numbers in genotypes to actual bases and ambiguities		
		for(a in 1:ncol(genotypes)) {
			# keep only genotype
			genotypes[,a] <- substr(genotypes[,a], 1, 3)
			# convert to bases
			genotypes[,a][genotypes[,a] == "0/0"] <- allele_info[genotypes[,a] == "0/0",1]
			genotypes[,a][genotypes[,a] == "1/1"] <- allele_info[genotypes[,a] == "1/1",2]
			genotypes[,a][genotypes[,a] == "./."] <- "?"
			genotypes[,a][genotypes[,a] == "0/1"] <- het[genotypes[,a] == "0/1"]
			# choose random allele for heterozygotes
			genotypes[,a][genotypes[,a] == "R"] <- sample(c("A", "G"), length(genotypes[,a][genotypes[,a] == "R"]), replace=T)
			genotypes[,a][genotypes[,a] == "K"] <- sample(c("T", "G"), length(genotypes[,a][genotypes[,a] == "K"]), replace=T)
			genotypes[,a][genotypes[,a] == "S"] <- sample(c("C", "G"), length(genotypes[,a][genotypes[,a] == "S"]), replace=T)
			genotypes[,a][genotypes[,a] == "Y"] <- sample(c("C", "T"), length(genotypes[,a][genotypes[,a] == "Y"]), replace=T)
			genotypes[,a][genotypes[,a] == "M"] <- sample(c("A", "C"), length(genotypes[,a][genotypes[,a] == "M"]), replace=T)
			genotypes[,a][genotypes[,a] == "W"] <- sample(c("A", "T"), length(genotypes[,a][genotypes[,a] == "W"]), replace=T)

		}
		
		n_inds <- ncol(genotypes)
		n_snps <- nrow(genotypes)
		
		to_write <- "#NEXUS"
		write(to_write, file=outname, ncolumns=1)
		to_write <- ""
		write(to_write, file=outname, ncolumns=1, append=T)
		to_write <- "begin data;"
		write(to_write, file=outname, ncolumns=1, append=T)
		to_write <- paste0("\tdimensions ntax=", n_inds, " nchar=", n_snps, ";")
		write(to_write, file=outname, ncolumns=1, append=T)
		to_write <- "\tformat datatype=DNA interleave=no missing=?;"
		write(to_write, file=outname, ncolumns=1, append=T)
		to_write <- "\tmatrix"
		write(to_write, file=outname, ncolumns=1, append=T)
		
		# write info for each individual
		for(a in 1:n_inds) {
			to_write <- paste0("\t\t", output_names_nexus[a], "\t", paste(genotypes[,a], collapse=""))
			write(to_write, file=outname, ncolumns=1, append=T)
		}
		to_write <- "\t;"
		write(to_write, file=outname, ncolumns=1, append=T)
		to_write <- "end;"
		write(to_write, file=outname, ncolumns=1, append=T)
		
		# mrbayes block
		to_write <- "begin mrbayes;"
		write(to_write, file=outname, ncolumns=1, append=T)
		to_write <- "\tset autoclose=yes nowarnings=yes;"
		write(to_write, file=outname, ncolumns=1, append=T)
		to_write <- "\tlset nst=6 rates=invgamma;"
		write(to_write, file=outname, ncolumns=1, append=T)
		to_write <- "\tmcmcp nchains=2 nruns=1 ngen=500000 samplefreq=2000 printfreq=10000 diagnfreq=10000 temp=0.4 swapfreq=10 burninfrac=0.25;"
		write(to_write, file=outname, ncolumns=1, append=T)
		to_write <- "\tmcmc;"
		write(to_write, file=outname, ncolumns=1, append=T)
		to_write <- "\tsumt;"
		write(to_write, file=outname, ncolumns=1, append=T)
		to_write <- "end;"
		write(to_write, file=outname, ncolumns=1, append=T)


	}
}



####################################
# make MrBayes file
####################################

x <- read_vcf(args[1])

# define output nexus name
output_nexus <- paste(strsplit(args[1], ".recode.vcf")[[1]][1], ".nex", sep="")
# create fasta sequence alignments for all files with a minimum number of sites
create_nexus_from_vcf(x, output_nexus, min_sites)







