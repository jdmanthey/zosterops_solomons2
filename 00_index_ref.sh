interactive -p nocona -c 4 -m 4G

cd references

source activate bcftools

samtools faidx GCA_965231275.1_bZosLat1.hap1.1_genomic.fna

java -jar picard.jar CreateSequenceDictionary \
R=/home/jmanthey/references/GCA_965231275.1_bZosLat1.hap1.1_genomic.fna \
O=/home/jmanthey/references/GCA_965231275.1_bZosLat1.hap1.1_genomic.dict

bwa-mem2 index GCA_965231275.1_bZosLat1.hap1.1_genomic.fna
