cd /lustre/scratch/jmanthey/02_zosterops/12_stats

# combine the output for different windows into a single file 
# first add a header for each file
grep 'pop1' OZ246474.1__10000001__10100000__stats.txt > ../window_100kbp_stats.txt

# add the stats to the file
for i in $( ls *__stats.txt ); do grep -v 'pop1' $i >> ../window_100kbp_stats.txt; done
