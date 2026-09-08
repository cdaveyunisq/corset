#!/bin/bash

echo corset_par check number of clusters with leq 1 read.
awk 'NR>1 {sum=0; for(i=2; i<=NF; i++) sum+=$i; if(sum <= 1) count++} END {print "Clusters with 0 or 1 total reads: " count}' examples/Cx_sitiens_par/Cx_sitiens_par-counts.txt

echo corset_par1 check number of clusters with leq 1 read.
awk 'NR>1 {sum=0; for(i=2; i<=NF; i++) sum+=$i; if(sum <= 1) count++} END {print "Clusters with 0 or 1 total reads: " count}' examples/Cx_sitiens_par1/Cx_sitiens_par1-counts.txt


echo corset_sync check number of clusters with leq 1 read.
awk 'NR>1 {sum=0; for(i=2; i<=NF; i++) sum+=$i; if(sum <= 1) count++} END {print "Clusters with 0 or 1 total reads: " count}' examples/Cx_sitiens_sync1/Cx_sitiens_sync-counts.txt


# visual comparison of counts distributions
echo corset_par CustomDSU counts summary
awk 'NR>1 {sum=0; for(i=2; i<=NF; i++) sum+=$i; print int(log(sum+1)/log(2))}' examples/Cx_sitiens_par/Cx_sitiens_par-counts.txt | sort -n | uniq -c

echo corset_par DSU counts summary
awk 'NR>1 {sum=0; for(i=2; i<=NF; i++) sum+=$i; print int(log(sum+1)/log(2))}' examples/Cx_sitiens_par1/Cx_sitiens_par1-counts.txt | sort -n | uniq -c


echo corset_sync counts summary
awk 'NR>1 {sum=0; for(i=2; i<=NF; i++) sum+=$i; print int(log(sum+1)/log(2))}' examples/Cx_sitiens_sync1/Cx_sitiens_sync-counts.txt | sort -n | uniq -c



echo "--------------------------------------------"


echo corset_par CustomDSU to 20 cluster sizes

awk '{print $2}' examples/Cx_sitiens_par/Cx_sitiens_par-clusters.txt | sort | uniq -c | sort -nr | head -n 20


echo corset_par DSU to 20 cluster sizes

awk '{print $2}' examples/Cx_sitiens_par1/Cx_sitiens_par1-clusters.txt | sort | uniq -c | sort -nr | head -n 20


echo corset_sync to 20 cluster sizes
awk '{print $2}' examples/Cx_sitiens_sync1/Cx_sitiens_sync-clusters.txt | sort | uniq -c | sort -nr | head -n 20