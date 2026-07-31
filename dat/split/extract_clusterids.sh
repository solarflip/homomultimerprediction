#!/bin/bash

FILE=$1 #datasets/homomer_pdbids_hash_clusterid_labels_fullset.csv
OUTPUT=homomer_cluster_fullset.txt

# extract rcsb ids
awk -F',' '{printf "%s%s", NR==1?"":",", $5} END{print ""}' "$FILE" > "$OUTPUT"