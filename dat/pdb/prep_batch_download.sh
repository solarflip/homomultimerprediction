#!/bin/bash

FILE=datasets/homomer_pdbids_hash_clusterid_labels_fullset.csv
OUTPUT=homomer_pdbids_fullset.txt

# extract rcsb ids
awk -F'[,_]' '{printf "%s%s", NR==1?"":",", $1} END{print ""}' "$FILE" > $OUTPUT