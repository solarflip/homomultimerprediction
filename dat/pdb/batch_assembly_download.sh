#!/bin/bash

input_file=dat/pdb/homomer_pdbids_fullset.txt

# Read the file, split on commas, and iterate
IFS=',' read -ra items < "$input_file"

for item in "${items[@]}"; do
    echo "Processing: $item"

    python dat/pdb/download_assemblies.py $item -o dat/pdb/assemblies
done