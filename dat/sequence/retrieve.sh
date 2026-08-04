#!/bin/bash

FILE=dat/datasets/homomer_pdbids_hash_clusterid_labels_fullset.csv
OUTPUT_DIR=dat/sequence/fullset

mkdir -p "$OUTPUT_DIR"

while IFS=, read -r id _; do
    [[ -z "$id" ]] && continue

    id=${id^^}

    pdb_id=${id%%_*}
    chain_id=${id##*_}

    echo "Processing ${pdb_id} chain ${chain_id}..."

    # Get polymer entity IDs and their author chains
    entity_id=$(curl -fsSL \
        "https://data.rcsb.org/rest/v1/core/polymer_entity/${pdb_id}/1" |
        jq -r --arg chain "$chain_id" '
            if .rcsb_polymer_entity_container_identifiers.auth_asym_ids
            | index($chain)
            then
                .rcsb_polymer_entity_container_identifiers.entity_id
            else
                empty
            end
        ')

    # If entity 1 was not the correct one, search all entities
    if [[ -z "$entity_id" ]]; then
        max_entity=$(curl -fsSL \
            "https://data.rcsb.org/rest/v1/core/entry/${pdb_id}" |
            jq -r '.rcsb_entry_container_identifiers.polymer_entity_ids[]' |
            tail -n 1)

        for e in $(seq 1 "$max_entity"); do
            entity_id=$(curl -fsSL \
                "https://data.rcsb.org/rest/v1/core/polymer_entity/${pdb_id}/${e}" |
                jq -r --arg chain "$chain_id" '
                    select(
                        .rcsb_polymer_entity_container_identifiers.auth_asym_ids
                        | index($chain)
                    )
                    | .rcsb_polymer_entity_container_identifiers.entity_id
                ')

            [[ -n "$entity_id" && "$entity_id" != "null" ]] && break
        done
    fi

    if [[ -z "$entity_id" || "$entity_id" == "null" ]]; then
        echo "Failed to get entity ID for ${id}"
        continue
    fi

    echo "Found entity ${entity_id}"

    # Retrieve sequence
    sequence=$(curl -fsSL \
        "https://data.rcsb.org/rest/v1/core/polymer_entity/${pdb_id}/${entity_id}" |
        jq -r '.entity_poly.pdbx_seq_one_letter_code_can')

    if [[ -z "$sequence" || "$sequence" == "null" ]]; then
        echo "Failed to get sequence for ${id}"
        continue
    fi

    sequence=$(echo "$sequence" | tr -d '[:space:]')

    outfile="${OUTPUT_DIR}/${id}.fasta"

    {
        echo ">${id}"
        echo "$sequence"
    } > "$outfile"

    echo "Saved ${outfile}"

done < "$FILE"