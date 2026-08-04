# Dataset
## Base and Structure Dataset
Original dataset from seq2symm (http://files.ipd.uw.edu/pub/seq2symm/datasets.zip)
```sh
curl http://files.ipd.uw.edu/pub/seq2symm/datasets.zip -o dataset.zip
unzip dataset.zip
```
Retrieve and download the `.zip` file. Intially use the `homomer_pdbids_hash_clusterid_labels_fullset.csv` file.
RCSB ids are extracted using:
```sh
./pdb/prep_batch_download.sh
```
and downloaded using 
```sh
./pdb/batch_download.sh -f homomer_pdbids_fullset.txt -o pdb/fullset -p
```

- [ ] unzip
- [ ] select single chain

## Sequence Dataset

- [ ] extract from pdb?

## ESM-2 Dataset

- [ ] embed using `esm-extract`

## Interface Residue Dataset

