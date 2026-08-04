"""
Custom dataset class for training, validation and testing.

Allows to select which types of input are to selected.
Class needs access to label file as well as the raw data dir.
"""
import torch
from torch.utils.data import Dataset
from pathlib import Path
import pandas as pd
import os
from src.utils.constants import (
    LABEL_COLUMN
)

class MultimodalDataset(Dataset):
    def __init__(self, esm2_dir : Path, metadata_file : Path, cluster_ids : Path):
        """
        PARAMS:
            - esm2_dir Path, Path to the directory where the esm2 embeddings are
            - metadata_file, Path to the file that contains a mapping from structure to label
            - cluster_ids, Path to the cluster_ids that are relevant for this dataset

        Introduces self.data:
        - must contain:
            - label
            - some feature
            - sequence?
            - pdb file?
        """
        # read dataset relevant cluster ids 
        if not os.path.exists(cluster_ids):
            raise ValueError(f"File containing cluster ids does not exist: {cluster_ids}")
        with open(cluster_ids, "r") as f:
            relevant_ids = f.readlines()
            relevant_ids = [int(x.strip()) for x in relevant_ids]

        # read metadata_file
        metadata = pd.read_csv(metadata_file)

        # create subset with relevant ones
        metadata_subset = metadata[metadata["CLUSTER"].isin(relevant_ids)]

        # create self.data
        self.data = []
        self.labels = []

        for index, row in metadata_subset.iterrows():
            self.labels.append(row["SYMM"])
            self.data.append(row["CHAINID"])

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        sample = self.data[idx]
        label = self.labels[idx]
        return sample, label