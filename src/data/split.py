"""
Receives a list of clusters and splits them according to a predefined distribution (e.g. 80/10/10 or 70/20/10).
"""
import os
import argparse
from pathlib import Path
from datetime import datetime
from sklearn.model_selection import train_test_split
from src.utils.constants import (
    TRAIN,
    VAL,
    TEST, 
    RANDOM_STATE
)

def split(input : Path, name : str, output_dir : Path):
    """
    """
    print(f"Performing split: {TRAIN*100}/{VAL*100}/{TEST*100}")

    # read all lines
    if not os.path.exists(input):
        raise ValueError(f"Invalid input path: {input}")

    with open(input, "r") as f:
        content = f.read()

    cluster_ids = content.split(",")
    
    cluster_ids = [x.strip() for x in cluster_ids]
    cluster_ids_unique = list(set(cluster_ids))
    print(f"Read {len(cluster_ids_unique)} cluster")


    # perform randomized split
    train_cl, test_cl = train_test_split(cluster_ids_unique, test_size=(VAL+TEST), random_state=RANDOM_STATE, shuffle=True)
    test_cl, val_cl = train_test_split(test_cl, test_size=0.33, random_state=RANDOM_STATE, shuffle=True)

    print(f"Clusters in train: {len(train_cl)}\nClusters in val: {len(val_cl)}\nClusters in test: {len(test_cl)}")

    # save
    output_dir = output_dir / Path(name)
    os.makedirs(output_dir, exist_ok=True)

    train_cl_out = "\n".join(train_cl)
    train_path = output_dir / Path("train.txt")
    with open(train_path, "a+") as f:
        f.write(train_cl_out)

    val_cl_out = "\n".join(val_cl)
    val_path = output_dir / Path("val.txt")
    with open(val_path, "a+") as f:
        f.write(val_cl_out)

    test_cl_out = "\n".join(test_cl)
    test_path = output_dir / Path("test.txt")
    with open(test_path, "a+") as f:
        f.write(test_cl_out)

    info_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    info_str += f"Performing split: {TRAIN*100}/{VAL*100}/{TEST*100}\n" 
    info_str += f"Read {len(cluster_ids_unique)} cluster\n"
    info_str += f"Clusters in train: {len(train_cl)}\nClusters in val: {len(val_cl)}\nClusters in test: {len(test_cl)}"
    info_path = output_dir / Path("INFO")
    with open(info_path, "a+") as f:
        f.write(info_str)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path, help="Path to a file containing a cluster id in each row.")
    parser.add_argument("name", type=str, help="Name of the folder, which will be created in order to put the final list inside of it.")
    parser.add_argument("-o", type=Path, help="Path to the output folder.", default=Path("."))
    args = parser.parse_args()

    split(args.input, args.name, args.o)
