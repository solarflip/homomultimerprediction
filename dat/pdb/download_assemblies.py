import os
from urllib import request, error
from pathlib import Path
import argparse

def fetch_bioassemblies(rcsb_id : str, outdir: Path):
    """
    Script from Laura P.

    Download all bioassembly files for given PDB entry and saves the paths as attribute

    Parameters:
        outdir (str): Directory where the files are stored

    Returns:
        list: List of file paths
    """
    assembly_id = 1
    bioassemblies = []
    id_specific_output = outdir / Path(rcsb_id)
    os.makedirs(id_specific_output, exist_ok=True)

    while True:
        url = f"https://files.rcsb.org/download/{rcsb_id}-assembly{assembly_id}.cif.gz"
        file_path = f"{id_specific_output}/{rcsb_id}-assembly{assembly_id}.cif.gz"
        try:
            # Download the file
            request.urlretrieve(url, file_path)
            bioassemblies.append(file_path)
            assembly_id += 1

        except error.HTTPError as e:
            if e.code == 404:
                break # Stop if 404 (Not Found) error occurs
            else:
                print(f"HTTP error: {e}")
                break # Stop on other HTTP errors

        except error.URLError as e:
            print(f"URL error: {e}")
            break # Stop on connection failure

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("rcsb_id", type=str)
    parser.add_argument("-o", "--output_dir", type=Path, default=".")
    args = parser.parse_args()

    fetch_bioassemblies(args.rcsb_id, args.output_dir)