from pathlib import Path
import os
import argparse
from Bio import PDB
from scipy.spatial import cKDTree
import json
import gzip
import itertools

STANDARD_AMINO_ACIDS = {
"ALA",
"ARG",
"ASN",
"ASP",
"CYS",
"GLN",
"GLU",
"GLY",
"HIS",
"ILE",
"LEU",
"LYS",
"MET",
"PHE",
"PRO",
"SER",
"THR",
"TRP",
"TYR",
"VAL",
}

def get_interfaces_for_all_assemblies(
        assembly_dir : Path,
        rcsb_id : str,
        distance_threshold: float,
        residue_count_threshold: int,
        outpath: str,
        rosetta: bool = False,
    ):
        """
        Identify all interfaces in all bio-assemblies and save results as xml file

        Parameters:
            distance_threshold (float): Distance cutoff in Ångstroms
            residue_count_threshold (int): Minimum number of residues per interface
            outpath (str): Directory where the final XML file and the Rosetta results are stored (if rosetta=True)
            rosetta (bool): Should the interfaces be analysed with Rosetta?

        Returns:
            dict: Interaction data dictionary
        """
        assembly_files = []
        for file in os.listdir(assembly_dir):
             assembly_files.append(str(assembly_dir / Path(file)))
        # Check that the structure files are present
        if not assembly_files:
            raise Exception("No structure files found")

        # Initialize dictionary
        interactions = {"pdb_code": rcsb_id.upper(), "method":"distance", "interface": []}
        interface_id = 1
        # Go through all assembly files and identify interfaces
        for infile in assembly_files:
            full_pdb_code = infile.split("/")[-1].split(".")[0]
            interactions, interface_id = _find_chain_interactions(
                complex_file=infile,
                full_pdb_code=full_pdb_code,
                interactions=interactions,
                interface_id=interface_id,
                outpath=outpath,
                rosetta=rosetta,
                distance_threshold=distance_threshold,
                residue_count_threshold=residue_count_threshold,
            )
        # Save results as json file
        with open(outpath / Path("interfaces.json"), "w") as f:
            json.dump(interactions, f, indent=2)

def _identify_contacts(
        complex_file: str,
        interface_id: int,
        full_pdb_code: str,
        outpath: str,
        chain_id1: str,
        chain_id2: str,
        atoms1: list,
        atoms2: list,
        tree1: cKDTree,
        tree2: cKDTree,
        distance_threshold: float,
        residue_count_threshold: int,
        rosetta: bool,
    ):
        """
        Identifies atom contacts and interface residues given a specific chain pair

        Parameters:
            complex_file (str): Path to complex file
            interface_id (int): Interface ID counter
            full_pdb_code (str): Full PDB code (e.g. 1234-assembly1)
            outpath (str): Directory where to store rosetta results if rosetta=True
            chain_id1 (str): ID of first chain
            chain_id2 (str): ID of second chain
            atoms1 (list): List with atoms of first chain
            atoms2 (list): List with atoms of second chain
            tree1 (cKDTree): Scipy cKDTree of first chain
            tree2 (cKDTree): Scipy cKDTree of second chain
            distance_threshold (float): Distance cutoff in Ångstroms
            residue_count_threshold (int): Minimum number of residues per interface
            rosetta (bool): Should Rosetta interface analysis be carried out?

        Returns:
            dict: Updated interaction data dictionary
            int: Updated interface ID counter
        """
        # Find atoms within threshold
        indices = tree1.query_ball_tree(tree2, distance_threshold)

        # Continue with next chain pair if no contacts were found
        if not any(sublist for sublist in indices):
            return None

        # Save interface data in dictionary
        molecules = [
            {
                "id": 1,
                "chain_id": chain_id1.split("-")[0],
                "chain_id_assembly": chain_id1,
                "class": "Protein",
                "int_nres": 0,
                "int_residues": {"residue": []},
            },
            {
                "id": 2,
                "chain_id": chain_id2.split("-")[0],
                "chain_id_assembly": chain_id2,
                "class": "Protein",
                "int_nres": 0,
                "int_residues": {"residue": []},
            },
        ]
        interface = {
            "id": interface_id,
            "assembly": full_pdb_code,
            "molecule": molecules,
        }

        # Get all interface residues
        interface, int_nres_1, int_nres_2 = _get_interface_residues(
            atom_indices=indices, atoms1=atoms1, atoms2=atoms2, interface=interface
        )
        interface["molecule"][0]["int_nres"] = int_nres_1
        interface["molecule"][1]["int_nres"] = int_nres_2

        # Keep interface if it contains at least residue_count_threshold interface residues
        if (int_nres_1 + int_nres_2) >= residue_count_threshold:
            return interface
        # Return None if interface too small
        return None

def _get_interface_residues(
        atom_indices: list, atoms1: list, atoms2: list, interface: dict
    ):
        """
        Get inteface residues for given interface

        Parameters:
            atom_indices (list): List of atom indices that are in contact
            atoms1 (list): List with atoms of first chain
            atoms2 (list): List with atoms of second chain
            interface (dict): Interface data dictionary

        Returns:
            tuple: Updated interface dictionary and number of interface residues
        """
        residues1 = []
        residues2 = []
        # Get interface residues
        for idx1, idx_list in enumerate(atom_indices):
            for idx2 in idx_list:
                _, _, _, res1 = atoms1[idx1]
                _, _, _, res2 = atoms2[idx2]

                residue1 = {"name": res1.get_resname(), "seq_num": res1.get_id()[1], "ins_code": res1.get_id()[2]}
                residue2 = {"name": res2.get_resname(), "seq_num": res2.get_id()[1], "ins_code":res2.get_id()[2]}

                if res1 not in residues1:
                    interface["molecule"][0]["int_residues"]["residue"].append(residue1)
                    residues1.append(res1)
                if res2 not in residues2:
                    interface["molecule"][1]["int_residues"]["residue"].append(residue2)
                    residues2.append(res2)

        # Sort residues
        interface["molecule"][0]["int_residues"]["residue"] = sorted(
            interface["molecule"][0]["int_residues"]["residue"],
            key=lambda x: x["seq_num"],
        )
        interface["molecule"][1]["int_residues"]["residue"] = sorted(
            interface["molecule"][1]["int_residues"]["residue"],
            key=lambda x: x["seq_num"],
        )

        # Count number of interface residues
        int_nres_1 = len(interface["molecule"][0]["int_residues"]["residue"])
        int_nres_2 = len(interface["molecule"][1]["int_residues"]["residue"])

        return interface, int_nres_1, int_nres_2

def _find_chain_interactions(
        complex_file: str,
        full_pdb_code: str,
        interactions: dict,
        interface_id: int,
        outpath: str,
        rosetta: bool = False,
        distance_threshold: float = 5.0,
        residue_count_threshold: int = 0,
    ):
        """
        Identifies interacting chains and corresponding residues in a protein bio-assembly structure.

        Parameters:
            complex_file (str): Path to the complex file
            full_pdb_code (str): Full PDB code (e.g. 1234-assembly1)
            interactions (dict): Interaction data
            interface_id (int): Interface ID counter
            outpath (str): Directory where to store Rosetta results if rosetta=True
            distance_threshold (float): Distance cutoff in Ångstroms
            residue_count_threshold (int): Minimum number of residues per interface

        Returns:
            dict: Updated interaction data dictionary
            int: Updated interface ID counter
        """

        # Determine the parser based on file extension
        if complex_file.endswith(".cif") or complex_file.endswith(".cif.gz"):
            parser = PDB.MMCIFParser(QUIET=True)
        elif complex_file.endswith(".pdb") or complex_file.endswith(".pdb.gz"):
            parser = PDB.PDBParser(QUIET=True)
        else:
            raise ValueError(
                "Unsupported file format. Please provide a PDB or mmCIF file."
            )

        # Load structure
        if complex_file.endswith(".gz"):
            with gzip.open(complex_file, "rt") as file:
                structure = parser.get_structure("protein", file)
        else:
            structure = parser.get_structure("protein", complex_file)

        # Check if structure contains more than one chain
        if len(list(structure.get_chains())) == 1:
            return interactions, interface_id

        # Extract all atoms, grouped by chain (exclude non-protein chains and hydrogen atoms)
        chain_atoms = {}
        chain_trees = {}
        for chain in structure.get_chains():
            residues = [
                res
                for res in chain
                if res.get_resname() in STANDARD_AMINO_ACIDS
            ]
            atoms = [
                atom
                for res in residues
                for atom in res.get_atoms()
                if not atom.get_name().startswith("H")
            ]
            if atoms:
                chain_atoms[chain.id] = [
                    (atom.get_coord(), atom, chain.id, atom.get_parent())
                    for atom in atoms
                ]
                chain_trees[chain.id] = cKDTree([atom.get_coord() for atom in atoms])  # Precompute KDTree

        interacting_chains = []

        # Identify contacts
        for (chain_id1, atoms1), (chain_id2, atoms2) in itertools.combinations(
            chain_atoms.items(), 2
        ):
            interface = _identify_contacts(
                complex_file=complex_file,
                interface_id=interface_id,
                full_pdb_code=full_pdb_code,
                outpath=outpath,
                chain_id1=chain_id1,
                chain_id2=chain_id2,
                atoms1=atoms1,
                atoms2=atoms2,
                tree1=chain_trees[chain_id1],
                tree2=chain_trees[chain_id2],
                distance_threshold=distance_threshold,
                residue_count_threshold=residue_count_threshold,
                rosetta=rosetta,
            )
            # If interface was detected between chains, add data to dict
            if interface:
                interactions["interface"].append(interface)
                interacting_chains.append((chain_id1, chain_id2))
                interface_id += 1

        return interactions, interface_id

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("assembly_dir", type=Path)
    parser.add_argument("-o", "--output_dir", type=Path)
    parser.add_argument("--distance_threshold", type=float, default=6.0)
    parser.add_argument("--residue_count_threshold", type=int, default=100)
    args = parser.parse_args()

    assembly_dir = args.assembly_dir
    rcsb_id = assembly_dir.stem
    distance_threshold = args.distance_threshold
    residue_count_threshold = args.residue_count_threshold
    outpath = args.output_dir / Path(rcsb_id)
    os.makedirs(outpath, exist_ok=True)

    get_interfaces_for_all_assemblies(assembly_dir=assembly_dir, rcsb_id=rcsb_id, distance_threshold=distance_threshold, residue_count_threshold=residue_count_threshold, outpath=outpath)