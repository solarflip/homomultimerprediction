
"""
Search PDB for homologous protein entries and extract RCSB bioassembly symmetry
information directly from the RCSB JSON API.
"""

import argparse
import re
import time
from typing import Dict, List, Tuple, Any
import requests
import pandas as pd
from Bio.Align import PairwiseAligner

# RCSB API endpoints
SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
ASSEMBLY_URL = "https://data.rcsb.org/rest/v1/core/assembly/{pdb_id}/{assembly_id}"
POLYMER_ENTITY_URL = (
    "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
)
# HTTP connections
SESSION = requests.Session()
SESSION.headers.update({
    "User-Agent": "PDB-Symmetry-Search/1.0"
})


def clean_protein_sequence(seq: str) -> str:
    """
    Clean protein sequence before submitting it to RCSB.

    Args:
        seq (str): protein sequence
    
    Returns:
        str: cleaned protein sequence
    """
    seq = seq.upper()
    seq = re.sub(r"[^ACDEFGHIKLMNPQRSTVWYXBZUO]", "", seq)
    return seq


def read_fasta(path: str) -> Dict[str, str]:
    """
    Read FASTA file into: protein_id -> sequence

    Args:
        path (str): path to FASTA file

    Returns:
        Dict[str, str]: mapping of protein_id to sequence
    """
    seqs = {}
    name = None
    chunks = []

    with open(path) as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    seqs[name] = clean_protein_sequence("".join(chunks))

                name = line[1:].split("|")[1] if "|" in line else line[1:]
                chunks = []
            else:
                chunks.append(line)

    if name is not None:
        seqs[name] = clean_protein_sequence("".join(chunks))

    return seqs


def get_json(url: str) -> dict:
    """
    GET JSON from RCSB API.

    Args:
        url (str): URL to fetch

    Returns:
        dict: parsed JSON response
    """
    response = SESSION.get(url, timeout=60)

    if response.status_code == 204:
        return {}

    response.raise_for_status()
    return response.json()


def rcsb_sequence_search(
    sequence: str,
    identity_cutoff: float,
    evalue_cutoff: float,
    max_hits: int,
) -> List[Tuple[str, str, float]]:
    """
    Search RCSB/PDB for homologous protein polymer entities.

    Args:
        sequence (str): protein sequence to search
        identity_cutoff (float): minimum sequence identity for hits
        evalue_cutoff (float): maximum e-value for hits
        max_hits (int): maximum number of hits to return

    Returns:
        List[Tuple[str, str, float]]: List of tuples: (pdb_id, entity_id, score)
    """
    query = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": evalue_cutoff,
                "identity_cutoff": identity_cutoff,
                "sequence_type": "protein",
                "target": "pdb_protein_sequence",
                "value": sequence,
            },
        },
        "request_options": {
            "paginate": {
                "start": 0,
                "rows": max_hits,
            },
            "scoring_strategy": "sequence",
        },
        "return_type": "polymer_entity",
    }

    response = SESSION.post(SEARCH_URL, json=query, timeout=60)

    if response.status_code == 204:
        return []

    if response.status_code != 200:
        raise RuntimeError(
            f"RCSB search failed with HTTP {response.status_code}: "
            f"{response.text[:500]}"
        )

    try:
        data = response.json()
    except Exception:
        raise RuntimeError(
            "RCSB did not return valid JSON. "
            f"HTTP status: {response.status_code}. "
            f"Response text: {response.text[:500]}"
        )

    hits = []

    # Extract PDB ID, entity ID, and score from the search results
    for item in data.get("result_set", []):
        identifier = item["identifier"]
        score = item.get("score", None)

        match = re.match(r"^([A-Za-z0-9]+)_(\d+)$", identifier)
        if not match:
            continue

        pdb_id = match.group(1).upper()
        entity_id = match.group(2)

        hits.append((pdb_id, entity_id, score))

    return hits


def get_pdb_entity_sequence_and_asym_ids(pdb_id: str, entity_id: str) -> Tuple[str, List[str]]:
    """
    Get canonical polymer entity sequence and asymmetric IDs from RCSB JSON API.

    Args:
        pdb_id (str): PDB ID (case-insensitive)
        entity_id (str): polymer entity ID (case-insensitive)

    Returns:
        Tuple[str, List[str]]: (sequence, asym_ids)
    """
    url = POLYMER_ENTITY_URL.format(
        pdb_id=pdb_id.lower(),
        entity_id=entity_id,
    )

    data = get_json(url)

    # Extract the canonical sequence
    seq = data.get("entity_poly", {}).get("pdbx_seq_one_letter_code_can", "")

    # Extract asymmetric unit IDs for the polymer entity
    identifiers = data.get(
        "rcsb_polymer_entity_container_identifiers",
        {},
    )
    asym_ids = identifiers.get("asym_ids", [])
    return clean_protein_sequence(seq), asym_ids


def calculate_identity_and_coverage(query_seq: str, hit_seq: str) -> Tuple[float, float, float]:
    """
    Calculate sequence identity and coverage between two sequences using global alignment.
    
    Args:
        query_seq (str): Query sequence
        hit_seq (str): Hit sequence

    Returns:
        Tuple[float, float, float]: (identity, query_coverage, target_coverage)
    """
    if not query_seq or not hit_seq:
        return 0.0, 0.0, 0.0

    # Use Biopython's PairwiseAligner for global alignment
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    alignment = aligner.align(query_seq, hit_seq)[0]

    covered_query_residues = 0
    covered_hit_residues = 0
    identical = 0

    for q_block, h_block in zip(alignment.aligned[0], alignment.aligned[1]):
        q_start, q_end = q_block
        h_start, h_end = h_block

        block_len = min(q_end - q_start, h_end - h_start)

        covered_query_residues += block_len
        covered_hit_residues += block_len

        for i in range(block_len):
            if query_seq[q_start + i] == hit_seq[h_start + i]:
                identical += 1

    if covered_query_residues == 0:
        return 0.0, 0.0, 0.0

    identity = identical / covered_query_residues
    query_coverage = covered_query_residues / len(query_seq)
    target_coverage = covered_hit_residues / len(hit_seq)

    return identity, query_coverage, target_coverage


def get_assembly_ids_for_entry(pdb_id: str) -> List[str]:
    """
    Get biological assembly IDs for a PDB entry.

    Args:
        pdb_id (str): PDB ID (case-insensitive)

    Returns:
        List[str]: List of assembly IDs
    """
    url = ENTRY_URL.format(pdb_id=pdb_id.lower())
    data = get_json(url)

    return (
        data.get("rcsb_entry_container_identifiers", {})
        .get("assembly_ids", [])
    )

def normalize_stoichiometry(stoichiometry) -> str:
    """
    Convert RCSB stoichiometry field to a clean string.

    Args:
        stoichiometry (str or list): RCSB stoichiometry field

    Returns:
        str: normalized stoichiometry string
    """

    if stoichiometry is None:
        return ""

    if isinstance(stoichiometry, list):
        if len(stoichiometry) == 0:
            return ""

        # Common case: ["A2"]
        if len(stoichiometry) == 1:
            return str(stoichiometry[0]).strip()

        # Fallback: join multiple entries
        return "".join(str(x).strip() for x in stoichiometry)

    return str(stoichiometry).strip()

def is_homomeric_stoichiometry(stoichiometry) -> bool:
    """
    Check whether stoichiometry looks homomeric.
    Kept:
        A, A2, A3, A12, A24, A60
    Excluded:
        AB, A2B2, ABC, A3B3

    Args:
        stoichiometry (str): RCSB stoichiometry field

    Returns:
        bool: True if homomeric, False otherwise
    """

    stoichiometry = normalize_stoichiometry(stoichiometry)

    if not stoichiometry or stoichiometry in [".", "?"]:
        return False

    return re.fullmatch(r"A\d*", stoichiometry) is not None


def extract_rcsb_bioassembly_symmetry_from_api(
    pdb_id: str,
    entity_id: str,
    asym_ids: List[str],
    verbose: bool = False,
) -> List[Dict[str, Any]]:
    """
    Extract bioassembly symmetry directly from RCSB JSON API.
    Keeps only homomeric assemblies based on Global Stoichiometry.

    Args:
        pdb_id (str): PDB ID (case-insensitive)
        entity_id (str): polymer entity ID (case-insensitive)
        asym_ids (List[str]): asymmetric unit IDs for the polymer entity
        verbose (bool): whether to print detailed information

    Returns:
        List[Dict[str, Any]]: List of dictionaries with symmetry information for each homomeric assembly.
    """
    output = []

    # Get all assembly IDs for the PDB entry
    assembly_ids = get_assembly_ids_for_entry(pdb_id)

    # Go through each assembly and extract symmetry information
    for assembly_id in assembly_ids:
        url = ASSEMBLY_URL.format(
            pdb_id=pdb_id.lower(),
            assembly_id=assembly_id,
        )

        assembly_json = get_json(url)

        if not assembly_json:
            continue

        # Extract asymmetric unit IDs for the assembly
        assembly_infos = assembly_json.get("pdbx_struct_assembly_gen", [])
        asym_ids_assembly = []
        for assembly_info in assembly_infos:
            current_assembly_id = assembly_info.get("assembly_id")
            if not current_assembly_id == assembly_id:
                continue
            asym_ids_assembly.extend(assembly_info.get("asym_id_list", []))
        
        # Safety check: Check if there is overlap between the asym_ids of the polymer entity and the asym_ids of the assembly
        if not set(asym_ids).intersection(set(asym_ids_assembly)):
            continue

        # Extract symmetry information from the assembly JSON
        symmetry_rows = assembly_json.get("rcsb_struct_symmetry", [])

        if isinstance(symmetry_rows, dict):
            symmetry_rows = [symmetry_rows]

        for sym in symmetry_rows:
            kind = sym.get("kind", "")
            symbol = sym.get("symbol", "")

            # Keep global symmetry annotation
            if kind and kind != "Global Symmetry":
                continue

            stoichiometry = normalize_stoichiometry(
                sym.get("stoichiometry")
                or sym.get("global_stoichiometry")
                or assembly_json.get("rcsb_assembly_info", {}).get("stoichiometry")
                or ""
            )
            global_symmetry = (
                symbol
                or sym.get("global_symmetry")
                or ""
            )

            # Check if the assembly is homomeric based on stoichiometry
            is_homomeric = is_homomeric_stoichiometry(stoichiometry)

            if verbose:
                print(
                    f"{pdb_id} entity {entity_id} assembly {assembly_id}: "
                    f"stoichiometry={stoichiometry}, "
                    f"symmetry={global_symmetry}, "
                    f"kind={kind}, "
                    f"homomeric={is_homomeric}"
                )

            # Only keep homomeric assemblies
            if not is_homomeric:
                continue

            output.append(
                {
                    "pdb_id": pdb_id,
                    "entity_id": entity_id,
                    "assembly_id": assembly_id,
                    "global_stoichiometry": stoichiometry,
                    "global_symmetry": global_symmetry,
                    "symmetry_kind": kind,
                    "is_homomeric": True,
                }
            )

    return output


def main():
    parser = argparse.ArgumentParser(
        description="Extract symmetry information from PDB bioassemblies."
    )

    parser.add_argument(
        "--fasta",
        required=True,
        help="Input FASTA file with query protein sequences.",
    )

    parser.add_argument(
        "--out",
        default="pdb_symmetry_comparison.tsv",
        help="Output TSV file.",
    )

    parser.add_argument(
        "--identity",
        type=float,
        default=0.30,
        help="Initial RCSB sequence-search identity cutoff.",
    )

    parser.add_argument(
        "--final-identity",
        type=float,
        default=0.35,
        help="Minimum final identity after global pairwise alignment.",
    )

    parser.add_argument(
        "--query-coverage",
        type=float,
        default=0.50,
        help="Minimum query coverage after global pairwise alignment.",
    )

    parser.add_argument(
        "--target-coverage",
        type=float,
        default=0.90,
        help="Minimum target/PDB entity coverage after global alignment.",
    )

    parser.add_argument(
        "--evalue",
        type=float,
        default=1e-3,
        help="Initial RCSB sequence-search e-value cutoff.",
    )

    parser.add_argument(
        "--max-hits",
        type=int,
        default=50,
        help="Maximum PDB homologs to retrieve per query protein.",
    )

    parser.add_argument(
        "--sleep",
        type=float,
        default=0.5,
        help="Delay between requests.",
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print detailed symmetry information for each bioassembly.",
    )

    args = parser.parse_args()

    sequences = read_fasta(args.fasta)

    all_rows = []
    no_pdb_homomer = []
    no_pdb_hit = []
    failed_proteins = []


    for protein_id, query_seq in sequences.items():

        print(f"Searching PDB homologs for {protein_id}")
        protein_has_homomer = False
        protein_has_hit = False

        try:
            hits = rcsb_sequence_search(
                sequence=query_seq,
                identity_cutoff=args.identity,
                evalue_cutoff=args.evalue,
                max_hits=args.max_hits,
            )
        except Exception as exc:
            print(f"WARNING: RCSB search failed for {protein_id}: {exc}")
            failed_proteins.append(protein_id)
            continue

        if not hits:
            no_pdb_hit.append(protein_id)
            continue

        for pdb_id, entity_id, search_score in hits:
            hit_id = f"{pdb_id}_{entity_id}"

            try:
                pdb_entity_seq, asym_ids = get_pdb_entity_sequence_and_asym_ids(pdb_id, entity_id)

                final_identity, query_coverage, target_coverage = calculate_identity_and_coverage(
                    query_seq=query_seq,
                    hit_seq=pdb_entity_seq,
                )

                if final_identity < args.final_identity:
                    continue

                if query_coverage < args.query_coverage:
                    continue

                if target_coverage < args.target_coverage:
                    continue

                protein_has_hit = True

                symmetry_rows = extract_rcsb_bioassembly_symmetry_from_api(
                    pdb_id=pdb_id,
                    entity_id=entity_id,
                    asym_ids=asym_ids,
                    verbose=args.verbose,
                )

                for sym_row in symmetry_rows:
                    protein_has_homomer = True
                    all_rows.append(
                        {
                            "protein_id": protein_id,
                            "pdb_hit": hit_id,
                            "final_identity": round(final_identity, 4),
                            "query_coverage": round(query_coverage, 4),
                            "target_coverage": round(target_coverage, 4),
                            **sym_row,
                        }
                    )

                time.sleep(args.sleep)

            except Exception as exc:
                print(f"WARNING: failed for {protein_id} / {hit_id}: {exc}")
                continue
        # Remember proteins with no PDB hits
        if not protein_has_hit:
            no_pdb_hit.append(protein_id)
        # Remember proteins with no homomeric PDB hits
        if not protein_has_homomer:
            no_pdb_homomer.append(protein_id)

    df = pd.DataFrame(all_rows)
    df.to_csv(args.out, sep="\t", index=False)

    log_file = args.out.replace(".tsv", "_missing.txt")
    with open(log_file, "w") as f:
        for protein_id in no_pdb_homomer:
            f.write(protein_id + "," + "No PDB homomer found" + "\n")
        for protein_id in no_pdb_hit:
            f.write(protein_id + "," + "No PDB hit found" + "\n")
        for protein_id in failed_proteins:
            f.write(protein_id + "," + "Failed to process" + "\n")

    print(f"Saved output to: {args.out}")
    print(f"Number of retained homomeric PDB bioassembly matches: {len(df)}")


if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.2f} seconds")