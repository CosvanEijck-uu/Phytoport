"""
fetch/proteins.py

Fetch FASTA sequences from UniProt for two given UniProt IDs,
rewrite headers, and save into a single FASTA file.

Example usage:
    python -m fetch.proteins P69905 P68871 -o sequences.fasta
"""

import argparse
import requests
import sys
import os
from itertools import product
from typing import List, Optional, Tuple


# ----------------------
# Structure protein logic
# ----------------------
def get_structure_protein_combinations_cross_species(
    targets_tsv: str, structure_proteins: List[str]
) -> List[Tuple[str, str]]:
    """
    Build cross-protein combinations (e.g., COP1×SPA1) within each species
    from the targets TSV file.

    Parameters:
        targets_tsv (str): Path to the TSV file containing target proteins.
        structure_proteins (List[str]): List of structure proteins to consider.

    Returns:
        List[Tuple[str, str]]: List of tuples representing protein pairs.

    Notes:
        - Only pairs detected within the same species are returned.
        - Returns an empty list if no pairs are found or the TSV is missing.
    """
    if not os.path.exists(targets_tsv):
        print(f"Warning: {targets_tsv} not found yet \
            — returning no structure pairs.")
        return []

    structure_proteins_lower = [sp.lower() for sp in structure_proteins]
    arabidopsis_map = {sp: set() for sp in structure_proteins_lower}
    solanum_map = {sp: set() for sp in structure_proteins_lower}

    with open(targets_tsv) as f:
        header = next(f)
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            _, arab_col, sol_col = parts

            detected_sps = set()
            for entry in arab_col.split(";"):
                entry_lower = entry.lower()
                for sp in structure_proteins_lower:
                    if sp in entry_lower:
                        uniprot_id = entry.split(
                            "|")[1] if "|" in entry else entry
                        arabidopsis_map[sp].add(uniprot_id)
                        detected_sps.add(sp)

            if detected_sps:
                for entry in sol_col.split(";"):
                    uniprot_id = entry.split("|")[1] if "|" in entry else entry
                    for sp in detected_sps:
                        solanum_map[sp].add(uniprot_id)

    arab_pairs = []
    sol_pairs = []
    for i, sp1 in enumerate(structure_proteins_lower):
        for sp2 in structure_proteins_lower[i + 1:]:
            arab_pairs.extend(
                product(arabidopsis_map[sp1], arabidopsis_map[sp2]))
            sol_pairs.extend(product(solanum_map[sp1], solanum_map[sp2]))

    all_pairs = sorted(arab_pairs) + sorted(sol_pairs)

    if not all_pairs:
        print(f"No structure protein pairs found in \
              {targets_tsv} for {structure_proteins}")
    else:
        print(f"Found structure protein pairs: {all_pairs}")

    return all_pairs


def rewrite_fasta_headers(fasta_file: str) -> None:
    """
    Rewrite FASTA headers to the format: >{last 5 letters of UniProt ID}|protein.

    Parameters:
        fasta_file (str): Path to the FASTA file to rewrite in-place.

    Raises:
        OSError: If the file cannot be read or written.
    """
    new_lines = []
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                uniprot_id = line.strip().split("|")[-1]
                new_header = f">{uniprot_id[-5:]}|protein"
                new_lines.append(new_header)
            else:
                new_lines.append(line.strip())
    with open(fasta_file, "w") as f:
        f.write("\n".join(new_lines) + "\n")


# ----------------------
# UniProt fetching logic
# ----------------------
def fetch_fasta(uniprot_id: str) -> Optional[str]:
    """
    Fetch a FASTA sequence from UniProt given a UniProt accession ID.

    Parameters:
        uniprot_id (str): UniProt accession ID (e.g., P69905).

    Returns:
        Optional[str]: FASTA text if successfully fetched; None on error.

    Notes:
        - Prints warnings to stderr on HTTP errors or request exceptions.
        - Times out after 10 seconds.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        if response.text.startswith(">"):
            return response.text.strip()
        print(f"Warning: Could not fetch sequence for \
            {uniprot_id} (HTTP {response.status_code})", file=sys.stderr)
        return None
    except requests.exceptions.Timeout:
        print(f"❌ Error: Request for {uniprot_id} timed out.", file=sys.stderr)
    except requests.exceptions.RequestException as e:
        print(f"❌ Error fetching {uniprot_id}: {e}", file=sys.stderr)
    return None


# ----------------------
# Main CLI logic
# ----------------------
def main(uniprot_ids: List[str], output_file: str) -> None:
    """
    Fetch FASTA sequences for two UniProt IDs and save them to one FASTA file.

    Parameters:
        uniprot_ids (List[str]): List of two UniProt IDs.
        output_file (str): Output FASTA file path.

    Raises:
        SystemExit: Exits if no sequences were fetched or on write error.
    """
    fasta_contents = []

    for uniprot_id in uniprot_ids:
        print(f"Fetching {uniprot_id} from UniProt...")
        fasta = fetch_fasta(uniprot_id)
        if fasta:
            fasta_contents.append(fasta)
        else:
            print(f"Skipping {uniprot_id} due to fetch error.",
                  file=sys.stderr)

    if not fasta_contents:
        print("❌ No sequences were fetched. Exiting.", file=sys.stderr)
        sys.exit(1)

    try:
        os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
        with open(output_file, "w") as f:
            f.write("\n\n".join(fasta_contents) + "\n")
        print(f"Successfully saved \
            {len(fasta_contents)} sequences to {output_file}")
    except OSError as e:
        print(f"❌ Error writing to {output_file}: {e}", file=sys.stderr)
        sys.exit(1)

    # Rewrite headers in-place
    rewrite_fasta_headers(output_file)
    print(f"FASTA headers rewritten in {output_file}")


# ----------------------
# CLI entry point
# ----------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fetch FASTA sequences for two UniProt IDs and save to one file."
    )
    parser.add_argument(
        "uniprot_ids",
        nargs=2,
        help="Two UniProt IDs to fetch (e.g., P69905 P68871)"
    )
    parser.add_argument(
        "-o", "--output",
        default="sequences.fasta",
        help="Output FASTA file (default: sequences.fasta)"
    )
    args = parser.parse_args()
    main(args.uniprot_ids, args.output)
