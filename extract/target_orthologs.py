"""
Extract target orthologs from an OrthoFinder TSV.

This script reads an OrthoFinder orthogroups TSV file, identifies orthogroups
that contain target gene symbols (supporting wildcard '*' at the end), and
writes the matching orthogroups to a new TSV file.

Usage:
    python extract.target_orthologs <input_file> -t <targets> -o <output_file>
"""

import sys
import argparse
from pathlib import Path
from typing import List, Tuple


def extract_symbol(gene_id: str) -> str:
    """
    Extract the gene symbol from a UniProt-style ID.

    Example:
        sp|Q94BM7|SPA4_ARATH -> SPA4

    Parameters
    ----------
    gene_id : str
        UniProt-style gene identifier.

    Returns
    -------
    str
        Extracted gene symbol; empty string if parsing fails.
    """
    try:
        last_part = gene_id.split("|")[2]  # e.g., SPA4_ARATH
        return last_part.split("_")[0]  # -> SPA4
    except IndexError:
        return ""


def parse_line(line: str) -> Tuple[str, List[List[str]]]:
    """
    Split a line from the orthogroups TSV into orthogroup and per-organism genes.

    Parameters
    ----------
    line : str
        A line from the TSV file.

    Returns
    -------
    Tuple[str, List[List[str]]]
        Orthogroup ID and list of columns with genes.
    """
    parts = line.strip().split("\t")
    orthogroup = parts[0]
    columns = [[g.strip() for g in col.split(",") if g.strip()]
               for col in parts[1:]]
    return orthogroup, columns


def match_symbol(symbol: str, patterns: List[str]) -> bool:
    """
    Check if a gene symbol matches any of the target patterns.

    Supports wildcard '*' at the end for prefix matching.

    Parameters
    ----------
    symbol : str
        Gene symbol to check.
    patterns : List[str]
        List of target patterns (e.g., ['HY5', 'SPA*']).

    Returns
    -------
    bool
        True if symbol matches any pattern; False otherwise.
    """
    symbol = symbol.upper()
    for pat in patterns:
        pat = pat.upper()
        if pat.endswith("*") and symbol.startswith(pat[:-1]):
            return True
        elif symbol == pat:
            return True
    return False


def get_organism_names(tsv_file: Path) -> List[str]:
    """
    Extract organism names from the header row of an OrthoFinder TSV.

    Parameters
    ----------
    tsv_file : Path
        Path to the TSV file.

    Returns
    -------
    List[str]
        List of organism names (columns excluding first 'Orthogroup').
    """
    try:
        with tsv_file.open("r", encoding="utf-8") as f:
            header = f.readline().strip().split("\t")
            return header[1:]
    except Exception as e:
        raise RuntimeError(f"Error reading header from {tsv_file}: {e}")


def main(args: argparse.Namespace) -> None:
    """
    Extract target orthologs and save to output TSV.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    """
    if not args.target:
        raise ValueError(
            "At least one target must be specified using -t/--target.")

    targets = args.target
    input_file = Path(args.file)
    output_path = Path(args.output)

    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        raise RuntimeError(f"Failed to create output directory ' \
        {output_path.parent}': {e}")

    try:
        organism_names = get_organism_names(input_file)
        if not organism_names:
            raise ValueError(f"No organism names found in \
                             {input_file} header.")
    except Exception as e:
        raise RuntimeError(f"Failed to read organism names: {e}")

    try:
        with input_file.open("r", encoding="utf-8") as in_file, \
                output_path.open("w", encoding="utf-8") as out_file:

            # Write header
            out_file.write("Orthogroup\t" + "\t".join(organism_names) + "\n")
            next(in_file)  # skip header line

            for line_number, line in enumerate(in_file, start=2):
                try:
                    orthogroup, columns = parse_line(line)
                except Exception as e:
                    print(f"Skipping line {line_number} due to parse error: \
                    {e}", file=sys.stderr)
                    continue

                # Check if any gene in any column matches a target symbol
                if any(
                    match_symbol(extract_symbol(gene), targets)
                    for col in columns
                    for gene in col
                ):
                    all_cols = [";".join(col) for col in columns]
                    out_file.write(orthogroup + "\t" +
                                   "\t".join(all_cols) + "\n")

        print(f"Results successfully saved to {output_path}")

    except FileNotFoundError:
        raise FileNotFoundError(f"Input file '{input_file}' not found.")
    except Exception as e:
        raise RuntimeError(f"An error occurred while processing the file: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract target orthologs from an OrthoFinder TSV."
    )
    parser.add_argument("file", help="Path to input TSV file")
    parser.add_argument(
        "-o", "--output", required=True,
        help="Full path to output TSV file (will be created if it doesn't exist)"
    )
    parser.add_argument(
        "-t", "--target", nargs="+",
        help="One or more target gene symbols (wildcard '*' at end supported, e.g., SPA*)"
    )
    args = parser.parse_args()

    try:
        main(args)
    except (FileNotFoundError, ValueError, RuntimeError) as e:
        print(f"Error: {e}", file=sys.stderr)
        exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        exit(1)
