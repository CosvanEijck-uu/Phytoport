"""
Phosphorylation Site Comparison Tool

This script compares predicted phosphorylation sites between two protein
sequences based on a TSV file of prediction results. It identifies
common and unique sites, applies score cutoffs, and writes a summary
report to an output TSV file.

Usage:
    python -m compare.ps_sites <tsv_file> <seq1> <seq2> [--cut-off 0.5] [-o output.tsv]
"""

import argparse
import csv
import os
import sys
from collections import defaultdict


def extract_sites(tsv_file: str, headers: list[str]) -> dict[str, list[tuple[int, str, float]]]:
    """
    Extract phosphorylation sites for specified sequence headers from a TSV file.

    Parameters
    ----------
    tsv_file : str
        Path to the TSV file containing phosphorylation predictions.
    headers : list[str]
        List of sequence headers to extract data for (last field after '|').

    Returns
    -------
    dict[str, list[tuple[int, str, float]]]
        Dictionary mapping sequence headers to lists of (position, residue, score).

    Raises
    ------
    FileNotFoundError
        If the TSV file cannot be found.
    ValueError
        If no phosphorylation sites are found for the given headers.
    """
    if not os.path.isfile(tsv_file):
        raise FileNotFoundError(f"TSV file not found: {tsv_file}")

    headers_set = set(headers)
    results = defaultdict(list)

    # Determine allowed residues based on filename
    filename = os.path.basename(tsv_file)
    if "_Y" in filename:
        allowed_residues = {"Y"}
    elif "_SorT" in filename:
        allowed_residues = {"S", "T"}
    else:
        allowed_residues = {"S", "T", "Y"}

    try:
        with open(tsv_file, "r", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                if len(row) < 4:
                    continue

                full_header = row[0]
                short_header = full_header.split("|")[-1]

                try:
                    position = int(row[1])
                    residue = row[2].upper()
                    score = float(row[3])
                except ValueError:
                    continue

                if residue not in allowed_residues:
                    continue
                if short_header in headers_set:
                    results[short_header].append((position, residue, score))
    except OSError as e:
        raise RuntimeError(f"Failed to read TSV file: {e}")

    if not results:
        raise ValueError(f"No phosphorylation sites found for: {headers}")

    return results


def filter_sites_by_cutoff(sites: list[tuple[int, str, float]], cutoff: float) -> list[tuple[int, str, float]]:
    """
    Filter phosphorylation sites by score cutoff.

    Parameters
    ----------
    sites : list[tuple[int, str, float]]
        List of phosphorylation sites (position, residue, score).
    cutoff : float
        Minimum score required to include a site.

    Returns
    -------
    list[tuple[int, str, float]]
        Filtered list of phosphorylation sites.
    """
    return [s for s in sites if s[2] >= cutoff]


def find_common_sites(
    sites1: list[tuple[int, str, float]],
    sites2: list[tuple[int, str, float]]
) -> tuple[list[tuple[int, str, float, float]], list[tuple[int, str, float]], list[tuple[int, str, float]]]:
    """
    Identify common and unique phosphorylation sites between two sequences.

    Parameters
    ----------
    sites1, sites2 : list[tuple[int, str, float]]
        Phosphorylation sites for the two sequences.

    Returns
    -------
    tuple
        (common_sites, unique_to_first, unique_to_second)
    """
    common, unique1 = [], []
    unique2 = sites2.copy()

    for pos1, res1, score1 in sites1:
        match_found = False
        for s2 in sites2:
            pos2, res2, score2 = s2
            if pos1 == pos2 and res1 == res2:
                common.append((pos1, res1, score1, score2))
                if s2 in unique2:
                    unique2.remove(s2)
                match_found = True
                break
        if not match_found:
            unique1.append((pos1, res1, score1))

    # Sort results for readability
    for lst in (common, unique1, unique2):
        lst.sort(key=lambda x: x[0])

    return common, unique1, unique2


def save_to_tsv(
    output_file: str,
    header1: str,
    header2: str,
    common: list[tuple[int, str, float, float]],
    unique1: list[tuple[int, str, float]],
    unique2: list[tuple[int, str, float]]
) -> None:
    """
    Save phosphorylation site comparison results to a TSV file.

    Parameters
    ----------
    output_file : str
        Path to output TSV file.
    header1, header2 : str
        Sequence headers used in the comparison.
    common, unique1, unique2 : list
        Common and unique site lists.
    """
    try:
        with open(output_file, "w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow([
                "Category", "Position", "Residue",
                f"{header1}_Score", f"{header2}_Score"
            ])

            for pos, res, score1, score2 in common:
                writer.writerow(["Common", pos, res, score1, score2])

            for pos, res, score in unique1:
                writer.writerow([f"Unique_to_{header1}", pos, res, score, ""])

            for pos, res, score in unique2:
                writer.writerow([f"Unique_to_{header2}", pos, res, "", score])
    except OSError as e:
        raise RuntimeError(f"Failed to write output file: {e}")


def compare_phospho_sites(
    tsv_file: str,
    headers: list[str],
    cutoff: float = 0.0,
    output_file: str = "phospho_comparison.tsv"
) -> None:
    """
    Compare phosphorylation sites between two sequences and save results.

    Parameters
    ----------
    tsv_file : str
        Input TSV file of phosphorylation predictions.
    headers : list[str]
        Two sequence headers to compare.
    cutoff : float, optional
        Minimum prediction score to include (default: 0.0).
    output_file : str, optional
        Output TSV file (default: 'phospho_comparison.tsv').
    """
    if len(headers) != 2:
        raise ValueError(
            "Please provide exactly two sequence headers for comparison.")

    site_dict = extract_sites(tsv_file, headers)
    header1, header2 = headers

    sites1 = filter_sites_by_cutoff(site_dict.get(header1, []), cutoff)
    sites2 = filter_sites_by_cutoff(site_dict.get(header2, []), cutoff)

    common, unique1, unique2 = find_common_sites(sites1, sites2)
    save_to_tsv(output_file, header1, header2, common, unique1, unique2)
    print(f"Phosphorylation comparison saved to: {output_file}")


def main() -> None:
    """Command-line entry point."""
    parser = argparse.ArgumentParser(
        description="Compare predicted phosphorylation sites between two protein sequences."
    )
    parser.add_argument(
        "tsv_file", help="Input TSV file of phosphorylation predictions.")
    parser.add_argument(
        "sequence_headers",
        nargs=2,
        help="Two sequence headers to compare (only the part after the last '|' is needed)."
    )
    parser.add_argument(
        "--cut-off",
        type=float,
        default=0.0,
        help="Minimum prediction score to include a site (default: 0.0)."
    )
    parser.add_argument(
        "-o", "--output",
        default="phospho_comparison.tsv",
        help="Output TSV file (default: phospho_comparison.tsv)."
    )
    args = parser.parse_args()

    try:
        compare_phospho_sites(
            args.tsv_file,
            args.sequence_headers,
            cutoff=args.cut_off,
            output_file=args.output
        )
    except (FileNotFoundError, ValueError, RuntimeError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
