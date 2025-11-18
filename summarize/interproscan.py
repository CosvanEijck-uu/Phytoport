"""
Summarize InterProScan domain comparison results across orthogroups.

Usage:
    python -m summarize.interpro <input_dirs_comma> <output_file>
"""

import os
import csv
import sys
from collections import defaultdict
from typing import List, Dict, Set


def load_feature_sets(dirs: List[str], mode: str = "domains") -> Dict[str, Set[str]]:
    """
    Load domain or feature sets from multiple directories containing TSV files.

    Each TSV file is expected to represent a pairwise comparison between two proteins,
    named in the format: `<protein1>__vs__<protein2>.tsv`.
    The function parses these files and constructs a mapping of each protein
    to its associated domains or features (from the 5th column of each row).

    Args:
        dirs (List[str]): List of directory paths containing TSV comparison files.
        mode (str): Feature type to load (default: "domains"). Currently used only for labeling purposes.

    Returns:
        Dict[str, Set[str]]: A dictionary mapping each protein to a set of domains/features.

    Notes:
        - Proteins from both columns of the TSV file are included.
        - Files not ending with `.tsv` or not following the naming convention are skipped.
    """
    data: Dict[str, Set[str]] = defaultdict(set)

    for dir_path in dirs:
        for fname in sorted(os.listdir(dir_path)):
            if not fname.endswith(".tsv"):
                continue
            parts = fname.replace(".tsv", "").split("__vs__")
            if len(parts) != 2:
                continue
            prot1, prot2 = parts

            with open(os.path.join(dir_path, fname), encoding="utf-8") as f:
                reader = csv.reader(f, delimiter="\t")
                for row in reader:
                    if len(row) > 4:
                        dom = row[4]
                        data[prot1].add(dom)
                        data[prot2].add(dom)
    return data


def write_overlap_matrix(feature_map: Dict[str, Set[str]], outfile: str, header_label: str = "Domains") -> None:
    """
    Write an overlap matrix representing shared features between all protein pairs.

    For each pair of proteins, the number of overlapping domains/features
    is calculated and written to a tab-separated output file.

    Args:
        feature_map (Dict[str, Set[str]]): Mapping of protein IDs to sets of features.
        outfile (str): Path to the output TSV file.
        header_label (str): Label for the first column header (default: "Domains").

    Returns:
        None

    Notes:
        - If the feature_map is empty, no file is created.
        - Ensures the output directory exists.
    """
    proteins = sorted(feature_map.keys())
    if not proteins:
        return

    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    with open(outfile, "w", newline="", encoding="utf-8") as out_f:
        writer = csv.writer(out_f, delimiter="\t")
        writer.writerow([header_label] + proteins)
        for p1 in proteins:
            row = [p1] + [len(feature_map[p1].intersection(feature_map[p2]))
                          for p2 in proteins]
            writer.writerow(row)


def summarize_overlap(dirs: List[str], outfile: str) -> None:
    """
    Summarize domain or feature overlaps across multiple comparison directories.

    This function loads domain/feature sets using `load_feature_sets()`
    and then writes an overlap matrix using `write_overlap_matrix()`.

    Args:
        dirs (List[str]): List of directories containing TSV comparison files.
        outfile (str): Path to the output TSV file summarizing overlaps.

    Returns:
        None
    """
    data = load_feature_sets(dirs, mode="domains")
    write_overlap_matrix(data, outfile, header_label="Domains")


def main(input_dirs_comma: str, output_file: str) -> None:
    """
    Entry point for summarizing InterProScan domain overlaps from the command line.

    Args:
        input_dirs_comma (str): Comma-separated string of directories containing TSV comparison files.
        output_file (str): Path to the output TSV file.

    Returns:
        None
    """
    dirs = input_dirs_comma.split(",")
    summarize_overlap(dirs, output_file)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python -m summarize.interpro <input_dirs_comma> <output_file>")
        sys.exit(1)

    input_dirs = sys.argv[1]
    output_file = sys.argv[2]
    main(input_dirs, output_file)
