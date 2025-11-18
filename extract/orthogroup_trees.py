"""
Extract target orthogroups from a resolved gene trees file.

This script reads a list of orthogroup IDs from a TSV file and extracts
their corresponding gene trees from a resolved tree file in Newick format.
Each extracted tree is saved as an individual text file in the specified
output directory.

Usage:
    python -m extract.orthogroup_trees <tsv_file> -r <resolved_trees_file> -o <output_dir>
"""

import argparse
from pathlib import Path
import csv
import sys
from typing import List, Dict


def load_orthogroups(tsv_file: Path) -> List[str]:
    """
    Load orthogroup IDs from a TSV file.

    Parameters
    ----------
    tsv_file : Path
        Path to the TSV file containing orthogroup IDs (one per line).

    Returns
    -------
    List[str]
        List of orthogroup IDs.

    Raises
    ------
    FileNotFoundError
        If the TSV file does not exist.
    ValueError
        If no orthogroup IDs are found in the TSV file.
    """
    if not tsv_file.is_file():
        raise FileNotFoundError(f"TSV file not found: {tsv_file}")

    orthogroups: List[str] = []
    try:
        with open(tsv_file, "r", newline="", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                if row and row[0].strip():
                    orthogroups.append(row[0].strip())
    except OSError as e:
        raise RuntimeError(f"Failed to read TSV file: {e}")

    if not orthogroups:
        raise ValueError(f"No orthogroup IDs found in TSV file: {tsv_file}")

    return orthogroups


def extract_trees(resolved_tree_file: Path, target_ogs: List[str]) -> Dict[str, str]:
    """
    Extract trees for target orthogroups from a resolved tree file.

    Parameters
    ----------
    resolved_tree_file : Path
        Path to resolved gene trees file in Newick format.
    target_ogs : List[str]
        List of orthogroup IDs to extract.

    Returns
    -------
    Dict[str, str]
        Mapping from orthogroup ID to Newick tree string.

    Raises
    ------
    FileNotFoundError
        If the resolved tree file does not exist.
    RuntimeError
        If the file cannot be read.
    """
    if not resolved_tree_file.is_file():
        raise FileNotFoundError(f"Resolved tree file not found: \
            {resolved_tree_file}")

    extracted: Dict[str, str] = {}
    try:
        with open(resolved_tree_file, "r", encoding="utf-8") as f:
            for line_num, line in enumerate(f, start=1):
                line = line.strip()
                if not line or ":" not in line:
                    continue
                try:
                    og_id, newick = line.split(":", 1)
                except ValueError:
                    print(
                        f"Warning: Skipping malformed line \
                    {line_num} in {resolved_tree_file}",
                        file=sys.stderr,
                    )
                    continue
                og_id = og_id.strip()
                if og_id in target_ogs:
                    extracted[og_id] = newick.strip()
    except OSError as e:
        raise RuntimeError(f"Failed to read resolved tree file: {e}")

    return extracted


def write_trees(trees: Dict[str, str], target_ogs: List[str], output_dir: Path) -> None:
    """
    Write extracted trees to individual text files.

    Parameters
    ----------
    trees : Dict[str, str]
        Mapping from orthogroup ID to Newick tree string.
    target_ogs : List[str]
        List of orthogroup IDs to write (order preserved).
    output_dir : Path
        Directory to write individual tree files.

    Raises
    ------
    RuntimeError
        If a tree file cannot be written.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    for og_id in target_ogs:
        if og_id not in trees:
            continue
        tree_file = output_dir / f"{og_id}.txt"
        try:
            with open(tree_file, "w", encoding="utf-8") as f:
                f.write(trees[og_id] + "\n")
        except OSError as e:
            raise RuntimeError(f"Failed to write tree file {tree_file}: {e}")


def main() -> None:
    """Command-line interface for extracting orthogroup trees."""
    parser = argparse.ArgumentParser(
        description="Extract target orthogroups from a resolved gene trees file."
    )
    parser.add_argument(
        "file", help="Path to input TSV file with orthogroup IDs.")
    parser.add_argument(
        "-r", "--resolved_trees", required=True, help="Resolved gene trees file (Newick format)."
    )
    parser.add_argument(
        "-o", "--output_dir", required=True, help="Directory to write individual tree files."
    )
    args = parser.parse_args()

    try:
        tsv_file = Path(args.file)
        resolved_tree_file = Path(args.resolved_trees)
        output_dir = Path(args.output_dir)

        target_ogs = load_orthogroups(tsv_file)
        trees = extract_trees(resolved_tree_file, target_ogs)

        missing = set(target_ogs) - set(trees.keys())
        if missing:
            print(
                f"Warning: \
                {len(missing)} orthogroups not found in resolved trees file",
                file=sys.stderr,
            )

        write_trees(trees, target_ogs, output_dir)
        print(f"Extracted {len(trees)} trees to {output_dir}")

    except (FileNotFoundError, ValueError, RuntimeError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
