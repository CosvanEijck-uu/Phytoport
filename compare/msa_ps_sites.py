"""phospho_msa_compare.py

Refactored Phosphorylation Site Comparison Tool with MSA support (Option A).

Features
- Reads a TSV of phosphorylation predictions (per-sequence positions).
- Reads an MSA (FASTA) and maps raw sequence positions to MSA column indices.
- Converts sites to MSA coordinates and compares them there (common vs unique).
- Robust error handling, clear docstrings, and a clean main() entrypoint.

Usage
    python phospho_msa_compare.py <predictions.tsv> <msa.fasta> <seq_header1> <seq_header2> \
        [--cut-off 0.5] [-o output.tsv]

Notes
- Sequence header matching uses the substring after the last '|' (same as the TSV parsing
  convention used previously).
- The script avoids external dependencies (simple FASTA parser) so it runs in plain Python.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict
from typing import Dict, List, Tuple

# ---- Types -----------------------------------------------------------------
Site = Tuple[int, str, float]  # (raw_pos or msa_col, residue, score)
MappedSite = Tuple[int, int, str, float]  # (msa_col, raw_pos, residue, score)


# ---- FASTA / MSA utilities -------------------------------------------------
def parse_fasta(filename: str) -> Dict[str, str]:
    """
    Minimal FASTA parser that returns a mapping header -> sequence (raw string,
    including gap characters if present in the file).

    Header canonicalization: we return only the part after the last '|' to match
    the TSV parsing used elsewhere.
    """
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"FASTA/MSA file not found: {filename}")

    records: Dict[str, List[str]] = {}
    header = None
    try:
        with open(filename, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                if line.startswith(">"):
                    header = line[1:].strip()
                    short_header = header.split("|")[-1]
                    if short_header in records:
                        # duplicate header in file
                        raise ValueError(f"Duplicate FASTA \
                            header found: {short_header}")
                    records[short_header] = []
                else:
                    if header is None:
                        raise ValueError("FASTA file does not \
                            start with a header (" > " line)")
                    records[short_header].append(line.strip())
    except OSError as e:
        raise RuntimeError(f"Failed to read FASTA/MSA file: {e}")

    return {h: "".join(seq_parts) for h, seq_parts in records.items()}


def map_seq_to_msa_positions(aligned_seq: str) -> Dict[int, int]:
    """
    Build a mapping from unaligned (raw) sequence 1-based positions -> MSA column 1-based
    indices.

    Gaps in the aligned sequence ('-' or '.') are treated as alignment gaps. The mapping
    increments the raw position counter only for non-gap characters.
    """
    mapping: Dict[int, int] = {}
    raw_pos = 0
    for msa_idx, aa in enumerate(aligned_seq, start=1):
        if aa not in ("-", "."):
            raw_pos += 1
            mapping[raw_pos] = msa_idx
    return mapping


# ---- TSV extraction and conversion ----------------------------------------
def extract_sites(tsv_file: str, headers: List[str]) -> Dict[str, List[Site]]:
    """
    Extract phosphorylation sites for the requested sequence headers from TSV.

    The TSV is expected to have at least four columns per row:
        full_header \t position \t residue \t score

    Header matching uses the substring after the last '|' (consistent with the MSA parser).
    """
    if not os.path.isfile(tsv_file):
        raise FileNotFoundError(f"TSV file not found: {tsv_file}")

    headers_set = set(headers)
    results: Dict[str, List[Site]] = defaultdict(list)

    # infer allowed residues from file name (keeps original behavior)
    filename = os.path.basename(tsv_file)
    if "_Y" in filename:
        allowed_residues = {"Y"}
    elif "_SorT" in filename:
        allowed_residues = {"S", "T"}
    else:
        allowed_residues = {"S", "T", "Y"}

    try:
        with open(tsv_file, "r", encoding="utf-8") as fh:
            reader = csv.reader(fh, delimiter="\t")
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
                    # skip malformed rows
                    continue
                if residue not in allowed_residues:
                    continue
                if short_header in headers_set:
                    results[short_header].append((position, residue, score))
    except OSError as e:
        raise RuntimeError(f"Failed to read TSV file: {e}")

    # If no results for a requested header, return empty lists (caller will handle)
    for h in headers:
        results.setdefault(h, [])

    return results


def convert_sites_to_msa(sites: List[Site], pos_map: Dict[int, int]) -> List[MappedSite]:
    """
    Convert raw-site tuples (raw_pos, residue, score) to mapped site tuples
    (msa_col, raw_pos, residue, score). Sites that cannot be mapped are skipped.
    """
    converted: List[MappedSite] = []
    for raw_pos, residue, score in sites:
        msa_col = pos_map.get(raw_pos)
        if msa_col is None:
            # site cannot be mapped to the MSA (e.g., position in a region trimmed out or
            # inconsistent alignment). We skip but could log or report these if desired.
            continue
        converted.append((msa_col, raw_pos, residue, score))
    return converted


# ---- Filtering & comparison -----------------------------------------------

def filter_sites_by_cutoff(sites: List[Site], cutoff: float) -> List[Site]:
    """Return sites with score >= cutoff (operates on raw-site tuples)."""
    return [s for s in sites if s[2] >= cutoff]


def find_common_sites_msa(
    mapped1: List[MappedSite], mapped2: List[MappedSite]
) -> Tuple[List[Tuple[int, int, int, str, float, float]], List[MappedSite], List[MappedSite]]:
    """
    Identify common and unique phosphorylation sites in MSA coordinates.

    Returns
    -------
    common: list of tuples (msa_col, raw1, raw2, residue, score1, score2)
    unique1: list of mapped sites unique to first sequence [(msa_col, raw, residue, score), ...]
    unique2: list of mapped sites unique to second sequence
    """
    # Build index by (msa_col, residue) for quick matching
    index2 = {(msa_col, residue): (raw_pos, score)
              for (msa_col, raw_pos, residue, score) in mapped2}

    common: List[Tuple[int, int, int, str, float, float]] = []
    unique1: List[MappedSite] = []
    matched_msa_keys = set()

    for msa_col, raw_pos1, residue1, score1 in mapped1:
        key = (msa_col, residue1)
        if key in index2:
            raw_pos2, score2 = index2[key]
            common.append((msa_col, raw_pos1, raw_pos2,
                          residue1, score1, score2))
            matched_msa_keys.add(key)
        else:
            unique1.append((msa_col, raw_pos1, residue1, score1))

    # unique2 are those in mapped2 whose (msa_col, residue) isn't matched
    unique2 = [s for s in mapped2 if (s[0], s[2]) not in matched_msa_keys]

    # sort for readability
    common.sort(key=lambda x: x[0])
    unique1.sort(key=lambda x: x[0])
    unique2.sort(key=lambda x: x[0])

    return common, unique1, unique2


# ---- Output ---------------------------------------------------------------

def save_to_tsv(
    output_file: str,
    header1: str,
    header2: str,
    common: List[Tuple[int, int, int, str, float, float]],
    unique1: List[MappedSite],
    unique2: List[MappedSite],
) -> None:
    """
    Save comparison to TSV. Columns are designed to be informative for both
    MSA-based and raw-position-based inspection.

    For common sites we write:
      Category	MSA_col	Residue	RawPos_header1	Score_header1	RawPos_header2	Score_header2

    For unique sites we write:
      Category	MSA_col	Residue	RawPos	Score
    """
    try:
        with open(output_file, "w", encoding="utf-8", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            # header
            writer.writerow([
                "Category", "MSA_col", "Residue",
                f"{header1}_RawPos", f"{header1}_Score",
                f"{header2}_RawPos", f"{header2}_Score",
            ])

            for msa_col, raw1, raw2, residue, score1, score2 in common:
                writer.writerow(["Common", msa_col, residue,
                                raw1, score1, raw2, score2])

            for msa_col, raw, residue, score in unique1:
                writer.writerow(
                    [f"Unique_to_{header1}", msa_col, residue, raw, score, "", ""])

            for msa_col, raw, residue, score in unique2:
                writer.writerow(
                    [f"Unique_to_{header2}", msa_col, residue, "", "", raw, score])
    except OSError as e:
        raise RuntimeError(f"Failed to write output file: {e}")


# ---- High-level orchestration ---------------------------------------------

def compare_phospho_sites(
    tsv_file: str,
    msa_file: str,
    headers: List[str],
    cutoff: float = 0.0,
    output_file: str = "phospho_comparison.tsv",
) -> None:
    """
    Main comparsion routine with MSA support.

    Steps:
    - parse MSA
    - build position maps for both sequences
    - extract sites from TSV
    - filter by cutoff
    - convert to MSA coordinates
    - find common/unique sites and save results
    """
    if len(headers) != 2:
        raise ValueError("Please provide exactly \
            two sequence headers for comparison.")

    # load MSA
    msa = parse_fasta(msa_file)
    missing = [h for h in headers if h not in msa]
    if missing:
        raise ValueError(f"The following headers \
            were not found in the MSA: {missing}")

    # build raw->msa position maps
    pos_maps = {h: map_seq_to_msa_positions(msa[h]) for h in headers}

    # extract sites from TSV
    site_dict = extract_sites(tsv_file, headers)
    header1, header2 = headers

    if not site_dict[header1] and not site_dict[header2]:
        raise ValueError("No phosphorylation sites \
            found for either header in the TSV.")

    # filter by cutoff
    raw_sites1 = filter_sites_by_cutoff(site_dict.get(header1, []), cutoff)
    raw_sites2 = filter_sites_by_cutoff(site_dict.get(header2, []), cutoff)

    # convert to MSA coordinates
    mapped1 = convert_sites_to_msa(raw_sites1, pos_maps[header1])
    mapped2 = convert_sites_to_msa(raw_sites2, pos_maps[header2])

    if not mapped1 and not mapped2:
        raise ValueError("After mapping to the MSA \
        and applying the cutoff no sites remain.")

    # find common and unique sites
    common, unique1, unique2 = find_common_sites_msa(mapped1, mapped2)

    # save
    save_to_tsv(output_file, header1, header2, common, unique1, unique2)

    print(f"Phosphorylation comparison (MSA-based) saved to: {output_file}")


# ---- CLI Entrypoint ------------------------------------------------------

def main(argv: List[str] | None = None) -> int:
    """
    CLI wrapper. Returns exit code (0 == success).
    """
    parser = argparse.ArgumentParser(
        description="Compare predicted phosphorylation sites between two sequences using an MSA (alignment coordinates)."
    )
    parser.add_argument(
        "tsv_file", help="Input TSV file of phosphorylation predictions.")
    parser.add_argument(
        "msa_file", help="MSA file in FASTA format (aligned sequences).")
    parser.add_argument("sequence_headers", nargs=2,
                        help="Two sequence headers to compare (use the part after the last '|').")
    parser.add_argument("--cut-off", type=float, default=0.0,
                        help="Minimum prediction score to include a site (default: 0.0).")
    parser.add_argument("-o", "--output", default="phospho_comparison.tsv",
                        help="Output TSV file (default: phospho_comparison.tsv).")

    args = parser.parse_args(argv)

    try:
        compare_phospho_sites(
            args.tsv_file,
            args.msa_file,
            args.sequence_headers,
            cutoff=args.cut_off,
            output_file=args.output,
        )
        return 0
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    except Exception as exc:  # defensive: unexpected errors
        print(f"Unexpected error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
