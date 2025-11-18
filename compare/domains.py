"""
Domain Comparison Tool

This script compares predicted domains or motifs between two protein sequences
based on a TSV file of domain predictions. It identifies common and unique
domains, computes delta start/end distances for common domains, and outputs
a detailed TSV report.

Usage:
    python -m compare.domains <tsv_file> <seq1> <seq2> [-o output.tsv]
"""

import argparse
import csv
import os
import sys
from collections import defaultdict
from typing import List, Tuple

Domain = Tuple[str, str, int, int]  # (predictor, domain_name, start, end)


def extract_headers(tsv_file: str, headers: List[str]) -> dict[str, List[Domain]]:
    """
    Extract domain data from TSV for given sequence headers.

    Parameters
    ----------
    tsv_file : str
        Path to TSV file containing domain predictions.
    headers : list[str]
        Sequence headers to extract domains for (after last '|').

    Returns
    -------
    dict[str, list[Domain]]
        Mapping from header to list of domains (predictor, domain, start, end).

    Raises
    ------
    FileNotFoundError
        If TSV file does not exist.
    ValueError
        If no domains are found for the provided headers.
    """
    if not os.path.isfile(tsv_file):
        raise FileNotFoundError(f"TSV file not found: {tsv_file}")

    results: dict[str, List[Domain]] = defaultdict(list)
    headers_set = set(headers)

    try:
        with open(tsv_file, "r", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                if len(row) < 8:
                    continue
                full_header = row[0]
                short_header = full_header.split("|")[-1]
                if short_header not in headers_set:
                    continue
                predictor, domain = row[3], row[5]
                try:
                    start, end = int(row[6]), int(row[7])
                except ValueError:
                    continue
                results[short_header].append((predictor, domain, start, end))
    except OSError as e:
        raise RuntimeError(f"Failed to read TSV file: {e}")

    if not results:
        raise ValueError(f"No domains found for headers: {headers}")

    return results


def find_common_domains(domains1: List[Domain], domains2: List[Domain]):
    """
    Find common and unique domains between two sequences.

    Domains are matched by domain name; for common domains, the nearest
    start position is used to pair.

    Parameters
    ----------
    domains1, domains2 : list[Domain]
        Lists of domains for the two sequences.

    Returns
    -------
    tuple
        (common_domains, unique_to_first, unique_to_second)
    """
    by_name2 = defaultdict(list)
    for d in domains2:
        by_name2[d[1]].append(d)

    common, unique1 = [], []
    used_d2 = set()

    for d1 in domains1:
        src1, name1, s1, e1 = d1
        if name1 not in by_name2:
            unique1.append(d1)
            continue

        candidates = [d2 for d2 in by_name2[name1] if d2 not in used_d2]
        if not candidates:
            unique1.append(d1)
            continue

        # Match nearest start position
        d2 = min(candidates, key=lambda x: abs(x[2] - s1))
        used_d2.add(d2)
        common.append((name1, (s1, e1), (d2[2], d2[3]), src1, d2[0]))

    unique2 = [d for dlist in by_name2.values()
               for d in dlist if d not in used_d2]
    return common, unique1, unique2


def summarize_distances(common: List[Tuple]) -> dict[str, dict[str, int]]:
    """
    Compute maximum Δstart and Δend per domain name for common domains.

    Parameters
    ----------
    common : list
        List of paired common domains.

    Returns
    -------
    dict
        Mapping from domain name to {"Δstart": int, "Δend": int}.
    """
    summary = defaultdict(lambda: {"Δstart": 0, "Δend": 0})
    for name, (s1, e1), (s2, e2), *_ in common:
        ds, de = abs(s1 - s2), abs(e1 - e2)
        summary[name]["Δstart"] = max(summary[name]["Δstart"], ds)
        summary[name]["Δend"] = max(summary[name]["Δend"], de)
    return summary


def save_to_tsv(
    output_file: str,
    header1: str,
    header2: str,
    common: List,
    unique1: List[Domain],
    unique2: List[Domain],
    summary: dict[str, dict[str, int]]
) -> None:
    """
    Save domain comparison results to a TSV file.

    Parameters
    ----------
    output_file : str
        Path to output TSV file.
    header1, header2 : str
        Sequence headers compared.
    common : list
        Common domains.
    unique1, unique2 : list[Domain]
        Unique domains for each sequence.
    summary : dict
        Summary of Δstart/Δend for common domains.

    Raises
    ------
    RuntimeError
        If writing the TSV file fails.
    """
    try:
        with open(output_file, "w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow([
                "Category", "Domain", "Protein1_Source", "Protein1_Start", "Protein1_End",
                "Protein2_Source", "Protein2_Start", "Protein2_End", "Delta_Start", "Delta_End"
            ])

            for name, (s1, e1), (s2, e2), src1, src2 in common:
                writer.writerow(["Common", name, src1, s1, e1,
                                src2, s2, e2, abs(s1 - s2), abs(e1 - e2)])

            for predictor, domain, start, end in unique1:
                writer.writerow(
                    [f"Unique_to_{header1}", domain, predictor, start, end, "", "", "", "", ""])
            for predictor, domain, start, end in unique2:
                writer.writerow(
                    [f"Unique_to_{header2}", domain, "", "", "", predictor, start, end, "", ""])

            for name, vals in summary.items():
                writer.writerow(["Summary", name, "", "", "",
                                "", "", "", vals["Δstart"], vals["Δend"]])
    except OSError as e:
        raise RuntimeError(f"Failed to write output file: {e}")


def compare_domains(tsv_file: str, headers: List[str], output_file: str) -> None:
    """
    Compare domains between two protein sequences and save results.

    Parameters
    ----------
    tsv_file : str
        Input TSV file of domain predictions.
    headers : list[str]
        Two sequence headers to compare.
    output_file : str
        Output TSV file path.
    """
    if len(headers) != 2:
        raise ValueError("Provide exactly two sequence headers.")

    domain_dict = extract_headers(tsv_file, headers)
    header1, header2 = headers
    domains1 = domain_dict.get(header1, [])
    domains2 = domain_dict.get(header2, [])

    common, only1, only2 = find_common_domains(domains1, domains2)
    summary = summarize_distances(common)

    save_to_tsv(output_file, header1, header2, common, only1, only2, summary)
    print(f"Comparison results saved to {output_file}")


def main() -> None:
    """Command-line interface for domain comparison."""
    parser = argparse.ArgumentParser(
        description="Compare domains/motifs between two protein headers.")
    parser.add_argument(
        "tsv_file", help="Input TSV file of domain predictions.")
    parser.add_argument("sequence_headers", nargs=2,
                        help="Two sequence headers to compare (after last '|').")
    parser.add_argument(
        "-o", "--output", default="comparison_results.tsv", help="Output TSV file.")
    args = parser.parse_args()

    try:
        compare_domains(args.tsv_file, args.sequence_headers, args.output)
    except (FileNotFoundError, ValueError, RuntimeError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
