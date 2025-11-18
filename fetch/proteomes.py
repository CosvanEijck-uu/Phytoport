"""
fetch/proteomes.py

Download reference proteomes from UniProt by proteome UPID.

Usage:
    python -m fetch.proteomes --upid UP000006548 UP000005640 --output_dir ./data
"""

import os
import gzip
import argparse
import shutil
import requests
from typing import Optional, List
from tqdm import tqdm


def download_gz_fasta(proteome_upid: str, output_path: str) -> None:
    """
    Download a compressed FASTA from UniProt and extract it.

    Parameters:
        proteome_upid (str): UniProt reference proteome ID.
        output_path (str): Path where the extracted FASTA file will be saved.

    Raises:
        OSError: If the file cannot be written or extracted.
        requests.RequestException: On download errors.
    """
    gz_path = output_path + ".gz"
    url = (
        f"https://rest.uniprot.org/uniprotkb/stream?"
        f"compressed=true&format=fasta&query=proteome: \
    {proteome_upid}&is_isoform:true"
    )

    if os.path.exists(output_path):
        print(f"{output_path} already exists. Skipping download.")
        return

    try:
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            total_size = int(r.headers.get("content-length", 0))
            with open(gz_path, "wb") as f, tqdm(
                total=total_size, unit="B", unit_scale=True, desc=f"Downloading {gz_path}"
            ) as t:
                for chunk in r.iter_content(1024):
                    f.write(chunk)
                    t.update(len(chunk))

        with gzip.open(gz_path, "rb") as f_in, open(output_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        os.remove(gz_path)
        print(f"Saved {output_path} and removed {gz_path}")

    except requests.RequestException as e:
        print(f"Failed to download {proteome_upid}: {e}")
    except OSError as e:
        print(f"Failed to extract {gz_path}: {e}")


def download_reference_proteome_for_upids(
    upids: List[str], output_dir: str, prefixes: Optional[List[str]] = None
) -> None:
    """
    Download reference proteomes for a list of UniProt proteome UPIDs.

    Parameters:
        upids (List[str]): List of UniProt proteome IDs.
        output_dir (str): Directory to save the FASTA files.
        prefixes (Optional[List[str]]): Optional prefixes for output filenames
            (matched by order to upids).
    """
    for i, upid in enumerate(upids):
        prefix = prefixes[i] if prefixes and i < len(prefixes) else upid
        output_path = os.path.join(output_dir, f"{prefix}_proteome.fasta")
        download_gz_fasta(upid, output_path)


# ----------------------
# Main
# ----------------------
def main(
    upids: Optional[List[str]],
    prefixes: Optional[List[str]],
    output_dir: str
) -> None:
    """
    Main logic to download reference proteomes by UPID.

    Parameters:
        upids (Optional[List[str]]): List of UniProt proteome UPIDs.
        prefixes (Optional[List[str]]): Optional filename prefixes for UPIDs.
        output_dir (str): Directory to save FASTA files.

    Raises:
        ValueError: If no UPIDs are provided.
    """
    os.makedirs(output_dir, exist_ok=True)

    if not upids:
        raise ValueError(
            "You must provide at least one UniProt proteome UPID.")

    download_reference_proteome_for_upids(upids, output_dir, prefixes)


# ----------------------
# CLI entry point
# ----------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Download reference proteomes from UniProt by UPID."
    )
    parser.add_argument("--upid", nargs="+", required=True,
                        help="UniProt proteome UPIDs (e.g., UP000006548)")
    parser.add_argument("--prefix", nargs="+",
                        help="Optional prefixes for output filenames (matched by order to --upid)")
    parser.add_argument("--output_dir", default=".",
                        help="Directory to save downloaded FASTA files")

    args = parser.parse_args()
    try:
        main(args.upid, args.prefix, args.output_dir)
    except ValueError as e:
        parser.error(str(e))
