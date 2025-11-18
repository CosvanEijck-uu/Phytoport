"""
Fetch the first instance of barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
from GEO (GSE or GSM) FTP directories, with an HTTPS fallback if FTP is blocked.

Usage:
    python -m fetch.geo_dataset GSM6529487
    python -m fetch.geo_dataset GSE123456 -o ./data/geo

Notes:
    - Tries FTP first with a 10-second timeout.
    - If FTP fails, automatically falls back to HTTPS.
    - Only downloads the first occurrence of each 10X-style file.
    - Saves files in the output directory under a subfolder named by the GEO ID.
"""

import sys
import os
import re
import ftplib
import urllib.request
import argparse
from typing import Optional


def find_geo_ftp_path(geo_id: str) -> str:
    """
    Return the relative GEO FTP directory path for the given GEO accession.

    Parameters
    ----------
    geo_id : str
        GEO identifier starting with GSE or GSM.

    Returns
    -------
    str
        Relative FTP path to the GEO supplementary files.

    Raises
    ------
    ValueError
        If the identifier does not start with GSE or GSM.
    """
    geo_num = re.findall(r'\d+', geo_id)[0]
    prefix = geo_num[:-3] + "nnn" if len(geo_num) > 3 else "0nnn"

    if geo_id.startswith("GSE"):
        return f"/geo/series/GSE{prefix}/{geo_id}/suppl/"
    elif geo_id.startswith("GSM"):
        return f"/geo/samples/GSM{prefix}/{geo_id}/suppl/"
    else:
        raise ValueError("Identifier must start with GSE or GSM.")


def fetch_first_matrix_files(geo_id: str, base_outdir: str) -> None:
    """
    Fetch the first occurrence of 10X-style GEO files (barcodes, features, matrix)
    using FTP first and falling back to HTTPS if FTP is blocked.

    Parameters
    ----------
    geo_id : str
        GEO accession (GSE or GSM).
    base_outdir : str
        Base output directory where the GEO subfolder will be created.

    Raises
    ------
    RuntimeError
        If neither FTP nor HTTPS can access the GEO folder.
    """
    outdir = os.path.abspath(os.path.join(base_outdir, geo_id))
    os.makedirs(outdir, exist_ok=True)

    ftp_server = "ftp.ncbi.nlm.nih.gov"
    ftp_path = find_geo_ftp_path(geo_id)
    files = []

    # --- Attempt FTP first ---
    try:
        ftp = ftplib.FTP(ftp_server, timeout=10)  # 10-second timeout
        ftp.login()
        ftp.cwd(ftp_path)
        files = ftp.nlst()  # list all files in directory
        ftp.quit()
        print(f"✅ FTP connection succeeded for {geo_id}")
        base_url = f"ftp://{ftp_server}{ftp_path}"
    except Exception as e:
        print(f"⚠️ FTP failed ({e}), falling back to HTTPS...")
        # Construct HTTPS URL for the same GEO folder
        base_url = f"https://ftp.ncbi.nlm.nih.gov{ftp_path}"
        try:
            import requests
            resp = requests.get(base_url)
            if resp.status_code != 200:
                raise RuntimeError(f"HTTPS access failed with status \
                    {resp.status_code}")
            # Extract file names from HTML directory listing
            files = re.findall(r'href="([^"]+)"', resp.text)
        except Exception as e2:
            raise RuntimeError(f"Cannot access GEO via FTP or HTTPS: {e2}")

    # --- Find target 10X files ---
    targets: dict[str, Optional[str]] = {
        "barcodes": None, "features": None, "matrix": None
    }

    for f in files:
        lower_f = f.lower()
        if targets["barcodes"] is None and "barcodes.tsv.gz" in lower_f:
            targets["barcodes"] = f
        elif targets["features"] is None and "features.tsv.gz" in lower_f:
            targets["features"] = f
        elif targets["matrix"] is None and "matrix.mtx.gz" in lower_f:
            targets["matrix"] = f
        if all(targets.values()):
            break

    # --- Download the files ---
    for key, filename in targets.items():
        if filename:
            dest = os.path.join(
                outdir, f"{key}.tsv.gz" if key != "matrix" else "matrix.mtx.gz")
            url = f"{base_url}{filename}" if base_url.endswith(
                "/") else f"{base_url}/{filename}"
            print(f"📥 Downloading {filename} → {dest}")
            try:
                urllib.request.urlretrieve(url, dest)
            except Exception as e:
                print(f"⚠️ Failed to download {filename}: {e}")
        else:
            print(f"⚠️ No {key} file found for {geo_id}")

    print(f"\n✅ Done. Files saved in: {outdir}")


def main(geo_id: str, outdir: str = ".") -> None:
    """
    Main logic to fetch GEO 10X files.

    Parameters
    ----------
    geo_id : str
        GEO accession (GSE or GSM)
    outdir : str
        Base output directory
    """
    if not re.match(r"^(GSM|GSE)\d+$", geo_id):
        raise ValueError(
            "Invalid GEO ID. Must look like GSM1234567 or GSE123456."
        )
    fetch_first_matrix_files(geo_id, outdir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fetch GEO 10X-style files (GSE or GSM).")
    parser.add_argument(
        "geo_id",
        help="GSE or GSM identifier (e.g., GSM6529487 or GSE123456)"
    )
    parser.add_argument(
        "-o", "--outdir",
        default=".",
        help="Base output directory to save GEO folders (default: current directory)"
    )
    args = parser.parse_args()

    try:
        main(args.geo_id, args.outdir)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
