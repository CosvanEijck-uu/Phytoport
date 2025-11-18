import requests
import json
import time
import argparse
import sys
from requests.exceptions import RequestException

API_SUBMIT_URL = "https://alignment.rcsb.org/api/v1/structures/submit"
API_RESULTS_URL = "https://alignment.rcsb.org/api/v1/structures/results"


def submit_alignment(file1: str, file2: str, method: str = "fatcat-flexible") -> str:
    """
    Submit two CIF/mmCIF files for pairwise structure alignment using RCSB FATCAT.

    Parameters:
        file1 (str): Path to the first CIF/mmCIF file.
        file2 (str): Path to the second CIF/mmCIF file.
        method (str): Alignment method ("fatcat-flexible" or "fatcat-rigid").

    Returns:
        str: Job ticket UUID.

    Raises:
        RuntimeError: If the submission fails.
    """
    query = {
        "context": {
            "mode": "pairwise",
            "method": {"name": method},
            "structures": [
                {"format": "mmcif"},
                {"format": "mmcif"}
            ]
        }
    }

    data = {"query": json.dumps(query)}

    try:
        files = [
            ("files", (file1, open(file1, "rb"))),
            ("files", (file2, open(file2, "rb")))
        ]
    except FileNotFoundError as e:
        raise RuntimeError(f"File not found: {e.filename}")

    print(f"Submitting alignment job ({method})...")
    try:
        response = requests.post(API_SUBMIT_URL, params=data, files=files)
        response.raise_for_status()
    except RequestException as e:
        raise RuntimeError(f"Failed to submit alignment job: {e}")
    finally:
        # Close file handles
        for _, (_, f) in files:
            f.close()

    ticket = response.text.strip()
    if not ticket:
        raise RuntimeError(
            f"Failed to obtain ticket. Response: {response.text}")

    print(f"Job submitted! Ticket: {ticket}")
    return ticket


def wait_for_results(ticket: str, poll_interval: int = 5) -> dict:
    """
    Poll the RCSB API until the alignment job is complete.

    Parameters:
        ticket (str): Job ticket UUID.
        poll_interval (int): Time in seconds between polling attempts.

    Returns:
        dict: Full JSON results from the alignment job.

    Raises:
        RuntimeError: If the job returns an error or unexpected response.
    """
    url = f"{API_RESULTS_URL}?uuid={ticket}"

    while True:
        try:
            response = requests.get(url)
            response.raise_for_status()
        except RequestException as e:
            raise RuntimeError(f"Failed to fetch results: {e}")

        try:
            result = response.json()
        except json.JSONDecodeError:
            text = getattr(response, "text", "<no response>")
            raise RuntimeError(f"Invalid JSON response: {text}")

        status = result.get("info", {}).get("status")
        if status == "RUNNING":
            print("Job still running... waiting...")
            time.sleep(poll_interval)
        elif status == "COMPLETE":
            print("Job complete!")
            return result
        elif status == "ERROR":
            message = result.get("info", {}).get("message", "No message")
            raise RuntimeError(f"Alignment job error: {message}")
        else:
            raise RuntimeError(f"Unexpected response: {result}")


def print_summary(results: dict):
    """
    Print concise summary statistics of the alignment.

    Parameters:
        results (dict): JSON results from the RCSB API.
    """
    if not results.get("results"):
        print("No alignment results found.")
        return

    summary = results["results"][0].get("summary", {})
    print("\n=== Alignment Summary ===")
    for score in summary.get("scores", []):
        print(f"{score['type']}: {score['value']}")
    print(f"Aligned residue pairs: \
        {summary.get('n_aln_residue_pairs', 'N/A')}")
    print(f"Modeled residues per structure: \
        {summary.get('n_modeled_residues', 'N/A')}")


def main():
    parser = argparse.ArgumentParser(
        description="Submit two CIF/mmCIF files for pairwise FATCAT alignment and summarize results."
    )
    parser.add_argument("file1", help="First CIF/mmCIF file")
    parser.add_argument("file2", help="Second CIF/mmCIF file")
    parser.add_argument(
        "--method",
        default="fatcat-flexible",
        choices=["fatcat-flexible", "fatcat-rigid", "tm-align"],
        help="Alignment method (default: fatcat-flexible)"
    )
    parser.add_argument(
        "--out",
        default="alignment_results.json",
        help="File to save full JSON results (default: alignment_results.json)"
    )
    args = parser.parse_args()

    try:
        ticket = submit_alignment(args.file1, args.file2, args.method)
        results = wait_for_results(ticket)
    except RuntimeError as e:
        print(f"Error: {e}")
        sys.exit(1)

    # Save full JSON results
    try:
        with open(args.out, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\nFull results written to {args.out}")
    except IOError as e:
        print(f"Failed to write results to {args.out}: {e}")

    # Print summary statistics
    print_summary(results)


if __name__ == "__main__":
    main()
