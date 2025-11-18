"""
Protein Sequence Visualizer
============================

Visualizes one or more protein sequences from a FASTA file, using per-residue
scores stored in a TSV file. Each residue is drawn as a colored square with
font size scaling according to the score.

Features
--------
- Supports multiple sequences in FASTA.
- TSV may contain scores for multiple identifiers.
- Residues without scores are drawn in a default gray.
- Font size scales with the score.
- Automatic row-wrapping for long sequences.
- Sequence label shown above each sequence (based on final `|` section).
- Colorbar legend included beneath the visualization, width matches sequence.
- Optional CLI subset selection based on sequence label.
"""

import argparse
import sys
import os
from typing import Dict, List

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors

os.environ.pop("QT_STYLE_OVERRIDE", None)


def read_fasta(file_path: str) -> Dict[str, str]:
    """Read a FASTA file and return sequences in a dict."""
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"FASTA file not found: {file_path}")

    sequences = {}
    seq_id = None
    seq_lines: List[str] = []

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_id:
                    sequences[seq_id] = "".join(seq_lines).upper()
                seq_id = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if seq_id:
            sequences[seq_id] = "".join(seq_lines).upper()
    return sequences


def read_tsv(tsv_path: str) -> Dict[str, Dict[int, float]]:
    """Read TSV file with columns: seq_id, position, ..., score."""
    if not os.path.isfile(tsv_path):
        raise FileNotFoundError(f"TSV file not found: {tsv_path}")

    score_dict: Dict[str, Dict[int, float]] = {}

    with open(tsv_path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            seq_id = parts[0].strip().strip('"').lstrip(">")
            try:
                pos = int(parts[1]) - 1
                score = float(parts[3])
            except ValueError:
                continue
            score_dict.setdefault(seq_id, {})[pos] = score

    return score_dict


def wrap_sequences(sequences: Dict[str, str], row_length: int) -> Dict[str, List[str]]:
    """Split sequences into chunks of fixed length."""
    return {
        seq_id: [seq[i:i + row_length] for i in range(0, len(seq), row_length)]
        for seq_id, seq in sequences.items()
    }


def compute_total_rows(wrapped: Dict[str, List[str]], separator_height: float) -> float:
    """Compute total vertical rows including separator space."""
    total_rows = sum(len(rows) for rows in wrapped.values())
    total_rows += separator_height * (len(wrapped) - 1)
    return total_rows


def create_color_norm(score_dict: Dict[str, Dict[int, float]]):
    """Create colormap normalization based on all scores."""
    all_scores = [score for seq_scores in score_dict.values()
                  for score in seq_scores.values()]
    if not all_scores:
        raise ValueError("No scores found in TSV input.")
    vmin = min(all_scores)
    vmax = max(all_scores)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.colormaps.get_cmap("coolwarm")
    return cmap, norm, vmin, vmax


def score_to_font(score: float, min_score: float, max_score: float,
                  min_font: int = 6, max_font: int = 16,
                  font_increase: float = 4.0) -> float:
    """Convert score to font size with linear scaling."""
    if max_score == min_score:
        return (min_font + max_font) / 2 + font_increase
    return min_font + (score - min_score) / (max_score - min_score) * (max_font - min_font) + font_increase


def draw_sequences(
    ax,
    wrapped_sequences: Dict[str, List[str]],
    score_dict: Dict[str, Dict[int, float]],
    total_rows: float,
    row_length: int,
    separator_height: float = 0.3,
    min_font: int = 6,
    max_font: int = 16,
    font_increase: float = 4,
    label_font: int = 10,
    number_font: int = 6,
    default_color: str = "#FFFFFF",
):
    """Draw sequences, residues, numbers, labels, separators, and score-based colors."""
    cmap, norm, min_score, max_score = create_color_norm(score_dict)
    y_offset = 0

    for seq_idx, (seq_id, rows) in enumerate(wrapped_sequences.items()):
        label = seq_id.split("|")[-1]
        first_row_y = total_rows - y_offset - 0.5
        ax.text(-1, first_row_y, label, ha="right", va="center",
                fontsize=label_font, weight="bold")

        for row_idx, chunk in enumerate(rows):
            for i, aa in enumerate(chunk):
                seq_pos = row_idx * row_length + i
                x = i
                y = total_rows - y_offset - 1

                if seq_pos in score_dict.get(seq_id, {}):
                    score = score_dict[seq_id][seq_pos]
                    color = cmap(norm(score))
                    fontsize = score_to_font(
                        score, min_score, max_score, min_font, max_font, font_increase)
                else:
                    color = default_color
                    fontsize = min_font + font_increase

                ax.add_patch(plt.Rectangle((x, y), 1, 1, color=color))
                ax.text(x + 0.5, y + 0.5, aa, ha="center", va="center",
                        fontsize=fontsize, color="black")
                ax.text(x + 0.5, y + 0.4, str(seq_pos + 1), ha="center", va="top",
                        fontsize=number_font, color="black", alpha=0.7)

            y_offset += 1

        if seq_idx < len(wrapped_sequences) - 1:
            sep_y = total_rows - y_offset - 0.5 * separator_height
            ax.hlines(y=sep_y, xmin=0, xmax=row_length,
                      color="black", linewidth=1)
            y_offset += separator_height


def add_colorbar(fig, ax, cmap, norm, row_length: int):
    """Add a horizontal colorbar beneath the sequences with width matching the sequence row."""
    ax_pos = ax.get_position()
    cbar_height = 0.02
    cbar_bottom = 0.02
    cbar_left = ax_pos.x0
    cbar_width = row_length / row_length * ax_pos.width  # full row width
    cax = fig.add_axes([cbar_left, cbar_bottom, ax_pos.width, cbar_height])

    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cbar.set_label("Score", fontsize=12)


def visualize_sequences(
    sequences: Dict[str, str],
    score_dict: Dict[str, Dict[int, float]],
    output_file: str = "protein.png",
    min_font: int = 6,
    max_font: int = 20,
    font_increase: float = 8,
    row_length: int = 100,
    separator_height: float = 0.3,
    label_font: int = 10,
):
    """Generate the full visualization."""
    if not sequences:
        raise ValueError("No sequences available to visualize.")

    wrapped = wrap_sequences(sequences, row_length)
    total_rows = compute_total_rows(wrapped, separator_height)

    fig_height = total_rows * 1.5
    fig_width = max(row_length / 5, 10)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    draw_sequences(ax, wrapped, score_dict, total_rows, row_length,
                   separator_height, min_font, max_font, font_increase, label_font)

    ax.set_xlim(-5, row_length)
    ax.set_ylim(0, total_rows)
    ax.axis("off")
    ax.set_position([0, 0.05, 1, 0.95])

    cmap, norm, _, _ = create_color_norm(score_dict)
    add_colorbar(fig, ax, cmap, norm, row_length=row_length)

    plt.savefig(output_file, dpi=300, transparent=True, bbox_inches="tight")
    plt.close()
    print(f"Protein visualization saved as '{output_file}'.")


def main():
    parser = argparse.ArgumentParser(
        description="Visualize multiple protein sequences with per-residue TSV scores."
    )
    parser.add_argument("fasta_file", help="Input FASTA file")
    parser.add_argument("-t", "--tsv_file", required=True,
                        help="TSV file containing per-residue scores")
    parser.add_argument("-o", "--output", default="protein.png",
                        help="Output PNG file")
    parser.add_argument("--row_length", type=int, default=100,
                        help="Residues per display row")
    parser.add_argument("--label_font", type=int, default=10,
                        help="Font size for sequence labels")
    parser.add_argument("-s", "--subset", type=str,
                        help="Comma-separated list of sequence labels to visualize")

    args = parser.parse_args()

    try:
        sequences = read_fasta(args.fasta_file)
        score_dict = read_tsv(args.tsv_file)

        if args.subset:
            subset_ids = set(args.subset.split(","))
            sequences = {
                seq_id: seq
                for seq_id, seq in sequences.items()
                if seq_id.split("|")[-1] in subset_ids
            }

        visualize_sequences(sequences, score_dict,
                            output_file=args.output,
                            row_length=args.row_length,
                            label_font=args.label_font)

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
