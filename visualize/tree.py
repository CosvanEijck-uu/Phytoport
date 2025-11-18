"""
Render Newick trees with optional condensed MSA using ETE3.

This script renders phylogenetic trees with optional domain headers and
sequence visualization. It supports an automatic, robust cap on branch
lengths to avoid skew from outliers, which can be overridden via
`--max-branch-length`.
"""

import argparse
from typing import Optional, List, Tuple
from ete3 import PhyloTree, TreeStyle, faces

# Okabe–Ito colorblind-safe palette (hex)
_COLORBLIND_PALETTE = [
    "#000000",  # black
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # bluish green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
]

_domain_colors = {}  # ensure this exists globally


def domain_to_color(domain: str) -> str:
    """
    Assign a deterministic color-blind friendly color to a domain.

    Cycles through the Okabe–Ito palette deterministically using the count
    of previously-seen domains.

    Parameters
    ----------
    domain : str
        Domain identifier.

    Returns
    -------
    str
        Hex color code.
    """
    if domain not in _domain_colors:
        index = len(_domain_colors) % len(_COLORBLIND_PALETTE)
        _domain_colors[domain] = _COLORBLIND_PALETTE[index]
    return _domain_colors[domain]


def aa_to_color(aa: str) -> str:
    """
    Map an amino-acid character to a color.

    Unknown characters map to a light gray.
    """
    colors = {
        "A": "#009E73",
        "C": "#0072B2",
        "D": "#D55E00",
        "E": "#E69F00",
        "F": "#56B4E9",
        "G": "#000000",
        "H": "#CC79A7",
        "I": "#999999",
        "K": "#009E73",
        "L": "#0072B2",
        "M": "#56B4E9",
        "N": "#F0E442",
        "P": "#CC79A7",
        "Q": "#E69F00",
        "R": "#000000",
        "S": "#56B4E9",
        "T": "#009E73",
        "V": "#0072B2",
        "W": "#D55E00",
        "Y": "#CC79A7",
        "-": "#999999",
    }
    return colors.get(aa.upper(), "#BBBBBB")


def parse_leaf_name(name: Optional[str]) -> str:
    """Simplify a leaf name by splitting on '|' and returning the last part."""
    if name is None:
        return ""
    return name.split("|")[-1]


# -------------------
# Argument parsing
# -------------------


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Render Newick trees with optional condensed MSA using ETE3."
    )
    parser.add_argument("tree_file", help="Input Newick tree file")
    parser.add_argument("-a", "--alignment",
                        help="Optional sequence alignment file")
    parser.add_argument(
        "-o",
        "--output",
        default="tree.png",
        help="Output image file (PNG, PDF, SVG, etc.)",
    )
    parser.add_argument(
        "--domains",
        help="Optional TSV file with domain annotations to display in header",
    )
    parser.add_argument(
        "-w", "--width", type=int, default=250, help="Output image width"
    )
    parser.add_argument(
        "--reroot", help="Optional leaf name to reroot the tree on")
    parser.add_argument(
        "--no-branch-lengths",
        action="store_true",
        help="Hide branch lengths in the plot",
    )
    parser.add_argument(
        "--max-branch-length",
        type=float,
        help="Prune nodes with branch length greater than this value",
    )
    parser.add_argument(
        "--horizontal-scale", type=float, help="Set horizontal scale of the tree"
    )
    parser.add_argument(
        "--max-domains",
        type=int,
        default=5,
        help="Maximum number of longest domains (per unique source) to display",
    )
    return parser.parse_args()


# -------------------
# File loaders
# -------------------


def load_domains(
    tsv_file: str,
    max_domains: int = 5
) -> List[Tuple[int, int, str, str]]:
    """
    Load domain annotations from a TSV file and keep only the longest
    domain per unique source. Domains named '-' are ignored.

    Parameters
    ----------
    tsv_file : str
        Path to the TSV file.
    max_domains : int, optional
        Maximum number of longest domains to return (default: 5).

    Returns
    -------
    list of (start, end, name, source)
        Sorted by domain length (longest first). Returns empty list on error.
    """
    source_to_domains = {}

    try:
        with open(tsv_file, "r") as f:
            for line_no, line in enumerate(f, start=1):
                if not line.strip() or line.startswith("#"):
                    continue

                row = line.rstrip("\n").split("\t")
                if len(row) < 9:
                    # ignore malformed rows but warn
                    print(f"WARNING: Skipping malformed row \
                    {line_no} in {tsv_file}")
                    continue

                source = row[3]
                name = row[5]

                # Skip domains with invalid or placeholder names
                if name.strip() == "-":
                    continue

                try:
                    start = int(row[6])
                    end = int(row[7])
                except ValueError:
                    print(
                        f"WARNING: Invalid start/end on row {line_no} in {tsv_file}; skipping")
                    continue

                # ensure start <= end
                if start > end:
                    start, end = end, start

                domain = (start, end, name, source)
                source_to_domains.setdefault(source, set()).add(domain)

    except FileNotFoundError:
        print(f"ERROR: Domain file not found: {tsv_file}")
        return []
    except PermissionError:
        print(f"ERROR: Permission denied reading domain file: {tsv_file}")
        return []
    except Exception as e:
        print(f"ERROR: Unexpected error while loading domains from \
            {tsv_file}: {e}")
        return []

    # Longest domain per source
    try:
        longest_per_source = [
            max(domains, key=lambda x: (x[1] - x[0]))
            for domains in source_to_domains.values()
            if domains
        ]
    except Exception as e:
        print(f"ERROR: Failed to determine longest domains: {e}")
        return []

    # Sort longest → shortest
    longest_per_source.sort(key=lambda x: (x[1] - x[0]), reverse=True)

    # Return up to max_domains
    return longest_per_source[:max_domains]


# -------------------
# Tree modifications
# -------------------


def cap_branch_lengths(tree: PhyloTree, max_branch_length: float):
    """
    Prune branch lengths greater than max_branch_length.

    This function:
    - Stores the original distance in node.orig_dist
    - Caps node.dist to max_branch_length
    - Appends '*' to the node.name to mark it

    Parameters
    ----------
    tree : PhyloTree
        The tree to modify (modified in place).
    max_branch_length : float
        Maximum allowed branch length. Must be non-negative.

    Raises
    ------
    ValueError
        If max_branch_length is negative.
    """
    if max_branch_length < 0:
        raise ValueError("max_branch_length must be non-negative")

    for node in tree.traverse("postorder"):
        # Skip root and nodes without numeric dist
        try:
            dist = float(node.dist)
        except Exception:
            continue

        if dist > max_branch_length and not node.is_root():
            node.add_feature("orig_dist", node.dist)
            node.dist = max_branch_length
            node.name = f"{node.name}*"


# -------------------
# Tree layout helpers
# -------------------


def layout(node):
    """Custom layout for leaf nodes with optional sequences."""
    if not node.is_leaf():
        return

    simple_name = parse_leaf_name(node.name)
    if hasattr(node, "orig_dist"):
        simple_name = f"{simple_name} (orig={node.orig_dist:.2f})"

    # Leaf label
    faces.add_face_to_node(
        faces.TextFace("   "), node, column=0, position="branch-right"
    )
    label_face = faces.TextFace(simple_name)
    label_face.fsize = 10
    label_face.fgcolor = "black"
    faces.add_face_to_node(label_face, node, column=1, position="branch-right")
    faces.add_face_to_node(faces.TextFace("      "),
                           node, column=2, position="aligned")

    # Sequence visualization (if linked)
    if hasattr(node, "sequence") and node.sequence:
        max_display = 5000
        for i, aa in enumerate(node.sequence[:max_display]):
            color = aa_to_color(aa)
            if aa == "-":
                rect = faces.RectFace(
                    width=2, height=1, fgcolor=color, bgcolor=color)
            else:
                rect = faces.RectFace(width=1, height=10,
                                      fgcolor=color, bgcolor=color)
            faces.add_face_to_node(
                rect, node, column=3 + i, position="aligned")


def add_header(
    ts: TreeStyle,
    domains: Optional[List[Tuple[int, int, str, str]]],
):
    """
    Add a colored domain header to the tree.

    Parameters
    ----------
    ts : TreeStyle
        TreeStyle instance whose aligned_header will be modified.
    domains : list or None
        List of (start, end, name, source) tuples, or None.
    """
    for col in range(3):
        ts.aligned_header.add_face(faces.TextFace(" "), column=col)

    if not domains:
        return

    for start, end, name, _ in domains:
        color = domain_to_color(name)
        for i in range(start, end + 1):
            rect = faces.RectFace(width=1, height=10,
                                  fgcolor=color, bgcolor=color)
            ts.aligned_header.add_face(rect, column=3 + i)


def add_legend(ts: TreeStyle, domains: Optional[List[Tuple[int, int, str, str]]]):
    """Add a legend for domain colors."""
    seen = {}

    if not domains:
        return

    for _, _, name, _ in domains:
        if name not in seen:
            seen[name] = domain_to_color(name)

    for domain, color in seen.items():
        ts.legend.add_face(
            faces.RectFace(width=10, height=10, fgcolor=color, bgcolor=color), column=0
        )
        ts.legend.add_face(faces.TextFace(f" {domain}", fsize=10), column=1)


# -------------------
# Tree rendering
# -------------------


def render_tree(
    tree: PhyloTree,
    output: str,
    width: int,
    show_branch_lengths: bool,
    horizontal_scale: float,
    domains: Optional[List[Tuple[int, int, str, str]]],
):
    """
    Render the tree to an image file.

    Parameters
    ----------
    tree : PhyloTree
        Tree to render.
    output : str
        Output filename.
    width : int
        Image width in mm.
    show_branch_lengths : bool
        Whether to show branch lengths.
    horizontal_scale : float
        Optional horizontal scale for the tree (passed to TreeStyle.scale).
    domains : list or None
        Domain header data to add (optional).
    """
    ts = TreeStyle()
    ts.show_branch_length = show_branch_lengths
    ts.show_leaf_name = False
    ts.layout_fn = layout
    ts.branch_vertical_margin = 20
    if horizontal_scale:
        ts.scale = horizontal_scale

    if domains:
        add_header(ts, domains)
        add_legend(ts, domains)

    try:
        tree.render(output, w=width, tree_style=ts, units="mm")
        print(f"Tree rendered to {output}")
    except FileNotFoundError:
        print(
            f"ERROR: Cannot write output file (not found or path invalid): {output}")
    except PermissionError:
        print(f"ERROR: Permission denied writing output file: {output}")
    except Exception as e:
        print(f"ERROR: Failed to render tree to {output}: {e}")


def compute_auto_branch_cap(tree: PhyloTree, min_fraction_of_max: float = 0.2) -> Optional[float]:
    """
    Compute an outlier-robust automatic maximum branch length.

    Uses Tukey fence (Q3 + 1.5*IQR), 95th percentile, and median+3*MAD,
    picks a conservative cap and enforces a lower bound based on a fraction
    of max branch length.

    Parameters
    ----------
    tree : PhyloTree
        Tree to analyze (branch lengths taken from node.dist).
    min_fraction_of_max : float
        Minimum fraction of the maximum observed branch length to allow the cap to be.

    Returns
    -------
    float or None
        Calculated cap, or None if no branch lengths are present.
    """
    import math
    from statistics import median

    lengths = []
    for n in tree.traverse("postorder"):
        if n.is_root():
            continue
        try:
            lengths.append(float(n.dist))
        except Exception:
            # skip nodes without numeric distance
            continue

    if not lengths:
        return None

    lengths_sorted = sorted(lengths)

    def quantile(sorted_list, q):
        pos = q * (len(sorted_list) - 1)
        lo = int(math.floor(pos))
        hi = int(math.ceil(pos))
        if lo == hi:
            return sorted_list[lo]
        w = pos - lo
        return sorted_list[lo] * (1 - w) + sorted_list[hi] * w

    try:
        q1 = quantile(lengths_sorted, 0.25)
        q3 = quantile(lengths_sorted, 0.75)
        iqr = q3 - q1
        tukey_upper = q3 + 1.5 * iqr
        p95 = quantile(lengths_sorted, 0.95)
        med = median(lengths_sorted)
        mad = median([abs(x - med) for x in lengths_sorted])
        mad_upper = med + 3 * mad
        candidates = [tukey_upper, p95, mad_upper]

        cap = max(med, min(candidates))
        max_len = max(lengths_sorted)
        cap = max(cap, min_fraction_of_max * max_len)
        return cap
    except Exception as e:
        print(f"WARNING: Failed to compute automatic branch cap: {e}")
        return None


# -------------------
# Main
# -------------------


def main():
    args = parse_arguments()

    # Load tree
    try:
        tree = PhyloTree(args.tree_file, format=1, quoted_node_names=True)
    except FileNotFoundError:
        print(f"ERROR: Tree file not found: {args.tree_file}")
        return
    except PermissionError:
        print(f"ERROR: Permission denied reading tree file: {args.tree_file}")
        return
    except Exception as e:
        print(f"ERROR: Could not load tree from {args.tree_file}: {e}")
        return

    # Link alignment if provided
    if args.alignment:
        try:
            tree.link_to_alignment(args.alignment, alg_format="fasta")
        except FileNotFoundError:
            print(f"WARNING: Alignment file not found: {args.alignment}")
        except Exception as e:
            print(f"WARNING: Could not link alignment {args.alignment}: {e}")

    # Reroot tree if requested
    if args.reroot:
        try:
            target_node = next(
                (
                    leaf
                    for leaf in tree.iter_leaves()
                    if parse_leaf_name(leaf.name) == args.reroot
                ),
                None,
            )
            if target_node:
                tree.set_outgroup(target_node)
            else:
                print(f"WARNING: Leaf '\
                {args.reroot}' not found for rerooting.")
        except Exception as e:
            print(f"WARNING: Could not reroot at '{args.reroot}': {e}")

    # Cap branch lengths (either user-specified or auto)
    if args.max_branch_length is not None:
        cap_value = args.max_branch_length
        print(f"Using user-defined max-branch-length = {cap_value}")
    else:
        cap_value = compute_auto_branch_cap(tree)
        if cap_value is not None:
            print(f"Auto-computed max-branch-length = {cap_value:.5f}")
        else:
            print(
                "WARNING: Could not compute branch-length cap (no numeric branch lengths found).")
            cap_value = None

    if cap_value is not None:
        try:
            before = [float(n.dist)
                      for n in tree.traverse() if not n.is_root()]
        except Exception:
            before = []
        try:
            cap_branch_lengths(tree, cap_value)
        except ValueError as e:
            print(f"ERROR: Invalid cap value: {e}")
        except Exception as e:
            print(f"ERROR: Failed to apply branch cap: {e}")
        finally:
            try:
                after = [float(n.dist)
                         for n in tree.traverse() if not n.is_root()]
            except Exception:
                after = []
        try:
            num_capped = sum(1 for b in before if b >
                             cap_value) if before else 0
            print(f"Capped {num_capped} branches above threshold \
            {cap_value:.5f}.")
        except Exception:
            pass

    # Load domains
    domains = None
    if args.domains:
        domains = load_domains(args.domains, max_domains=args.max_domains)
        if not domains:
            print(f"NOTE: No valid domains loaded from \
                {args.domains} (file empty or malformed).")

    # Render tree
    try:
        render_tree(
            tree,
            output=args.output,
            width=args.width,
            show_branch_lengths=not args.no_branch_lengths,
            horizontal_scale=args.horizontal_scale,
            domains=domains,
        )
    except Exception as e:
        print(f"ERROR: Rendering failed: {e}")


if __name__ == "__main__":
    main()
