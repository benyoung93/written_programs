#!/usr/bin/env python3
"""
parse_busco.py: Parse BUSCO short_summary.txt files into a long-format TSV table.

Usage:
    ./parse_busco.py -i /path/to/busco_results -o busco_summary.tsv

Columns in output:
    Sample   Database   Version   Mode   Measure   Value
"""

import argparse
from pathlib import Path
import re
import sys

# Regex patterns for BUSCO summary lines
PATTERNS = {
    "Complete": r"^\s*(\d+)\s+Complete BUSCOs",
    "Complete_Single_Copy": r"^\s*(\d+)\s+Complete and single-copy BUSCOs",
    "Complete_Duplicated": r"^\s*(\d+)\s+Complete and duplicated BUSCOs",
    "Fragmented": r"^\s*(\d+)\s+Fragmented BUSCOs",
    "Missing": r"^\s*(\d+)\s+Missing BUSCOs",
    "Total": r"^\s*(\d+)\s+Total BUSCO groups searched",
}


def parse_busco_summary(file_path, sample):
    """
    Parse a BUSCO summary file.

    Returns:
        List of tuples: (sample, database, version, mode, measure, value)
    """
    rows = []
    mode = "unknown"
    db = "unknown"
    version = "unknown"

    with open(file_path) as f:
        for line in f:
            # BUSCO version
            if line.startswith("# BUSCO version is:"):
                version_num = line.strip().split(":")[1].strip()
                version = f"v{version_num}"

            # Lineage dataset / database
            if line.startswith("# The lineage dataset is:"):
                db = line.strip().split(":")[1].split()[0].strip()

            # Run mode
            if "BUSCO was run in mode:" in line:
                if "proteins" in line:
                    mode = "protein"
                elif "euk_genome_met" in line:
                    mode = "nucleotide"

            # Metrics
            for measure, pattern in PATTERNS.items():
                match = re.match(pattern, line)
                if match:
                    value = match.group(1)
                    rows.append((sample, db, version, mode, measure, value))
    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Parse BUSCO short_summary.txt files into a long-format TSV table."
    )
    parser.add_argument(
        "-i", "--indir",
        required=True,
        help="Base input directory containing BUSCO result folders (e.g., seedXXX_asco)."
    )
    parser.add_argument(
        "-o", "--outfile",
        default="-",
        help="Output TSV file. Default: stdout"
    )
    args = parser.parse_args()

    indir = Path(args.indir)
    out = sys.stdout if args.outfile == "-" else open(args.outfile, "w")

    # Write header
    print("Sample\tDatabase\tVersion\tMode\tMeasure\tValue", file=out)

    # Loop through subdirectories
    for sample_dir in sorted(indir.iterdir()):
        if not sample_dir.is_dir():
            continue

        sample = sample_dir.name

        # Find all short_summary.txt files in this sample directory
        for summary_file in sample_dir.rglob("short_summary.*.txt"):
            rows = parse_busco_summary(summary_file, sample)
            for row in rows:
                print("\t".join(row), file=out)

    if out is not sys.stdout:
        out.close()


if __name__ == "__main__":
    main()
