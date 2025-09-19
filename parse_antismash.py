#!/usr/bin/env python3
"""
parse_clusterblast.py

Parses ClusterBlast text outputs inside sample directories and produces:
 - per-txt extracted results (long table of blast hits + BGC metadata + merged query-cluster annotations)
 - per-sample merged TSV combining all txt files for that sample

Usage:
  python parse_clusterblast.py --input-dir ./samples --pattern "barcode*" --out-dir ./results --subdir knownclusterblastdirectory --verbose
"""
import argparse
import glob
import os
import re
import pandas as pd
from pathlib import Path

def normalize_gene_id(s: str) -> str:
    """Normalize gene id strings for consistent merging."""
    if pd.isna(s):
        return s
    s = str(s)
    # strip whitespace, invisible chars, trailing punctuation
    s = s.strip()
    s = re.sub(r'[\r\n\t]+', '', s)
    s = s.replace('"', '').replace("'", "")
    # sometimes extra punctuation (like trailing commas/semicolons) appears
    s = s.rstrip(' ,;:')
    return s

def parse_query_cluster_block(lines):
    """
    Given lines (list of strings) that follow the "Table of genes, locations, strands and annotations of query cluster:" header,
    parse until blank line (or until 'Significant hits:' etc).
    Returns a DataFrame with columns: query_gene, start, stop, strand, annotation
    """
    rows = []
    for ln in lines:
        ln = ln.rstrip("\n")
        if ln.strip() == "" or ln.strip().startswith("Significant hits:") or ln.strip().startswith("Details:"):
            break
        # Expect lines like: BacillusspCRS201_03720\t3658517\t3660080\t+\tRibonuclease Y\t
        parts = ln.split('\t')
        if len(parts) >= 5:
            gene = normalize_gene_id(parts[0])
            start = parts[1].strip()
            stop = parts[2].strip()
            strand = parts[3].strip()
            ann = parts[4].strip() if parts[4] != '' else None
            rows.append({
                "query_gene": gene,
                "start": start,
                "stop": stop,
                "strand": strand,
                "annotation": ann
            })
        else:
            # fallback: try whitespace split
            parts_ws = ln.split()
            if len(parts_ws) >= 5:
                gene = normalize_gene_id(parts_ws[0])
                start = parts_ws[1]
                stop = parts_ws[2]
                strand = parts_ws[3]
                ann = " ".join(parts_ws[4:])
                rows.append({
                    "query_gene": gene,
                    "start": start,
                    "stop": stop,
                    "strand": strand,
                    "annotation": ann
                })
            else:
                # skip odd lines
                continue
    df = pd.DataFrame(rows)
    if not df.empty:
        df["query_gene"] = df["query_gene"].astype(str).apply(normalize_gene_id)
    return df

def parse_file(filepath):
    """
    Parse a single ClusterBlast .txt file and return:
      - query_df: DataFrame of query cluster genes (query_gene, start, stop, strand, annotation)
      - hits_df: long DataFrame of all blast hits with appended BGC metadata columns:
        ['BCG_identifier','BCG_name','BCG_type','proteins_with_blast_hits','cum_BLAST_Score',
         'query_gene','subject_gene','pct_identity','blast_score','pct_coverage','evalue','region']
    """
    text = Path(filepath).read_text(encoding='utf-8', errors='ignore')
    lines = text.splitlines()

    # 1) find query cluster block
    query_df = pd.DataFrame()
    for i, ln in enumerate(lines):
        if ln.strip().startswith("Table of genes, locations, strands and annotations of query cluster:"):
            # collect following lines
            tail = lines[i+1:]
            query_df = parse_query_cluster_block(tail)
            break

    # 2) find all '>>' cluster blocks and parse metadata + table of blast hits
    hits_rows = []
    # split on lines that start with ">>" OR on the numeric ordered ">>\n1. BGC..."
    # But in the file pattern used, each cluster block begins with ">>" line or with ">>" later (we'll search for ">>\nN. BGC")
    # We'll iterate searching for lines starting with ">>"
    idxs = [i for i, ln in enumerate(lines) if ln.strip().startswith(">>")]
    # If ">>" markers used differently, also find lines that match r"^\d+\.\s+BGC"
    if not idxs:
        idxs = [i for i, ln in enumerate(lines) if re.match(r'^\d+\.\s+BGC', ln.strip())]
    # For each start index, find the next '>>' or end of file
    for j, start_idx in enumerate(idxs):
        # determine end
        end_idx = len(lines)
        if j+1 < len(idxs):
            end_idx = idxs[j+1]
        block = lines[start_idx:end_idx]
        # Parse BGC id (first non-empty line in block often contains "1. BGC0001089.5")
        bgc_id = None
        source_name = None
        bgc_type = None
        proteins_with_hits = None
        cum_score = None
        # parse header lines inside block
        for ln in block[:15]:  # metadata usually near top of block
            ln_stripped = ln.strip()
            m_bgc = re.match(r'^\d+\.\s+([A-Z]{3}\d+[\d\.\-A-Za-z]*)', ln_stripped)
            if m_bgc:
                bgc_id = m_bgc.group(1)
            if ln_stripped.startswith("Source:"):
                # Source: bacillaene
                source_name = ln_stripped.split("Source:",1)[1].strip()
            if ln_stripped.startswith("Type:"):
                bgc_type = ln_stripped.split("Type:",1)[1].strip()
            if ln_stripped.startswith("Number of proteins with BLAST hits to this cluster:"):
                proteins_with_hits = ln_stripped.split(":",1)[1].strip()
            if ln_stripped.startswith("Cumulative BLAST score:"):
                cum_score = ln_stripped.split(":",1)[1].strip()
            # also support lines like: "1. BGC0001089.5" alone
            m2 = re.match(r'^\d+\.\s+(BGC[0-9\.\-]+)', ln_stripped)
            if m2 and not bgc_id:
                bgc_id = m2.group(1)

        # clean the BGC identifier (user asked to keep the BGC bit — we'll remove trailing dot-suffix if present)
        if bgc_id:
            bgc_id_clean = re.sub(r'\..*$', '', bgc_id)  # strip ".5" etc -> BGC0001089
        else:
            bgc_id_clean = None

        # find the "Table of Blast hits" header inside this block
        hits_start = None
        for k, ln in enumerate(block):
            if ln.strip().startswith("Table of Blast hits"):
                hits_start = k + 1
                break
        if hits_start is None:
            # no hits table here
            continue

        # collect hits lines until blank line or next header
        for ln in block[hits_start:]:
            if ln.strip() == "" or ln.strip().startswith(">>") or re.match(r'^\d+\.\s+BGC', ln.strip()):
                break
            # each hits line expected as tab-separated fields:
            parts = ln.split('\t')
            if len(parts) >= 6:
                q = normalize_gene_id(parts[0])
                subj = parts[1].strip()
                pct_id = parts[2].strip()
                bscore = parts[3].strip()
                pct_cov = parts[4].strip()
                ev = parts[5].strip()
                hits_rows.append({
                    "BCG_identifier": bgc_id_clean,
                    "BCG_name": source_name,
                    "BCG_type": bgc_type,
                    "proteins_with_blast_hits": proteins_with_hits,
                    "cum_BLAST_Score": cum_score,
                    "query_gene": q,
                    "subject_gene": subj,
                    "pct_identity": pct_id,
                    "blast_score": bscore,
                    "pct_coverage": pct_cov,
                    "evalue": ev
                })
            else:
                # try whitespace-split fallback
                parts_ws = ln.split()
                if len(parts_ws) >= 6:
                    q = normalize_gene_id(parts_ws[0])
                    subj = parts_ws[1]
                    pct_id = parts_ws[2]
                    bscore = parts_ws[3]
                    pct_cov = parts_ws[4]
                    ev = parts_ws[5]
                    hits_rows.append({
                        "BCG_identifier": bgc_id_clean,
                        "BCG_name": source_name,
                        "BCG_type": bgc_type,
                        "proteins_with_blast_hits": proteins_with_hits,
                        "cum_BLAST_Score": cum_score,
                        "query_gene": q,
                        "subject_gene": subj,
                        "pct_identity": pct_id,
                        "blast_score": bscore,
                        "pct_coverage": pct_cov,
                        "evalue": ev
                    })
                else:
                    continue

    hits_df = pd.DataFrame(hits_rows)
    # normalize gene ids
    if not hits_df.empty:
        hits_df["query_gene"] = hits_df["query_gene"].astype(str).apply(normalize_gene_id)
    if not query_df.empty:
        query_df["query_gene"] = query_df["query_gene"].astype(str).apply(normalize_gene_id)

    return query_df, hits_df

def process_sample(sample_dir, subdir_name, out_dir, verbose=True):
    sample_name = os.path.basename(sample_dir.rstrip("/"))
    known_dir = os.path.join(sample_dir, subdir_name)
    if not os.path.isdir(known_dir):
        # try lowercase variant
        known_dir = os.path.join(sample_dir, subdir_name.lower())
        if not os.path.isdir(known_dir):
            if verbose:
                print(f"[WARN] {sample_name}: subdir '{subdir_name}' not found -> skipping sample")
            return

    txt_files = sorted(glob.glob(os.path.join(known_dir, "*.txt")))
    if not txt_files:
        if verbose:
            print(f"[WARN] {sample_name}: no .txt files found in {known_dir}")
        return

    sample_out_dir = os.path.join(out_dir, sample_name)
    os.makedirs(sample_out_dir, exist_ok=True)

    per_file_paths = []
    merged_rows = []  # accumulate for sample-level merge

    for txt in txt_files:
        fname = os.path.basename(txt)
        region = re.sub(r'.*?(_c\d+|_C\d+|_c\D*\d+).*', r'\1', fname)
        # region may be like '_c10' or similar; make consistent (strip leading underscore)
        region_clean = region.lstrip('_') if region else ""
        if verbose:
            print(f"[INFO] {sample_name}: processing {fname} ...")

        query_df, hits_df = parse_file(txt)

        # Attach region to hits
        if hits_df.empty:
            if verbose:
                print(f"[INFO] {fname}: no blast-hits parsed, skipping writing for this file.")
            continue

        hits_df["region"] = region_clean

        # Merge hits with query cluster annotation on query_gene (robust)
        # Keep only needed columns from query_df
        if not query_df.empty:
            query_sub = query_df.drop_duplicates(subset=["query_gene"])[["query_gene", "start", "stop", "strand", "annotation"]]
        else:
            query_sub = pd.DataFrame(columns=["query_gene", "start", "stop", "strand", "annotation"])

        # ensure normalization
        for col in ["query_gene"]:
            if col in hits_df.columns:
                hits_df[col] = hits_df[col].astype(str).apply(normalize_gene_id)
            if col in query_sub.columns:
                query_sub[col] = query_sub[col].astype(str).apply(normalize_gene_id)

        merged = pd.merge(hits_df, query_sub, on="query_gene", how="left", indicator=True)

        # Diagnostics: which rows failed to merge?
        unmatched = merged[merged["_merge"] == "left_only"]
        if not unmatched.empty:
            unmatched_path = os.path.join(sample_out_dir, fname + ".unmatched_rows.csv")
            unmatched.to_csv(unmatched_path, sep="\t", index=False)
            if verbose:
                print(f"[DIAG] {fname}: {len(unmatched)} hits had no matching query-cluster row. Saved unmatched sample to {unmatched_path}")
        else:
            if verbose:
                print(f"[DIAG] {fname}: all hits matched query-cluster genes")

        # drop merge indicator column
        merged = merged.drop(columns=["_merge"])

        # Save per-file result
        out_file = os.path.join(sample_out_dir, fname + ".results.tsv")
        merged.to_csv(out_file, sep="\t", index=False)
        if verbose:
            print(f"[INFO] {fname}: wrote parsed+merged results to {out_file}")

        per_file_paths.append(out_file)
        merged_rows.append(merged)

    # after all files, create a sample-level merged file
    if merged_rows:
        combined = pd.concat(merged_rows, sort=False).reset_index(drop=True)
        sample_merged_path = os.path.join(out_dir, f"{sample_name}.merged_results.tsv")
        combined.to_csv(sample_merged_path, sep="\t", index=False)
        if verbose:
            print(f"[SUCCESS] {sample_name}: sample merged file written to {sample_merged_path}")
    else:
        if verbose:
            print(f"[WARN] {sample_name}: no parsed data to merge for sample.")

def main():
    p = argparse.ArgumentParser(description="Parse ClusterBlast .txt files across multiple samples and merge results.")
    p.add_argument("--input-dir", "-i", required=True, help="Path to directory containing sample directories")
    p.add_argument("--pattern", "-p", default="barcode*", help="Glob pattern for sample directory names (e.g. barcode* or sample*)")
    p.add_argument("--out-dir", "-o", required=True, help="Output base directory (script will create per-sample subfolders)")
    p.add_argument("--subdir", default="knownclusterblastdirectory", help="Subdirectory inside each sample where txt files are located (default: knownclusterblastdirectory)")
    p.add_argument("--verbose", "-v", action="store_true", help="Print processing checkpoints")
    args = p.parse_args()

    input_dir = args.input_dir
    pat = args.pattern
    out_dir = args.out_dir
    subdir = args.subdir

    os.makedirs(out_dir, exist_ok=True)
    # find sample directories matching pattern
    search_pattern = os.path.join(input_dir, pat)
    sample_dirs = [d for d in glob.glob(search_pattern) if os.path.isdir(d)]
    if args.verbose:
        print(f"[START] Found {len(sample_dirs)} sample dirs matching pattern '{pat}' in {input_dir}")

    for sample in sample_dirs:
        process_sample(sample, subdir, out_dir, verbose=args.verbose)

    if args.verbose:
        print("[DONE] All samples processed.")

if __name__ == "__main__":
    main()
