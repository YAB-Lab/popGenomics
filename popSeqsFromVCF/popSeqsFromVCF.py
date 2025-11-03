#!/usr/bin/env python3
"""
Extract a user-specified feature from a GTF, filter an all-sites VCF to those
intervals, run vcf2phylip to build an alignment, and (if minus strand) reverse
the final PHYLIP alignment (optionally reverse-complement).

Requirements:
  - Python 3.8+
  - bcftools in PATH
  - tabix in PATH (bcftools will use it to index the VCF)
  - vcf2phylip in PATH (or provide full path via --vcf2phylip)

Notes:
  - GTF is assumed to be standard 1-based inclusive coordinates.
  - BED output is 0-based half-open, merged and sorted.
  - Reversal currently assumes SEQUENTIAL PHYLIP output from vcf2phylip.
    (We can extend to interleaved in a later iteration if needed.)
"""

import argparse
import os
import sys
import gzip
import shutil
import subprocess
import time
from pathlib import Path
from collections import defaultdict
from typing import List, Tuple, Dict, Optional

# ---------- GTF parsing & interval handling ----------

def parse_gtf_for_feature(
    gtf_path: str,
    feature_type: str,
    feature_id: str
) -> Tuple[List[Tuple[str, int, int]], str]:
    """
    Return intervals (chrom, start0, end) for the requested (feature_type, feature_id),
    merging across records that match either gene_id==feature_id OR transcript_id==feature_id.
    Also return the strand ('+' or '-') inferred from the first matching record.

    GTF columns: chrom, source, feature, start, end, score, strand, frame, attributes
    Coordinates in GTF are 1-based inclusive; convert to 0-based half-open for BED.
    """
    intervals: List[Tuple[str, int, int]] = []
    strand: Optional[str] = None

    def attrs_to_dict(attr_field: str) -> Dict[str, str]:
        d = {}
        # attributes like: key "value"; key2 "value2";
        for part in attr_field.strip().strip(';').split(';'):
            part = part.strip()
            if not part:
                continue
            # allow both key "value" and key=value styles
            if ' ' in part:
                k, v = part.split(' ', 1)
                d[k] = v.strip().strip('"')
            elif '=' in part:
                k, v = part.split('=', 1)
                d[k] = v.strip().strip('"')
        return d

    opener = gzip.open if gtf_path.endswith(('.gz', '.bgz')) else open
    with opener(gtf_path, 'rt') as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            chrom, _source, feat, start, end, _score, sstrand, _frame, attrs = fields
            if feat != feature_type:
                continue
            ad = attrs_to_dict(attrs)
            gid = ad.get('gene_id')
            tid = ad.get('transcript_id')
            if feature_id != gid and feature_id != tid:
                continue
            try:
                start0 = int(start) - 1
                end1 = int(end)  # half-open in BED
            except ValueError:
                continue
            if strand is None:
                strand = sstrand
            intervals.append((chrom, start0, end1))

    if not intervals:
        raise RuntimeError(
            f"No intervals found for feature_type '{feature_type}' with "
            f"gene_id/transcript_id '{feature_id}' in {gtf_path}"
        )

    if strand is None:
        strand = '+'

    # merge per-chrom
    intervals = merge_intervals(intervals)
    return intervals, strand


def merge_intervals(intervals: List[Tuple[str, int, int]]) -> List[Tuple[str, int, int]]:
    """Merge overlapping/adjacent intervals within each chromosome."""
    by_chr: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for chrom, s, e in intervals:
        by_chr[chrom].append((s, e))

    merged: List[Tuple[str, int, int]] = []
    for chrom, spans in by_chr.items():
        spans.sort()
        cur_s, cur_e = spans[0]
        for s, e in spans[1:]:
            if s <= cur_e:                # overlap or adjacency
                if e > cur_e:
                    cur_e = e
            else:
                merged.append((chrom, cur_s, cur_e))
                cur_s, cur_e = s, e
        merged.append((chrom, cur_s, cur_e))

    # sort lexicographically: chrom, start
    merged.sort(key=lambda x: (x[0], x[1]))
    return merged


def write_bed(intervals: List[Tuple[str, int, int]], out_bed: str) -> None:
    with open(out_bed, 'w') as out:
        for chrom, s, e in intervals:
            out.write(f"{chrom}\t{s}\t{e}\n")


# ---------- External tools ----------

def run_cmd(cmd: List[str], check: bool = True) -> None:
    try:
        subprocess.run(cmd, check=check)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Command failed ({e.returncode}): {' '.join(cmd)}") from e


def ensure_executable_in_path(exe: str, friendly: str) -> None:
    if shutil.which(exe) is None:
        raise RuntimeError(
            f"Could not find '{exe}' in PATH. Please install or provide its full path. "
            f"(Needed for {friendly})"
        )


# ---------- PHYLIP handling ----------

def reverse_sequence(seq: str) -> str:
    return seq[::-1]


def revcomp_sequence(seq: str) -> str:
    # Full IUPAC DNA complement (upper & lower)
    comp_map = str.maketrans(
        "ACGTRYMKSWBDHVNacgtrymkswbdhvn-?.",
        "TGCAYRKMWSVHDBNtgcayrkmwsvhdbn-?."
    )
    return seq.translate(comp_map)[::-1]


def read_sequential_phylip(path: str) -> Dict[str, str]:
    """
    Very simple reader for *sequential* PHYLIP format.
    Assumes each taxon has one line with its full sequence after the header.
    """
    with open(path, 'r') as fh:
        header = fh.readline()
        if not header:
            raise RuntimeError("Empty PHYLIP file?")
        parts = header.strip().split()
        if len(parts) < 2:
            raise RuntimeError(f"PHYLIP header malformed: {header.strip()}")
        try:
            nseq = int(parts[0])
            _length = int(parts[1])
        except ValueError:
            raise RuntimeError(f"PHYLIP header malformed: {header.strip()}")

        data: Dict[str, str] = {}
        for _ in range(nseq):
            line = fh.readline()
            if not line:
                raise RuntimeError("Unexpected end of PHYLIP while reading sequences.")
            # PHYLIP "relaxed" often uses name then spaces then sequence
            name_seq = line.rstrip('\n')
            if not name_seq.strip():
                # tolerate blank lines
                continue
            # split on whitespace once
            parts = name_seq.split(None, 1)
            if len(parts) == 1:
                raise RuntimeError(
                    "Could not parse sequence line (expected 'name sequence'): "
                    f"'{name_seq}'"
                )
            name, seq = parts[0], parts[1].replace(" ", "")
            data[name] = seq
        return data


def write_sequential_phylip(path: str, data: Dict[str, str]) -> None:
    names = list(data.keys())
    if not names:
        raise RuntimeError("No sequences to write.")
    Ls = {len(s) for s in data.values()}
    if len(Ls) != 1:
        raise RuntimeError("Sequences are not aligned to equal length.")
    L = Ls.pop()
    with open(path, 'w') as out:
        out.write(f"{len(names)} {L}\n")
        # pad names to a reasonable width (max name length or 10)
        name_width = max(10, max(len(n) for n in names))
        for n in names:
            out.write(f"{n.ljust(name_width)} {data[n]}\n")


def discover_phylip_output(candidates_dirs, t0, preferred_stem=None):
    """
    Search candidate dirs for PHYLIP files created/modified after t0.
    Prefer files whose names include preferred_stem if provided.
    Returns a single resolved path or None.
    """
    exts = {".phy", ".phylip"}
    found = []
    for d in candidates_dirs:
        d = Path(d)
        if not d.exists():
            continue
        for p in d.glob("*"):
            if p.suffix.lower() in exts and p.is_file():
                try:
                    mtime = p.stat().st_mtime
                except OSError:
                    continue
                if mtime >= t0 - 1:  # small slack
                    found.append(p.resolve())

    if not found:
        return None

    if preferred_stem:
        pf = [p for p in found if preferred_stem in p.name]
        if len(pf) == 1:
            return pf[0]
        if len(pf) > 1:
            return max(pf, key=lambda p: p.stat().st_mtime)

    return max(found, key=lambda p: p.stat().st_mtime)

# ---------- FASTA/PHYLIP handling ----------#

from pathlib import Path
import time

def read_fasta(path: str) -> dict:
    seqs = {}
    name = None
    parts = []
    with open(path, "r") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(parts).replace(" ", "").replace("\t", "")
                name = line[1:].strip().split()[0]
                parts = []
            else:
                parts.append(line.strip())
        if name is not None:
            seqs[name] = "".join(parts).replace(" ", "").replace("\t", "")
    return seqs

def write_fasta(path: str, data: dict, wrap: int = 0):
    with open(path, "w") as out:
        for name, seq in data.items():
            out.write(f">{name}\n")
            if wrap and wrap > 0:
                for i in range(0, len(seq), wrap):
                    out.write(seq[i:i+wrap] + "\n")
            else:
                out.write(seq + "\n")

def fasta_to_phylip(data: dict) -> dict:
    # Already “aligned” if coming from vcf2phylip; ensure equal length.
    Ls = {len(s) for s in data.values()}
    if len(Ls) != 1:
        raise RuntimeError("FASTA sequences are not equal length; cannot write PHYLIP.")
    return data

def phylip_to_fasta(data: dict) -> dict:
    return data

def revcomp_sequence(seq: str) -> str:
    comp_map = str.maketrans(
        "ACGTRYMKSWBDHVNacgtrymkswbdhvn-?.",
        "TGCAYRKMWSVHDBNtgcayrkmwsvhdbn-?."
    )
    return seq.translate(comp_map)[::-1]

def discover_outputs(candidates_dirs, t0, preferred_stem=None):
    """
    Return dict of discovered newest files by ext: {'.phy': Path|None, '.phylip': Path|None, '.fa': Path|None, '.fasta': Path|None}
    Only files with mtime >= t0-1 are considered. Preference to names containing preferred_stem.
    """
    exts = {".phy", ".phylip", ".fa", ".fasta"}
    found = {e: [] for e in exts}
    for d in candidates_dirs:
        d = Path(d)
        if not d.exists():
            continue
        for p in d.glob("*"):
            if p.suffix.lower() in exts and p.is_file():
                try:
                    mtime = p.stat().st_mtime
                except OSError:
                    continue
                if mtime >= t0 - 1:
                    found[p.suffix.lower()].append(p.resolve())

    chosen = {}
    for e, items in found.items():
        if not items:
            chosen[e] = None
            continue
        if preferred_stem:
            pf = [p for p in items if preferred_stem in p.name]
            if pf:
                chosen[e] = max(pf, key=lambda p: p.stat().st_mtime)
                continue
        chosen[e] = max(items, key=lambda p: p.stat().st_mtime)
    return chosen


# ---------- Main pipeline ----------

def main():
    ap = argparse.ArgumentParser(
        description="Extract feature intervals from GTF, filter VCF with bcftools, "
                    "run vcf2phylip, and reverse final PHYLIP if minus strand."
    )
    ap.add_argument("--vcf", required=True, help="All-sites VCF (bgzipped or plain).")
    ap.add_argument("--gtf", required=True, help="Reference GTF with annotations.")
    ap.add_argument("--feature-id", required=True,
                    help="ID to match against gene_id or transcript_id.")
    ap.add_argument("--feature-type", required=True,
                    help="GTF feature type to extract (e.g., exon, CDS, UTR, transcript).")
    ap.add_argument("--out-prefix", required=True,
                    help="Prefix for outputs (BED, filtered VCF, PHYLIP).")
    ap.add_argument("--vcf2phylip", default="vcf2phylip",
                    help="Path/command for vcf2phylip. Default: vcf2phylip in PATH.")
    ap.add_argument("--vcf2phylip-args", default="", help="Extra args to pass through.")
    ap.add_argument("--keep-temp", action="store_true",
                    help="Keep intermediate BED and filtered VCF.")
    ap.add_argument(
        "--out-formats",
        choices=["phylip", "fasta", "both"],
        default="phylip",
        help="Which final outputs to produce (after minus-strand RC if needed)."
    )
    ap.add_argument("--phylip-out", default=None,
                    help="Path for final PHYLIP (default: {out-prefix}.phy)")
    ap.add_argument("--fasta-out", default=None,
                    help="Path for final FASTA (default: {out-prefix}.fa)")
    args = ap.parse_args()

    # Resolve outputs
    bed_path = f"{args.out_prefix}.feature.bed"
    filtered_vcf = f"{args.out_prefix}.filtered.vcf.gz"
    phylip_out = args.phylip_out or f"{args.out_prefix}.phy"
    fasta_out = args.fasta_out or f"{args.out_prefix}.fa"

    # Dependencies
    ensure_executable_in_path("bcftools", "VCF filtering")
    ensure_executable_in_path("tabix", "VCF indexing")
    ensure_executable_in_path(args.vcf2phylip.split()[0], "PHYLIP conversion")

    # 1) Pull intervals & strand from GTF
    print(f"[INFO] Parsing GTF: {args.gtf}", file=sys.stderr)
    intervals, strand = parse_gtf_for_feature(args.gtf, args.feature_type, args.feature_id)
    print(f"[INFO] Found {len(intervals)} merged intervals; strand = {strand}", file=sys.stderr)
    print(f"[INFO] Feature: id={args.feature_id} type={args.feature_type} strand={strand}", file=sys.stderr)

    # 2) Write BED
    print(f"[INFO] Writing BED: {bed_path}", file=sys.stderr)
    write_bed(intervals, bed_path)

    # 3) Filter VCF to BED regions
    print(f"[INFO] Filtering VCF with bcftools -> {filtered_vcf}", file=sys.stderr)
    run_cmd(["bcftools", "view", "-R", bed_path, "-Oz", "-o", filtered_vcf, args.vcf])
    run_cmd(["tabix", "-f", "-p", "vcf", filtered_vcf])

    # 4) Run vcf2phylip (your desired CLI: vcf2phylip -i <vcf> -f ...)
    print(f"[INFO] Running vcf2phylip (-f) on filtered VCF", file=sys.stderr)
    extra = args.vcf2phylip_args.strip().split() if args.vcf2phylip_args.strip() else []

    t0 = time.time()
    cwd = os.getcwd()
    out_dir = os.path.dirname(os.path.abspath(args.out_prefix)) or cwd
    preferred_stem = os.path.basename(args.out_prefix)

    cmd = [args.vcf2phylip, "-i", filtered_vcf, "-f", "-r", "-w"] + extra
    run_cmd(cmd)

    candidates_dirs = {cwd, out_dir}
    for ix, tok in enumerate(extra):
        if tok.lower() in {"--output-folder", "--output_dir", "--outdir", "-d"} and ix + 1 < len(extra):
            candidates_dirs.add(os.path.abspath(extra[ix + 1]))

    disc = discover_outputs(candidates_dirs=list(candidates_dirs), t0=t0, preferred_stem=preferred_stem)
    phy_src = disc.get(".phy") or disc.get(".phylip")
    fa_src  = disc.get(".fa")  or disc.get(".fasta")

    if args.out_formats in {"phylip", "both"} and phy_src is None and fa_src is None:
        raise RuntimeError("No PHYLIP/FASTA file was produced by vcf2phylip (checked newest after run).")
    if args.out_formats == "phylip" and (phy_src is None and fa_src is None):
        raise RuntimeError("Requested PHYLIP but no PHYLIP/FASTA discovered to convert from.")
    if args.out_formats == "fasta" and (fa_src is None and phy_src is None):
        raise RuntimeError("Requested FASTA but no PHYLIP/FASTA discovered to convert from.")

    print(f"[INFO] Discovered outputs — PHYLIP: {phy_src if phy_src else 'none'}, FASTA: {fa_src if fa_src else 'none'}", file=sys.stderr)
    print(f"[INFO] Feature: id={args.feature_id} type={args.feature_type} strand={strand}", file=sys.stderr)

    # Load one or both into memory, RC if needed, and write the requested formats.
    phylip_data = None
    fasta_data  = None

    # Prefer loading from native source; otherwise convert.
    if phy_src:
        phylip_data = read_sequential_phylip(str(phy_src))
        fasta_data  = phylip_to_fasta(phylip_data)
    elif fa_src:
        fasta_data  = read_fasta(str(fa_src))
        phylip_data = fasta_to_phylip(fasta_data)

    # Apply minus-strand reverse-complement once on the in-memory sequences
    if strand == '-':
        print("[INFO] Minus strand detected: reverse-complementing final alignments.", file=sys.stderr)
        phylip_data = {k: revcomp_sequence(v) for k, v in phylip_data.items()} if phylip_data else None
        fasta_data  = {k: revcomp_sequence(v) for k, v in fasta_data.items()}  if fasta_data  else None

    # Write requested outputs deterministically to user-specified paths
    if args.out_formats in {"phylip", "both"} and phylip_data:
        write_sequential_phylip(phylip_out, phylip_data)
        print(f"[DONE] PHYLIP written: {phylip_out}", file=sys.stderr)

    if args.out_formats in {"fasta", "both"} and fasta_data:
        write_fasta(fasta_out, fasta_data, wrap=0)
        print(f"[DONE] FASTA written:  {fasta_out}", file=sys.stderr)

    # Clean discovered raw files unless keeping temps and they differ from finals
    if not args.keep_temp:
        for src in [phy_src, fa_src]:
            if src is None:
                continue
            abs_src = os.path.abspath(str(src))
            if (args.out_formats in {"phylip", "both"}) and abs_src == os.path.abspath(phylip_out):
                continue
            if (args.out_formats in {"fasta", "both"}) and abs_src == os.path.abspath(fasta_out):
                continue
            try:
                os.remove(abs_src)
            except OSError:
                pass


    # 5) Minus strand? Reverse-complement columns and write to phylip_out.
    # if strand == '-':
    #     print("[INFO] Minus strand detected: reverse-complementing final PHYLIP alignment.", file=sys.stderr)
    #     data = read_sequential_phylip(str(resolved_src))
    #     data = {k: revcomp_sequence(v) for k, v in data.items()}
    #     write_sequential_phylip(phylip_out, data)
    #     # optionally clean up the discovered source if not keeping temps
    #     if not args.keep_temp and os.path.abspath(str(resolved_src)) != os.path.abspath(phylip_out):
    #         try:
    #             os.remove(str(resolved_src))
    #         except OSError:
    #             pass
    # else:
    #     # Plus strand: move the discovered file to the requested output path (for determinism).
    #     if os.path.abspath(str(resolved_src)) != os.path.abspath(phylip_out):
    #         shutil.move(str(resolved_src), phylip_out)


    # # 5) Minus-strand? Reverse (or reverse-complement) columns and write to phylip_out.
    # if strand == '-':
    #     print("[INFO] Minus strand: reversing final PHYLIP alignment columns.", file=sys.stderr)
    #     data = read_sequential_phylip(str(resolved_src))
    #     if args.reverse_complement:
    #         data = {k: revcomp_sequence(v) for k, v in data.items()}
    #     else:
    #         data = {k: reverse_sequence(v) for k, v in data.items()}
    #     write_sequential_phylip(phylip_out, data)
    #     # optionally clean up the discovered source if not keeping temps
    #     if not args.keep_temp and os.path.abspath(str(resolved_src)) != os.path.abspath(phylip_out):
    #         try:
    #             os.remove(str(resolved_src))
    #         except OSError:
    #             pass
    # else:
    #     # Plus strand: move the discovered file to the requested output path (for determinism).
    #     if os.path.abspath(str(resolved_src)) != os.path.abspath(phylip_out):
    #         shutil.move(str(resolved_src), phylip_out)


    # 6) Cleanup
    if not args.keep_temp:
        try:
            os.remove(bed_path)
        except OSError:
            pass
        try:
            os.remove(filtered_vcf)
        except OSError:
            pass
        try:
            os.remove(filtered_vcf + ".tbi")
        except OSError:
            pass
        # try:
        #     if os.path.exists(phylip_out_tmp):
        #         os.remove(phylip_out_tmp)
        except OSError:
            pass

    print(f"[DONE] PHYLIP written: {phylip_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
