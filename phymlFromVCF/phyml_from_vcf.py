#!/usr/bin/env python3
"""
Build a phylogenetic tree with PhyML from (optionally filtered) VCF.

Pipeline:
  1) (optional) GTF/GFF -> BED (feature-based extraction)
  2) bcftools view ...   -> subset by region/sample/expr or BED
  3) vcf2phylip.py       -> PHYLIP alignment
  4) phyml               -> tree

"""

import argparse
import os
import sys
import shlex
import subprocess
import tempfile
from pathlib import Path
from collections import defaultdict

def run(cmd, cwd=None, log_prefix=""):
    """Run a shell command, stream output, and raise on failure."""
    print(f"{log_prefix}$ {' '.join(shlex.quote(c) for c in cmd)}", flush=True)
    proc = subprocess.Popen(cmd, cwd=cwd)
    proc.wait()
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd)
    return proc.returncode

def parse_args():
    p = argparse.ArgumentParser(
        description="Create a PhyML tree from a VCF with optional region/sample/expr filtering, and optional feature-based region selection from GTF/GFF.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # I/O
    p.add_argument("-i", "--vcf", required=True, help="Input VCF/BCF (bgzipped+indexed recommended).")
    p.add_argument("-o", "--out-prefix", required=True,
                   help="Output prefix (e.g., Results/trees/mytree). Files will be mytree.phy, mytree.phy_phyml_tree.txt, etc.")
    p.add_argument("--workdir", default=None, help="Working directory (defaults to a temporary directory).")
    p.add_argument("--keep-intermediates", action="store_true", help="Keep intermediate files.")

    # Regions (raw)
    p.add_argument("-r", "--region", action="append",
                   help="Region(s) like chr:start-end or contig. Can be given multiple times.")
    p.add_argument("-R", "--regions-bed",
                   help="BED file of regions (passed to bcftools -R).")

    # Regions via GTF/GFF features
    ganno = p.add_argument_group("Feature-based region selection (GTF/GFF -> BED)")
    ganno.add_argument("--anno", help="GTF or GFF file.")
    ganno.add_argument("--feature-type", help="Feature type to extract (e.g., gene, exon, CDS).")
    ganno.add_argument("--feature-attr", default="gene_id",
                       help="Attribute key to match (GTF example: gene_id/transcript_id; GFF examples: ID/Name).")
    ganno.add_argument("--feature-id", action="append",
                       help="One or more attribute values to select (can be given multiple times).")
    ganno.add_argument("--feature-id-file", help="File with one attribute value per line to select.")
    ganno.add_argument("--merge-bed", action="store_true",
                       help="Merge overlapping intervals before using bcftools -R.")
    ganno.add_argument("--pad", type=int, default=0,
                       help="Pad merged intervals by N bp on both sides (>=0). Applied after merging.")

    # Sample selection
    g = p.add_argument_group("Sample (sub)set")
    g.add_argument("--samples", nargs="+", help="List of samples to INCLUDE (passed to bcftools -s).")
    g.add_argument("--samples-file", help="File with one sample per line to INCLUDE (passed to bcftools -S).")
    g.add_argument("--exclude-samples", nargs="+", help="List of samples to EXCLUDE.")
    g.add_argument("--exclude-samples-file", help="File with one sample per line to EXCLUDE.")
    g.add_argument("--force-samples", action="store_true", help="Pass --force-samples to bcftools view.")

    # bcftools filtering
    g2 = p.add_argument_group("bcftools filters")
    g2.add_argument("-in", "--include-expr",
                    help="bcftools -in EXPRESSION (e.g., 'QUAL>30 && DP>10').")
    g2.add_argument("-ex", "--exclude-expr",
                    help="bcftools -ex EXPRESSION.")
    g2.add_argument("--bcftools-args", default="",
                    help="Extra args to append to 'bcftools view' (raw string, split with shlex).")

    # vcf2phylip
    p.add_argument("--vcf2phylip-path", default="vcf2phylip.py",
                   help="Path to your vcf2phylip.py script.")
    p.add_argument("--vcf2phylip-extra", default="",
                   help="Extra args for vcf2phylip.py (raw string, split with shlex).")

    # phyml
    p.add_argument("--phyml-args", default="-d nt -b 0",
                   help="Extra args to pass to phyml (raw string, split with shlex). Example: '-d nt -m GTR -b 100'.")

    # misc
    p.add_argument("--dry-run", action="store_true", help="Print commands and exit.")
    return p.parse_args()

def parse_attr_field(attr_str):
    """
    Parse GTF/GFF attribute column into a dict.
    Handles:
      - GTF style: key "value"; key "value";
      - GFF3 style: key=value;key2=value2;
    """
    d = {}
    if not attr_str or attr_str == ".":
        return d
    # Try GFF3 key=value; first
    semi = [x for x in attr_str.strip().split(";") if x.strip()]
    for item in semi:
        item = item.strip()
        if "=" in item:
            k, v = item.split("=", 1)
            d[k.strip()] = v.strip()
        else:
            # Possibly GTF style: key "value"
            parts = item.split()
            if len(parts) >= 2:
                k = parts[0].strip()
                v = " ".join(parts[1:]).strip()
                v = v.strip('"')
                d[k] = v
    return d

def load_feature_ids(feature_id_list, feature_id_file):
    ids = set()
    if feature_id_list:
        for x in feature_id_list:
            ids.add(x.strip())
    if feature_id_file:
        with open(feature_id_file) as fh:
            for line in fh:
                s = line.strip()
                if s:
                    ids.add(s)
    if not ids:
        sys.exit("ERROR: No feature IDs provided. Use --feature-id and/or --feature-id-file.")
    return ids

def gtf_gff_to_bed(anno_path, feature_type, feature_attr, ids_set, bed_out_path,
                   merge=False, pad=0):
    """
    Create BED from GTF/GFF selecting rows where column 3 == feature_type
    and attributes[feature_attr] in ids_set.

    Writes a BED with columns: chrom, start(0-based), end, name, score '.', strand
    """
    intervals = defaultdict(list)  # key: (chrom, strand) -> list of (start, end, name)
    with open(anno_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, src, ftype, start, end, score, strand, phase, attrs = parts
            if feature_type and ftype != feature_type:
                continue
            ad = parse_attr_field(attrs)
            val = ad.get(feature_attr)
            if val is None or val not in ids_set:
                continue
            try:
                s = int(start)
                e = int(end)
            except ValueError:
                continue
            # Convert to BED 0-based, half-open: GTF/GFF are 1-based inclusive
            bed_s = max(0, s - 1)
            bed_e = e
            name = val
            # store raw; merging later if requested
            intervals[(chrom, strand)].append((bed_s, bed_e, name))

    # Optionally merge per contig/strand
    if merge:
        merged = defaultdict(list)
        for key, lst in intervals.items():
            lst_sorted = sorted(lst, key=lambda x: x[0])
            cur_s, cur_e, cur_name = None, None, None
            for s, e, name in lst_sorted:
                if cur_s is None:
                    cur_s, cur_e, cur_name = s, e, name
                else:
                    if s <= cur_e:  # overlap/adjacent -> merge
                        if e > cur_e:
                            cur_e = e
                        # name could be concatenated; to keep BED clean, keep the first name
                    else:
                        merged[key].append((cur_s, cur_e, cur_name))
                        cur_s, cur_e, cur_name = s, e, name
            if cur_s is not None:
                merged[key].append((cur_s, cur_e, cur_name))
        intervals = merged

    # Optional padding
    if pad and pad > 0:
        padded = defaultdict(list)
        for (chrom, strand), lst in intervals.items():
            for s, e, name in lst:
                ps = max(0, s - pad)
                pe = e + pad
                padded[(chrom, strand)].append((ps, pe, name))
        intervals = padded

    # Write BED
    n_written = 0
    with open(bed_out_path, "w") as out:
        for (chrom, strand), lst in intervals.items():
            for s, e, name in sorted(lst, key=lambda x: (x[0], x[1])):
                out.write(f"{chrom}\t{s}\t{e}\t{name}\t.\t{strand}\n")
                n_written += 1

    if n_written == 0:
        sys.exit("ERROR: No BED intervals were produced from the given feature filters. "
                 "Check --feature-type / --feature-attr / --feature-id inputs.")
    return bed_out_path

def build_bcftools_cmd(args, out_vcf, bed_from_features=None):
    cmd = ["bcftools", "view", "--output-file", out_vcf, "--output-type", "v"]

    # Regions precedence: if feature-BED provided, we must not also use -r/-R.
    if bed_from_features:
        if args.region or args.regions_bed:
            sys.exit("ERROR: Provide EITHER feature-based selection (via --anno/--feature-*) OR -r/-R, not both.")
        cmd += ["-R", bed_from_features]
    else:
        if args.region and args.regions_bed:
            sys.exit("ERROR: Use either --region (-r) OR --regions-bed (-R), not both.")
        if args.region:
            for r in args.region:
                cmd += ["-r", r]
        if args.regions_bed:
            cmd += ["-R", args.regions_bed]

    # Include samples
    if args.samples and args.samples_file:
        sys.exit("ERROR: Use either --samples or --samples-file, not both.")
    if args.samples:
        cmd += ["-s", ",".join(args.samples)]
    if args.samples_file:
        cmd += ["-S", args.samples_file]

    # Exclude samples
    if args.exclude_samples and args.exclude_samples_file:
        sys.exit("ERROR: Use either --exclude-samples or --exclude-samples-file, not both.")
    if args.exclude_samples:
        excl = "^" + ",".join(args.exclude_samples)
        cmd += ["-s", excl]
    if args.exclude_samples_file:
        cmd += ["-S", f"^{args.exclude_samples_file}"]

    # Expressions
    if args.include_expr and args.exclude_expr:
        sys.exit("ERROR: Use either --include-expr (-i) OR --exclude-expr (-e), not both.")
    if args.include_expr:
        cmd += ["-i", args.include_expr]
    if args.exclude_expr:
        cmd += ["-e", args.exclude_expr]

    # Force samples
    if args.force_samples:
        cmd += ["--force-samples"]

    # Extra user args
    if args.bcftools_args.strip():
        cmd += shlex.split(args.bcftools_args)

    # Input VCF last
    cmd.append(args.vcf)
    return cmd

def main():
    args = parse_args()
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    # workspace
    tempdir_obj = None
    if args.workdir:
        workdir = Path(args.workdir)
        workdir.mkdir(parents=True, exist_ok=True)
    else:
        tempdir_obj = tempfile.TemporaryDirectory(prefix="phyml_from_vcf_")
        workdir = Path(tempdir_obj.name)

    # Optional: build BED from GTF/GFF
    bed_from_features = None
    if any([args.anno, args.feature_type, args.feature_id, args.feature_id_file]):
        # require anno, feature_type, and some IDs
        if not args.anno or not args.feature_type:
            sys.exit("ERROR: Feature selection requires --anno and --feature-type plus --feature-id/--feature-id-file.")
        ids = load_feature_ids(args.feature_id, args.feature_id_file)
        bed_from_features = str(workdir / "features.bed")
        gtf_gff_to_bed(
            anno_path=args.anno,
            feature_type=args.feature_type,
            feature_attr=args.feature_attr,
            ids_set=ids,
            bed_out_path=bed_from_features,
            merge=args.merge_bed,
            pad=max(0, args.pad or 0)
        )
        print(f"[info] Built feature BED: {bed_from_features}")

    filtered_vcf = str(workdir / "filtered.vcf")
    phylip_out  = str(out_prefix.with_suffix(".phy"))

    # 1) bcftools view
    bcft_cmd = build_bcftools_cmd(args, filtered_vcf, bed_from_features=bed_from_features)

    # 2) vcf2phylip.py
    v2p_cmd = [args.vcf2phylip_path, "-i", filtered_vcf, "-f"]
    if args.vcf2phylip_extra.strip():
        v2p_cmd += shlex.split(args.vcf2phylip_extra)

    if args.dry_run:
        print("Dry run commands:\n")
        print("WORKDIR:", workdir)
        print("$ " + " ".join(shlex.quote(c) for c in bcft_cmd))
        if bed_from_features:
            print(f"# (feature BED from {args.anno}): {bed_from_features}")
        print("$ (cd", shlex.quote(str(workdir)), "&&", " ".join(shlex.quote(c) for c in v2p_cmd), ")")
        print("$ phyml -i", shlex.quote(str(out_prefix.with_suffix('.phy'))), args.phyml_args)
        return

    # Execute steps
    run(bcft_cmd, log_prefix="[bcftools] ")
    run(v2p_cmd, cwd=str(workdir), log_prefix="[vcf2phylip] ")

    # Find a PHYLIP file in workdir (prefer the newest *.phy)
    phy_candidates = sorted(workdir.glob("*.phy"), key=lambda p: p.stat().st_mtime, reverse=True)
    if not phy_candidates:
        sys.exit("ERROR: Could not find a .phy file after running vcf2phylip.py")
    phy_src = phy_candidates[0]

    # Copy/move to desired output name
    if Path(phylip_out) != phy_src:
        Path(phylip_out).write_bytes(phy_src.read_bytes())

    # 3) PhyML
    phyml_full_cmd = ["phyml", "-i", phylip_out] + shlex.split(args.phyml_args)
    run(phyml_full_cmd, log_prefix="[phyml] ")

    # Clean up (optionally)
    if tempdir_obj and not args.keep_intermediates:
        tempdir_obj.cleanup()
    elif tempdir_obj and args.keep_intermediates:
        print(f"[info] Intermediate files kept at: {workdir}")

    print("\nDone.")
    print(f"PHYLIP alignment: {phylip_out}")
    print("PhyML outputs (typical):")
    print(f"  {phylip_out}_phyml_tree.txt")
    print(f"  {phylip_out}_phyml_stats.txt")
    print(f"  {phylip_out}_phyml_boot_trees.txt  (if bootstraps requested)")

if __name__ == "__main__":
    main()
