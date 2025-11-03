#!/usr/bin/env python3

import argparse
import subprocess
import os
import pandas as pd

def run_command(command, step):
    """Run a shell command with step description and handle errors."""
    print(f"[INFO] Running step: {step}")
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing command in step '{step}': {e}")
        exit(1)

def gzip_and_index_vcf(input_vcf):
    """Compress and index the VCF file if it is not already gzipped."""
    if not input_vcf.endswith(".gz"):
        print(f"[INFO] Gzipping and indexing VCF file: {input_vcf}")
        gzipped_vcf = f"{input_vcf}.gz"
        run_command(f"bgzip -c {input_vcf} > {gzipped_vcf}", "Gzip VCF")
        run_command(f"tabix -p vcf {gzipped_vcf}", "Index VCF")
        return gzipped_vcf
    return input_vcf

def filter_vcf_coding(input_vcf, output_vcf, region=None, bed_file=None, bcftools_filter=None):
    """Filter the VCF file using bcftools based on the given region, BED file, or additional filter options."""
    print("[INFO] Filtering VCF file for coding regions...")
    input_vcf = gzip_and_index_vcf(input_vcf)
    command = f"bcftools view {input_vcf}"
    if region:
        command += f" -r {region}"
    if bed_file:
        command += f" -R {bed_file}"
    if bcftools_filter:
        command += f" {bcftools_filter}"
    command += f" -o {output_vcf}"
    run_command(command, "VCF Filtering")

def filter_vcf_noncoding(input_vcf, output_vcf1, output_vcf2, region=None, region2=None, bed_file=None, bed_file2=None, bcftools_filter=None):
    """Filter the VCF file using bcftools based on the given region, BED file, or additional filter options."""
    print("[INFO] Filtering VCF file for noncoding regions...")
    input_vcf = gzip_and_index_vcf(input_vcf)
    command1 = f"bcftools view {input_vcf}"
    command2 = f"bcftools view {input_vcf}"
    if region and region2:
        command1 += f" -r {region}"
        command2 += f" -r {region2}"
    elif bed_file and bed_file2:
        command1 += f" -R {bed_file}"
        command2 += f" -R {bed_file2}"
    if bcftools_filter:
        command1 += f" {bcftools_filter}"
        command2 += f" {bcftools_filter}"
    command1 += f" -o {output_vcf1}"
    command2 += f" -o {output_vcf2}"
    run_command(command1, "VCF Filtering")
    run_command(command2, "VCF Filtering")

def filter_variants_coding(input_vcf, output_d, output_d0, output_p, output_p0):
    """Extract fixed and segregating sites using SnpSift."""
    print("[INFO] Extracting variant types using SnpSift...")
    run_command(f"SnpSift filter \"(ANN[*].EFFECT has 'missense_variant') && (AF = 1.0)\" {input_vcf} > {output_d}", "Filtering missense fixed variants")
    run_command(f"SnpSift filter \"(ANN[*].EFFECT has 'synonymous_variant') && (AF = 1.0)\" {input_vcf} > {output_d0}", "Filtering synonymous fixed variants")
    run_command(f"SnpSift filter \"(ANN[*].EFFECT has 'missense_variant') && (AF < 1.0)\" {input_vcf} > {output_p}", "Filtering missense segregating variants")
    run_command(f"SnpSift filter \"(ANN[*].EFFECT has 'synonymous_variant') && (AF < 1.0)\" {input_vcf} > {output_p0}", "Filtering synonymous segregating variants")

def filter_variants_noncoding(input_vcf1, input_vcf2, output_d, output_d0, output_p, output_p0):
    """Extract fixed and segregating sites using SnpSift."""
    print("[INFO] Extracting variant types using SnpSift...")
    run_command(f"SnpSift filter \"AF = 1.0\" {input_vcf1} > {output_d}", "Filtering test region fixed variants")
    run_command(f"SnpSift filter \"AF = 1.0\" {input_vcf2} > {output_d0}", "Filtering neutral fixed variants")
    run_command(f"SnpSift filter \"AF < 1.0\" {input_vcf1} > {output_p}", "Filtering test region segregating variants")
    run_command(f"SnpSift filter \"AF < 1.0\" {input_vcf2} > {output_p0}", "Filtering neutral segregating variants")

def count_variants(vcf_file):
    """Count the number of lines in a VCF file excluding header lines."""
    print(f"[INFO] Counting variants in {vcf_file}...")
    command = f"grep -vc '^#' {vcf_file}"
    return int(subprocess.check_output(command, shell=True).strip())

def extract_allele_frequencies(input_vcf, output_file):
    """Extract allele frequency data using SnpSift extractFields."""
    print(f"[INFO] Extracting allele frequencies from {input_vcf}...")
    command = (
        f"SnpSift extractFields {input_vcf} "
        f"CHROM POS AF ANN[*].GENE ANN[*].EFFECT ANN[*].IMPACT ANN[*].AA_POS "
        f"-s ',' -e '.' > {output_file}"
    )
    run_command(command, "Extracting allele frequencies")

def bin_frequencies(freq_file, bin_size=0.05):
    """Bin allele frequencies into intervals and count occurrences, keeping the highest alternative frequency."""
    print(f"[INFO] Binning allele frequencies from {freq_file}...")
    df = pd.read_csv(freq_file, sep='\t')
    
    # Handle cases with multiple alternative frequencies
    df["AF"] = df["AF"].apply(lambda x: max(map(float, str(x).split(','))) if ',' in str(x) else float(x))
    
    df["bin"] = (df["AF"] // bin_size) * bin_size
    binned_counts = df.groupby("bin").size().reset_index(name="count")
    return binned_counts

def run_asymptotic_mk(r_script_path, d0, d, mk_table, xlow, xhigh, output_file):
    """Run the asymptoticMK R script with the given inputs."""
    print("[INFO] Running asymptoticMK in R...")
    command = (
        f"Rscript {r_script_path} {d0} {d} {xlow} {xhigh} {mk_table} {output_file}"
    )
    run_command(command, "Run asymptoticMK in R")

def main():
    parser = argparse.ArgumentParser(description="Wrapper script for running the asymptotic McDonald-Kreitman test using a variant-annotated VCF.")
    parser.add_argument("--mode", choices=["coding", "noncoding"], required=True, help="Analysis mode: 'coding' or 'noncoding'")
    parser.add_argument("-i", "--input_vcf", required=True, help="Input annotated VCF file")
    parser.add_argument("-o", "--output_prefix", required=True, help="Prefix for output files")
    parser.add_argument("-r", "--region", help="Region to analyze ([CHR:min-max]; optional for 'coding' mode; required for 'noncoding' mode)")
    parser.add_argument("-r2", "--region2", help="Second region for filtering in noncoding mode (only required for 'noncoding' mode)")
    parser.add_argument("-R", "--bed_file", help="BED file for filtering (optional; required if mode = 'noncoding')")
    parser.add_argument("-R2", "--bed_file2", help="Second BED file for filtering in noncoding mode (required if mode = 'noncoding')")
    parser.add_argument("--bcftools_filter", help="Additional bcftools filtering options (optional)")
    parser.add_argument("--keep_intermediates", default=False, action="store_true", help="Keep intermediate files")
    # parser.add_argument("--snpsift_path", default=".", help="Path to SnpSift directory (default: current directory)")
    parser.add_argument("--r_script_path", help="Path to the asymptoticMK R script (optional)")
    parser.add_argument("--xlow", type=float, default=0.1, help="Lower bound for frequency spectrum (default: 0.1)")
    parser.add_argument("--xhigh", type=float, default=0.9, help="Upper bound for frequency spectrum (default: 0.9)")
    parser.add_argument("--run_asymptoticMK", action="store_true", help="Run asymptoticMK R script if specified")
    args = parser.parse_args()


    print("[INFO] Starting McDonald-Kreitman test filtering script...")
    
    results_d = f"{args.output_prefix}_results_d.vcf"
    results_d0 = f"{args.output_prefix}_results_d0.vcf"
    results_p = f"{args.output_prefix}_results_p.vcf"
    results_p0 = f"{args.output_prefix}_results_p0.vcf"
    
    if args.mode == "coding":
        print("[INFO] Running in coding sequence mode...")
        filtered_vcf = f"{args.output_prefix}_filtered.vcf"
        filter_vcf_coding(args.input_vcf, filtered_vcf, args.region, args.bed_file)

        filter_variants_coding(filtered_vcf, results_d, results_d0, results_p, results_p0)

    elif args.mode == "noncoding":
        print("[INFO] Running in non-coding sequence mode...")

        filtered_vcf1 = f"{args.output_prefix}_filtered_test_region.vcf"
        filtered_vcf2 = f"{args.output_prefix}_filtered_neutral_region.vcf"
        filter_vcf_noncoding(args.input_vcf, filtered_vcf1, filtered_vcf2, args.region, args.region2, args.bed_file, args.bed_file2)

        filter_variants_noncoding(filtered_vcf1, filtered_vcf2, results_d, results_d0, results_p, results_p0)

    d = count_variants(results_d)
    d0 = count_variants(results_d0)
    
    with open(f"{args.output_prefix}_variant_counts.txt", "w") as f:
        f.write(f"d: {d}\nd0: {d0}\n")

    extract_allele_frequencies(results_p, f"{args.output_prefix}_freq_p.txt")
    extract_allele_frequencies(results_p0, f"{args.output_prefix}_freq_p0.txt")

    binned_p = bin_frequencies(f"{args.output_prefix}_freq_p.txt")
    binned_p0 = bin_frequencies(f"{args.output_prefix}_freq_p0.txt")
    
    final_table = pd.merge(binned_p, binned_p0, on="bin", how="outer").fillna(0)
    final_table.columns = ["x", "p", "p0"]
    mk_table = f"{args.output_prefix}_mk_table.txt"
    final_table.to_csv(mk_table, sep='\t', index=False)

    if args.run_asymptoticMK and args.r_script_path:
        output_file = f"{args.output_prefix}_mk_results.txt"
        run_asymptotic_mk(args.r_script_path, d0, d, mk_table, args.xlow, args.xhigh, output_file)

    print("[INFO] McDonald-Kreitman test filtering script completed successfully.")

if __name__ == "__main__":
    main()