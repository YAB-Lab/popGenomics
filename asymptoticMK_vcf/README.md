
# MK Test Wrapper

This Python script is a flexible wrapper for preparing VCF files and associated variant data for McDonald-Kreitman (MK) tests. It supports both traditional **coding region** analysis and **noncoding region** analysis such as enhancers or UTRs.

## Features

- Filters VCF files using `bcftools`, including optional `--region`, `--bed` inputs, and additional `bcftools` filters.
- Handles compressed and uncompressed VCF files automatically.
- Supports two analysis modes:
  - `coding`: Extracts synonymous and missense mutations using `SnpSift`
  - `noncoding`: User supplies two genomic regions to define p/d vs p0/d0 comparisons
- Calculates allele frequencies and optionally runs the asymptotic MK test using an external R script.
- Keeps intermediate files only if specified with a flag.

## Requirements

- Python 3
- `bcftools`, `bgzip`, `tabix`, `SnpSift`, `grep`
- `Rscript` and the `asymptoticMK` R function (optional)

## Usage

```bash
python mk_wrapper.py \
    --input_vcf <input.vcf> \
    --output_prefix <prefix> \
    --mode <coding|noncoding> \
    [--region <chr:start-end>] [--region2 <chr:start-end>] \
    [--bed_file <regions.bed>] [--bed_file2 <regions2.bed>] \
    [--bcftools_filter "<bcftools expression>"] \
    [--snpsift_path <path/to/SnpSift>] \
    [--keep_intermediates] \
    [--run_asymptotic_mk] [--r_script <asymptoticMK.R>]
```

### Example for Coding Mode

```bash
python mk_wrapper.py \
    --input_vcf data.vcf \
    --output_prefix test_run \
    --mode coding \
    --snpsift_path ~/tools/snpEff/scripts/
```

### Example for Noncoding Mode

```bash
python mk_wrapper.py \
    --input_vcf data.vcf \
    --output_prefix enhancer_test \
    --mode noncoding \
    --bed_file enhancer_regions.bed \
    --bed_file2 flanking_regions.bed
```

### Running the Asymptotic MK Test

Add the `--run_asymptotic_mk` flag and specify the R script:

```bash
--run_asymptotic_mk --r_script asymptoticMK_local.R
```

## Output

- Filtered and allele-frequency-counted VCFs
- Frequency tables for each region (for input to asymptotic MK test)
- Optional MK results table if R script is run

## Author

Created and maintained by Yahya Ahmed-Braimah Lab (YAB Lab). For support or contributions, open an issue or submit a PR on the GitHub repository.
