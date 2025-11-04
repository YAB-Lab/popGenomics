# GTF → VCF → PHYLIP (with strand-aware reversal)

Extract a user-specified feature from a **GTF**, filter an all-sites **VCF** to those intervals with **bcftools**, build an alignment via **vcf2phylip**, and—if the feature is on the minus strand—**reverse‑complement** the final alignment. Optionally emit **PHYLIP**, **FASTA**, or **both**.

> **Language:** Python 3.8+  

---

## Overview

1. Parse a GTF for a given `--feature-type` and `--feature-id` (matches `gene_id` or `transcript_id`).  
2. Convert matched records to **BED (0‑based, half‑open)**, **merge** overlaps, and **sort**.  
3. `bcftools view -R` to subset the input VCF to those intervals; then `tabix` index.  
4. Run `vcf2phylip -i <filtered.vcf.gz> -f` to produce an alignment.  
5. Detect strand of the feature from the GTF; if **‘-’**, reverse‑complement the alignment sequences.  
6. Write deterministic outputs to `{out-prefix}.phy` and/or `{out-prefix}.fa`.

> Assumes **sequential PHYLIP** output from `vcf2phylip`. (Interleaved support can be added later.)

---

## Requirements

- Python ≥ 3.8
- `bcftools` in `PATH`
- `tabix` in `PATH` (for VCF indexing)
- `vcf2phylip` in `PATH` (or supply full path via `--vcf2phylip`)

### Quick conda setup (example)

```bash
mamba create -n gtf2phy -c conda-forge -c bioconda python=3.11 bcftools htslib
conda activate gtf2phy
# install vcf2phylip (pick one)
# pip:    pip install vcf2phylip   # if packaged in your index
# source: git clone https://github.com/edgardomortiz/vcf2phylip && export PATH=$PWD/vcf2phylip:$PATH
```

---

## Usage

```bash
python gtf_vcf2phylip.py   --vcf cohort.allsites.vcf.gz   --gtf genes.gtf   --feature-id FBtr0123456   --feature-type CDS   --out-prefix Results/trees/myfeature   [--vcf2phylip vcf2phylip]   [--vcf2phylip-args "--min-taxa 4 --max-missing 0.2"]   [--out-formats phylip|fasta|both]   [--phylip-out Results/trees/myfeature.phy]   [--fasta-out Results/trees/myfeature.fa]   [--keep-temp]
```

**Required**
- `--vcf` : all‑sites VCF/BCF (bgzipped recommended)  
- `--gtf` : GTF with `gene_id`/`transcript_id` attributes  
- `--feature-id` : value to match against `gene_id` *or* `transcript_id`  
- `--feature-type` : GTF feature type to extract (e.g., `exon`, `CDS`, `UTR`, `transcript`)  
- `--out-prefix` : prefix for outputs (`.bed`, filtered VCF, `.phy`/`.fa`)

**Output control**
- `--out-formats` : `phylip` (default), `fasta`, or `both`  
- `--phylip-out` / `--fasta-out` : explicit file paths (override defaults)

**vcf2phylip**
- `--vcf2phylip` : path/command to invoke `vcf2phylip`  
- `--vcf2phylip-args` : raw extra args forwarded to `vcf2phylip` (e.g., `--min-taxa 4 --max-missing 0.2`). The script will try to discover the created alignment file(s).

**Housekeeping**
- `--keep-temp` : keep intermediate BED/VCF files (and the original vcf2phylip outputs)

---

## Coordinate conventions

- **GTF input:** 1‑based, **inclusive**.  
- **BED output:** 0‑based, **half‑open**.  
- Overlapping/adjacent intervals are **merged per chromosome**, then sorted.

---

## Behavior details

- The script matches your `--feature-id` if it equals **either** `gene_id` **or** `transcript_id` for rows whose `feature` equals `--feature-type`.
- Strand is taken from the first matching GTF record. If none is present, **‘+’** is assumed.
- After `vcf2phylip`, the script **discovers** the newest PHYLIP/FASTA file(s) produced and rewrites them deterministically to the requested output paths.
- If the feature is on the **minus strand**, the final alignment sequences are **reverse‑complemented** before writing outputs (PHYLIP and/or FASTA).

---

## Examples

### 1) CDS of a transcript → PHYLIP only
```bash
python gtf_vcf2phylip.py   --vcf cohort.allsites.vcf.gz   --gtf genes.gtf   --feature-id FBtr0330000   --feature-type CDS   --out-prefix Results/trees/FBtr0330000
```

### 2) Exons of a gene → PHYLIP + FASTA (with extra vcf2phylip args)
```bash
python gtf_vcf2phylip.py   --vcf cohort.allsites.vcf.gz   --gtf genes.gtf   --feature-id FBgn0000008   --feature-type exon   --out-prefix Results/trees/FBgn0000008   --out-formats both   --vcf2phylip-args "--min-taxa 4 --max-missing 0.3"
```

### 3) Custom file names; keep intermediates
```bash
python gtf_vcf2phylip.py   --vcf cohort.allsites.vcf.gz   --gtf genes.gtf   --feature-id FBtr0123456   --feature-type transcript   --out-prefix Results/trees/FBtr0123456   --phylip-out Results/trees/FBtr0123456.aln.phy   --fasta-out  Results/trees/FBtr0123456.aln.fa   --out-formats both   --keep-temp
```

---

## Outputs

Given `--out-prefix Results/trees/myfeature`:

- `Results/trees/myfeature.feature.bed` — merged BED intervals  
- `Results/trees/myfeature.filtered.vcf.gz` (+ `.tbi`) — region‑filtered VCF  
- `Results/trees/myfeature.phy` — final PHYLIP (sequential) if `phylip` or `both`  
- `Results/trees/myfeature.fa` — final FASTA if `fasta` or `both`

> If `--keep-temp` is **not** set, the intermediate BED/VCF and raw vcf2phylip outputs are removed after the final files are written.

---

## Troubleshooting

- **“Could not find ‘bcftools’ / ‘tabix’ / ‘vcf2phylip’ in PATH.”**  
  Ensure the tools are installed and on your PATH, or pass the full path to `--vcf2phylip`.

- **“No intervals found …”**  
  Confirm `--feature-type` matches the GTF column 3 (e.g., `CDS`, `exon`) and that your `--feature-id` equals a `gene_id` or `transcript_id` in the attributes column.

- **“PHYLIP header malformed / unequal lengths.”**  
  `vcf2phylip` should produce aligned sequences of identical length. If you pass additional filters in `--vcf2phylip-args`, verify they don’t yield empty/degenerate alignments.

- **Interleaved PHYLIP**  
  This script currently assumes **sequential** PHYLIP. If your `vcf2phylip` emits interleaved, let’s extend the reader.

---

## Reproducibility

Capture the exact tool versions and command:
```bash
bcftools --version
tabix --version
vcf2phylip --help | head -n 5
python --version
python gtf_vcf2phylip.py ... 2>&1 | tee run.log
```

---

## License

Add a `LICENSE` file (MIT, BSD‑3, GPL‑3, etc.) to your repository.

---

## Acknowledgements

- **bcftools/htslib** teams  
- **vcf2phylip** (Edgar E. Ortiz)  
- Everyone maintaining GTF/GFF annotation standards
