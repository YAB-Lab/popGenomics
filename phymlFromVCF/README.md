# phyml-from-vcf

Build a phylogenetic tree with **PhyML** from an (optionally filtered) **VCF**.

This script chains together:

1. *(optional)* **GTF/GFF → BED** (feature-based region selection)
2. **bcftools view** → subset by region / BED / sample / expressions
3. **vcf2phylip.py** → PHYLIP alignment
4. **phyml** → tree inference


---

## Highlights

- **Feature-aware region selection** directly from a GTF/GFF (by feature type and attribute ID), with optional **merge** and **padding**.
- Flexible **region** selection via `-r` (chr:start-end) or `-R` (BED).
- Clean **sample include/exclude** mechanics (by list or file).
- Pass-through **bcftools expressions** for powerful filtering.
- Works in a **temporary workdir** by default; preserve intermediates with `--keep-intermediates`.
- **Dry-run** mode to preview all commands.

---

## Requirements

- Python ≥3.8
- [bcftools](http://www.htslib.org/doc/bcftools.html)
- [vcf2phylip](https://github.com/edgardomortiz/vcf2phylip) (script path configurable)
- [PhyML](http://www.atgc-montpellier.fr/phyml/)

Recommended: bgzipped + indexed VCF/BCF (`.vcf.gz` + `.tbi` or `.csi`).

### Quick Conda setup (example)

```bash
mamba create -n phyml_from_vcf -c conda-forge -c bioconda python=3.11 bcftools phyml
conda activate phyml_from_vcf
# Install vcf2phylip (choose one)
# 1) pip (if available): pip install vcf2phylip
# 2) clone script: git clone https://github.com/edgardomortiz/vcf2phylip && export PATH=$PWD/vcf2phylip:$PATH
```

---

## Usage

```bash
python phyml_from_vcf.py   -i INPUT.vcf.gz   -o Results/trees/mytree   [--workdir WORKDIR] [--keep-intermediates]   [-r chr1:1-100000 -r contig2]   [-R regions.bed]   [--anno genes.gtf --feature-type gene --feature-attr gene_id --feature-id FBgn0000008]   [--merge-bed --pad 200]   [--samples S1 S2 ... | --samples-file samples.txt]   [--exclude-samples X Y ... | --exclude-samples-file exclude.txt]   [--force-samples]   [-in 'QUAL>30 && DP>10' | -ex 'AC==0']   [--bcftools-args "--min-ac 1"]   [--vcf2phylip-path vcf2phylip.py]   [--vcf2phylip-extra "--min-taxa 4 --max-missing 0.2"]   [--phyml-args "-d nt -m GTR -b 100"]   [--dry-run]
```

Outputs (with `-o Results/trees/mytree`):
- `mytree.phy` — PHYLIP alignment
- `mytree.phy_phyml_tree.txt` — PhyML ML tree (and optionally `*_boot_trees.txt` if bootstraps requested)
- `mytree.phy_phyml_stats.txt` — PhyML statistics

> By default, intermediates live in a temporary directory and are removed unless `--keep-intermediates` is set.

---

## Common Workflows

### 1) Basic tree (no extra filtering)
```bash
python phyml_from_vcf.py   -i cohort.vcf.gz   -o Results/trees/basic_tree
```

### 2) Select genomic windows or contigs
```bash
# Single region and entire contig
python phyml_from_vcf.py   -i cohort.vcf.gz   -o Results/trees/chr2_window   -r chr2:1-2500000 -r chr3

# Multiple windows via BED
python phyml_from_vcf.py   -i cohort.vcf.gz   -o Results/trees/bed_subset   -R windows.bed
```

### 3) Feature-based selection from GTF/GFF
Extract intervals for specific features (e.g., genes), optionally **merge** overlapping intervals and **pad**.

```bash
# Select one or more genes by ID from a GTF
python phyml_from_vcf.py   -i cohort.vcf.gz   -o Results/trees/gene_panel   --anno genes.gtf   --feature-type gene   --feature-attr gene_id   --feature-id FBgn0000008 --feature-id FBgn0000017   --merge-bed --pad 200
```

Or provide a file of IDs (one per line):
```bash
python phyml_from_vcf.py   -i cohort.vcf.gz   -o Results/trees/gene_list   --anno genes.gtf   --feature-type exon   --feature-attr transcript_id   --feature-id-file transcripts.txt
```

> **Mutual exclusivity:** Choose **either** feature-based selection (`--anno/--feature-*`) **or** raw `-r/-R` region flags (not both).

### 4) Sample (sub)sets
```bash
# Include only specific samples
python phyml_from_vcf.py -i cohort.vcf.gz -o Results/trees/subset --samples A B C

# From file
python phyml_from_vcf.py -i cohort.vcf.gz -o Results/trees/subset --samples-file keep.txt

# Exclude a set
python phyml_from_vcf.py -i cohort.vcf.gz -o Results/trees/no_outgroups --exclude-samples OG1 OG2
```

### 5) bcftools expressions and extra args
```bash
# Include variants with quality/depth thresholds
python phyml_from_vcf.py   -i cohort.vcf.gz -o Results/trees/highconf   -in 'QUAL>30 && INFO/DP>10'

# Exclude monomorphic / failed sites, plus extra bcftools args
python phyml_from_vcf.py   -i cohort.vcf.gz -o Results/trees/polymorphic   -ex 'AC==0 || F_MISSING>0.5'   --bcftools-args "--min-ac 1"
```

### 6) Control PhyML
```bash
# Nucleotide model with bootstraps
python phyml_from_vcf.py   -i cohort.vcf.gz -o Results/trees/gtr100   --phyml-args "-d nt -m GTR -b 100"

# Amino-acid mode (if input alignment is AA; vcf2phylip emits NT by default)
python phyml_from_vcf.py   -i cohort.vcf.gz -o Results/trees/lg   --phyml-args "-d aa -m LG -b 0"
```

### 7) Dry-run to preview commands
```bash
python phyml_from_vcf.py   -i cohort.vcf.gz   -o Results/trees/test   --anno genes.gtf --feature-type gene --feature-id FBgn0000008   --dry-run
```

---

## Notes & Behavior

- **Intermediates** are written to a temporary work directory unless `--workdir` is provided. Use `--keep-intermediates` to retain them.
- The script **prefers** the newest `*.phy` produced by `vcf2phylip.py` in the workdir.
- Feature attributes:
  - GTF-style: `key "value";`
  - GFF3-style: `key=value;`
  - Choose with `--feature-attr` (e.g., `gene_id`, `transcript_id`, `ID`, `Name`).

**Mutual exclusivity / conflicts handled with explicit errors:**
- `--samples` **vs** `--samples-file`
- `--exclude-samples` **vs** `--exclude-samples-file`
- `-i/-ex` include/exclude expressions
- Feature-based selection **vs** `-r/-R` region flags

---

## Troubleshooting

- **"Could not find a .phy file after running vcf2phylip.py"**  
  Ensure `vcf2phylip.py` is on your `PATH` or set `--vcf2phylip-path`, and that your bcftools filtering isn’t removing all sites.

- **bcftools errors about samples / regions**  
  Verify mutual exclusivity rules and that sample names match the VCF header.

- **Empty BED from features**  
  Check `--feature-type`, `--feature-attr`, and the IDs you supplied (`--feature-id` or `--feature-id-file`). Use `--merge-bed`/`--pad` if appropriate.

- **PhyML model/argument issues**  
  Review `phyml --help` and ensure your `--phyml-args` match the desired datatype and model.

---

## Reproducibility tip

Capture your exact command and software versions:

```bash
bcftools --version
phyml --version
python --version
python phyml_from_vcf.py ... --dry-run 2>&1 | tee run.log
```

---

## License

Choose a license you prefer (e.g., MIT, BSD-3, GPL-3). Add a `LICENSE` file to your repository.

---

## Acknowledgements

- **bcftools** (HTSlib team)
- **vcf2phylip** (Edgar E. Ortiz)
- **PhyML** (Guindon *et al.*)
