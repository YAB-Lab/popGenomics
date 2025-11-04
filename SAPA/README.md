# SAPA v2.2 — Sequence Alignment & PAML Analysis

SAPA (**S**equence **A**lignment and **P**AML **A**nalysis) aligns coding sequences, translates them to protein, computes pairwise **Ka/Ks**, and optionally runs **PAML codeml** and **HyPhy** selection tests. It supports restricting samples, adding outgroups via an orthology map, and exporting/viewing alignments.

> **Language:** Perl  
> **Version:** 2.2

---

## Features

- **Input once, analyze many**: Provide a CDS FASTA and either a single `--gene_id` or a list via `--gene_list`.
- **Automated alignment pipeline**: Protein alignment (*ClustalW*) → codon-aware back-translation to DNA via BioPerl (`aa_to_dna_aln`).
- **Ka/Ks**: Pairwise dN/dS using `Bio::Tools::Run::Phylo::PAML::Codeml`.
- **PAML codeml**: Null model by default; optional **branch** and **branch-site** tests per selected foreground taxa.
- **HyPhy**: Run a specified HyPhy method on the DNA alignment (requires a pruned tree).
- **Sample control**: Include/exclude sample sets; orthology-based outgroup inclusion.
- **Convenience**: Save or auto-clean intermediate alignment files; optional macOS Geneious integration for quick viewing.

---

## Requirements

### Perl modules (BioPerl)
- `Bio::Tools::Run::Phylo::PAML::Codeml`
- `Bio::Tools::Run::Alignment::Clustalw`
- `Bio::Align::Utilities` (for `aa_to_dna_aln`)
- `Bio::SeqIO`, `Bio::AlignIO`

Install via CPAN (example):
```bash
cpanm Bio::Perl Bio::Tools::Run::Phylo::PAML::Codeml Bio::Tools::Run::Alignment::Clustalw Bio::Align::Utilities
```

### External tools / scripts (must be on `PATH`)
- **ClustalW** (for protein alignment)
- **PAML** (codeml)
- **PHAST** `tree_doctor` (tree pruning)
- **HyPhy** (optional; for `--run_hyphy`)
- `fastaSortByName.pl` (header-preserving sort of FASTA)
- `fastagrep.pl` (sequence grep; used for sample inclusion/exclusion)
- `paml_prep.pl` (FASTA→PHYLIP & mapping; ask maintainer)
- **Geneious** (optional; only if using `--view_*` on macOS)

### Input FASTA header format

CDS file must contain all sequences. Headers must include the species tag:
```
>GENE_ID ... species=<SPECIES_ID> [line=<LINE_ID>]
ATG...
```
The script normalizes headers internally (replacing spaces and using `species=...` and optional `line=...`).

---

## Usage

```bash
perl sapa.pl   --CDS_file CDS.fasta   (--gene_list genes.txt | --gene_id FBtr0000001)   --output Results/SAPA   [--show_samples]   [--include_samples | --exclude_samples --samples_file samples.txt]   [--include_outgroups --orthology_file orthology.tsv]   [--calculate_KaKs]   [--run_PAML --tree_file species.tree [--branch | --branchSite] --restrict_samples foreground.txt]   [--run_hyphy --tree_file species.tree --hyphy_method FEL]   [--save_alignments]   [--view_DNA_alignment] [--view_Protein_alignment]
```

**Required**
- `--CDS_file|-c` : FASTA of CDS from multiple species  
- One of `--gene_list|-g` (file with one gene_id per line) **or** `--gene_id|-i`  
- `--output|-o` : Output directory

**Helpful utilities**
- `--show_samples` : Print discovered species IDs from `CDS_file` and exit

**Sample control**
- `--include_samples` or `--exclude_samples` with `--samples_file <file>`
  - `samples.txt` should list either line IDs (preferred if present in headers) or species IDs (one per line)

**Outgroups / orthology**
- `--include_outgroups` : Include orthologs that may use different gene IDs
- `--orthology_file <tsv>` : Two-column file: `transcript_id<TAB>orthologue_id`

**PAML & HyPhy**
- `--run_PAML` requires `--tree_file`
  - Defaults to **NULL model** (`model=0`, `NSsites=0`)
  - Enable **branch model** with `--branch`
  - Enable **branch-site** models with `--branchSite` (H0/H1 will be tested)
  - `--restrict_samples <file>` : Foreground taxa list; IDs must match tree & FASTA
- `--run_hyphy` requires `--tree_file` and `--hyphy_method`
  - Example methods: `FEL`, `MEME`, `BUSTED`, etc. (method names depend on your HyPhy installation)

**Viewing / outputs**
- `--save_alignments` : Keep `.afa` alignments
- `--view_DNA_alignment`, `--view_Protein_alignment` : Open in Geneious (macOS)

---

## What the script does

1. **Choose transcripts:** from `--gene_id` or line-wise in `--gene_list`
2. **Extract sequences:**
   - If `--include_outgroups`: use orthology map to collect orthologs
   - Else: use `fastagrep.pl` to grab target transcript
   - Apply `--include_samples` or `--exclude_samples` using `--samples_file`
3. **Align:**
   - Translate CDS to protein → align with **ClustalW**  
   - Back-project alignment to codon space via `aa_to_dna_aln`
   - Export `*.aln_aa.afa` and `*.aln_dna.afa`
4. **Analyses (optional):**
   - **Ka/Ks**: pairwise dN/dS across taxa; aggregated into `KaKs.tmp/` and then `KaKs.txt`
   - **PAML**:
     - Prune `--tree_file` to present taxa using `tree_doctor`
     - Prepare `.phy` via `paml_prep.pl`
     - Run codeml:
       - NULL model (always)
       - **Branch** model if `--branch`
       - **Branch-site** models H0/H1 if `--branchSite`
     - Collect outputs in `PAML.output/`
   - **HyPhy**:
     - Prune tree and run `hyphy <method>` on DNA alignment
     - Collect outputs in `HyPhy.output/`
5. **Finalize outputs:** move to the requested `--output` directory; clean temporary files unless requested to keep.

---

## Examples

### 1) Show available samples and exit
```bash
perl sapa.pl -c CDS.fasta --show_samples
```

### 2) Single gene with Ka/Ks only
```bash
perl sapa.pl   -c CDS.fasta -i FBtr0330000 -o Results/SAPA   --calculate_KaKs --save_alignments
```

### 3) Gene list; include a sample panel; run PAML null + branch
```bash
perl sapa.pl   -c CDS.fasta -g gene_ids.txt -o Results/SAPA   --include_samples --samples_file keep.txt   --run_PAML --tree_file species.tree --branch   --restrict_samples foreground.txt   --save_alignments
```

### 4) Include outgroups via orthology map; run branch-site and HyPhy
```bash
perl sapa.pl   -c CDS.fasta -i FBtr0330000 -o Results/SAPA   --include_outgroups --orthology_file orthology.tsv   --run_PAML --tree_file species.tree --branchSite --restrict_samples fg.txt   --run_hyphy --hyphy_method MEME --tree_file species.tree   --save_alignments
```

---

## Outputs

Inside your `--output` directory:

- **Alignments (if saved)**
  - `GENE.DNA_alignment.afa`
  - `GENE.Protein_alignment.afa`
- **Pairwise Ka/Ks**
  - `KaKs.txt` (columns: `TRANSCRIPT SEQ1 SEQ2 Ka Ks Ka/Ks PROT_PERCENTID CDNA_PERCENTID`)
- **PAML**
  - `PAML.output/GENE.out` (NULL)
  - `PAML.output/GENE.<strain>.br.out` (branch models if requested)
  - `PAML.output/GENE.<strain>.brSt.H0.out/H1.out` (branch-site)
- **HyPhy**
  - `HyPhy.output/GENE.<METHOD>.result.txt`

Temporary/intermediate files are routinely cleaned (e.g., `.phy`, pruned trees), except when needed for downstream steps.

---

## Tips & Notes

- **Tree matching**: Species IDs in `--tree_file` **must** match FASTA (post-normalization); headers are normalized to `>species_or_lineID`.
- **Foreground list**: `--restrict_samples` must contain IDs that exist **both** in the pruned tree and the FASTA alignment.
- **Stop codons**: Proteins with internal `*` are warned; trailing `*` are stripped before alignment.
- **ClustalW**: Ensure `clustalw` is installed and executable on `PATH`.

---

## Troubleshooting

- **“Need at least 2 CDS sequences to proceed”**  
  Verify `CDS_file` contains at least two taxa for the given gene.

- **PAML fails or empty outputs**  
  Check that `paml_prep.pl` is available; confirm `species.tree` contains the same taxa names as your alignment headers **after** script header normalization.

- **HyPhy errors**  
  Confirm method name and that HyPhy can read your alignment format; ensure `tree_doctor` pruned tree includes the same taxa as the alignment.

- **Ka/Ks all NA**  
  Often due to too-short alignments or extreme divergence. Inspect `*.aln_dna.afa` and pairwise identities reported in `KaKs.txt`.

---

## Reproducibility

Record exact versions and commands:
```bash
perl -v
clustalw -version
codeml 2>&1 | head -n 5
hyphy --version
tree_doctor -h | head -n 1
```
Keep your run command in a log:
```bash
perl sapa.pl ... |& tee run.log
```

---

## License

Add a license file (MIT, BSD-3, GPL-3, etc.) to the repository.

---

## Acknowledgements

- BioPerl developers and maintainers  
- PAML (codeml), PHAST, HyPhy teams  
- Utility scripts authors: `fastagrep.pl`, `fastaSortByName.pl`, `paml_prep.pl`
