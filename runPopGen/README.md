# MK–PopGen Alignment Pipeline (Perl)

This script builds codon-aware alignments from CDS FASTA, then runs **McDonald–Kreitman (MK) tests** and/or **population genetic statistics** (Tajima’s D, Fu & Li’s D*, F*) with optional synonymous vs. nonsynonymous partitioning.

> **Language:** Perl  
---

## What it does

1. **Extract sequences** for a target transcript (or list) from a CDS FASTA; optionally include orthologs via a mapping file.  
2. **Translate → align proteins (ClustalW) → back-project to codon alignment** using BioPerl (`aa_to_dna_aln`).  
3. **Outputs alignments** (`*.aln_aa.afa`, `*.aln_dna.afa`).  
4. (Optional) **MK test** (polarized or unpolarized) per transcript and collates a final `MKout.txt`.  
5. (Optional) **Population genetics**: computes pi, theta, Tajima’s D, Fu & Li’s D*, F* on all sites; optionally split by **synonymous vs. nonsynonymous** sites, producing `PopGen.all.out.txt`, `PopGen.Syn.out.txt`, `PopGen.nonSyn.out.txt`.

---

## Requirements

### Perl modules (BioPerl)
- `Bio::Tools::Run::Alignment::Clustalw`
- `Bio::Align::Utilities` (for `aa_to_dna_aln`)
- `Bio::PopGen::IO`, `Bio::PopGen::Statistics`, `Bio::PopGen::Utilities`
- `Bio::SeqIO`, `Bio::AlignIO`

Install via CPAN (example):
```bash
cpanm Bio::Perl Bio::Tools::Run::Alignment::Clustalw Bio::Align::Utilities Bio::PopGen Bio::SeqIO Bio::AlignIO
```

### External tools / helper scripts (must be on `PATH`)
- **ClustalW** (alignment)
- `fastagrep.pl` (sequence extraction by ID; supports `-X`, `-f`, `-v`)
- `MK.pl` and `MK_trim.pl` (MK test + output cleanup; MK.pl originally by Alisha Holloway)
- `tajimasD_calculator.pl` (compute pi, theta, Tajima’s D, Fu & Li’s D*, F*)
- `splitSynNonsyn.py` (split alignment into synonymous vs. nonsynonymous sites)
- `fastaSortByName.pl` (optional convenience script)
- **Geneious** (optional; macOS only, for `--view_*`)

### Input FASTA header format
All CDS sequences must contain `species=` and optionally `line=` in the header:
```
>TRANSCRIPT_ID ... species=<SPECIES_ID> line=<LINE_ID>
ATG...
```
The script normalizes headers (e.g., replaces spaces; retains `species=` and `line=` where present).

---

## Usage

```bash
perl mk_popgen_pipeline.pl   --CDS_file CDS.fasta   (--gene_list genes.txt | --gene_id FBtr0000001)   --output Results/MKPopGen   [--show_samples]   [--include_samples | --exclude_samples --samples_file samples.txt]   [--include_outgroups --orthology_file orthology.tsv]   [--save_alignments]   [--view_DNA_alignment] [--view_Protein_alignment]   [--MK_test --pol pol --outgroup1 SP1 --outgroup2 SP2 --ingroup ING]   [--MK_test --pol unpol --outgroup OUT --ingroup ING]   [--PopGen [--split_codon_sites]]
```

### Required
- `--CDS_file|-c` : CDS FASTA containing all sequences  
- One of `--gene_list|-g` (file; one transcript per line) **or** `--gene_id|-i`  
- `--output|-o` : Output directory

### Helpful utilities
- `--show_samples` : Prints discovered `species=`/`line=` identifiers from headers and exits

### Sample selection
- `--include_samples` or `--exclude_samples` used with `--samples_file <file>`  
  `samples.txt` should list **line IDs** (preferred) or **species IDs**, one per line

### Orthologs
- `--include_outgroups` uses `--orthology_file` (TSV: `transcript_id<TAB>orthologue_id`) to add orthologs that might use different IDs within/between species

### MK test
Choose **one**:
- Polarized: `--MK_test --pol pol --outgroup1 SP1 --outgroup2 SP2 --ingroup ING`
- Unpolarized: `--MK_test --pol unpol --outgroup OUT --ingroup ING`

Outputs per transcript are written into `OUTPUT/MK.tmp/` and collated into `OUTPUT/MKout.txt`.

### PopGen stats
- `--PopGen` : Computes pi, theta, Tajima’s D, Fu & Li’s D*, F* on **all sites**
- `--split_codon_sites` : Additionally compute stats on **Syn** and **nonSyn** partitions, producing `PopGen.Syn.out.txt` and `PopGen.nonSyn.out.txt`

---

## Examples

### 1) Show samples and exit
```bash
perl mk_popgen_pipeline.pl -c CDS.fasta --show_samples
```

### 2) MK (unpolarized) for a single gene; keep alignments
```bash
perl mk_popgen_pipeline.pl   -c CDS.fasta -i FBtr0330000 -o Results/MKPopGen   --MK_test --pol unpol --outgroup Dvir --ingroup Dnov   --save_alignments
```

### 3) MK (polarized) + PopGen for a gene list; include only a panel of lines
```bash
perl mk_popgen_pipeline.pl   -c CDS.fasta -g gene_ids.txt -o Results/MKPopGen   --include_samples --samples_file keep.txt   --MK_test --pol pol --outgroup1 Dvir --outgroup2 Dlum --ingroup Dnov   --PopGen --split_codon_sites
```

### 4) Add orthologs via map; PopGen only
```bash
perl mk_popgen_pipeline.pl   -c CDS.fasta -i FBtr0330000 -o Results/MKPopGen   --include_outgroups --orthology_file orthology.tsv   --PopGen
```

---

## Outputs

Inside your `--output` directory:

- **Alignments (if saved)**
  - `TRANSCRIPT.DNA_alignment.afa`
  - `TRANSCRIPT.Protein_alignment.afa`
- **MK test**
  - `MKout.txt` (columns: `TRANSCRIPT  NS_POLY  S_POLY  NS_FIX  S_FIX  codons  final_NI  alpha  final_FET`)
- **Population genetics**
  - `PopGen.all.out.txt` — all sites
  - `PopGen.Syn.out.txt` — synonymous sites (if `--split_codon_sites`)
  - `PopGen.nonSyn.out.txt` — nonsynonymous sites (if `--split_codon_sites`)

Temporary directories used during runs:
- `MK.tmp/` (removed after collation)
- `PopGen.tmp/` (removed after collation)

---

## Notes & Behavior

- **At least two sequences** per transcript are required to align and analyze.  
- Header cleanup removes path-like suffixes from FASTA IDs in alignment outputs.  
- When splitting sites (`--split_codon_sites`), temporary files `TRANSCRIPT.Syn.fasta` and `TRANSCRIPT.nonSyn.fasta` are cleaned after use.  
- macOS-only Geneious viewing is available via `--view_*`.

---

## Troubleshooting

- **“Need at least 2 CDS sequences to proceed”**  
  Your CDS FASTA has only one line/species for that transcript. Add more sequences or skip it.

- **MK outputs empty**  
  Check your polarization settings match your available outgroups (`--pol pol` needs **two** outgroups; `--pol unpol` needs **one**).

- **PopGen stats NA or extreme**  
  May occur with very short alignments or low diversity; inspect `*.aln_dna.afa`.

- **Sample include/exclude not working**  
  Ensure `samples.txt` contains exactly the `line` IDs (or species IDs) present in the normalized headers (the script turns `line=` into `_` in names).

---

## Reproducibility

Record versions and keep the run command:
```bash
perl -v
clustalw -version
MK.pl -h | head -n 5
tajimasD_calculator.pl -h | head -n 5
python --version && splitSynNonsyn.py -h 2>/dev/null | head -n 1
perl mk_popgen_pipeline.pl ... |& tee run.log
```

---

## License

Add a `LICENSE` file (MIT, BSD-3, GPL-3, etc.) to your repo.

---

## Acknowledgements

- BioPerl developers and maintainers  
- MK.pl & MK_trim.pl authors (original MK.pl by **Alisha Holloway**)  
- Authors of `fastagrep.pl`, `tajimasD_calculator.pl`, `splitSynNonsyn.py`  
