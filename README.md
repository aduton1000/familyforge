# FamilyForge

**Pedigree-aware Mendelian + ACMG-style interpretation and multi-audience reporting on a family callset.**

FamilyForge is a pure-Python clinical variant-interpretation tool. Given one family's
annotated multi-sample VCF, a pedigree, and a curated gene panel, it produces
segregation-checked candidate variants, a pragmatic ACMG-style A/B/C tiering, a
copy-number (CNV) cross-check, genotype-based relatedness/sex QC, and a set of
audience-specific reports (clinician, publication, and letterhead clinical).

It is disease-agnostic — the **panel and pedigree are inputs, never hardcoded** — so any
family runs the same pipeline.

## The clinical problem it solves

A joint-genotyped, annotated VCF tells you *what* variants exist, but not *which* explain a
family's disease. Answering that requires reasoning over the pedigree: does a variant
segregate with the affected members under a plausible inheritance model? Is it rare enough
for that model? Does the gene, consequence, and in-silico evidence support pathogenicity?
FamilyForge automates that first pass, ranks the candidates, shows its evidence trail, and
produces documents you can take to a multidisciplinary review or a case report — while
being honest about what is *inferred* (e.g. phase or de-novo status when a parent is not
sequenced) rather than confirmed.

## Where it sits — the *forge* family

> **CaptureForge** (probe *design*) → **CallForge** (sequence *analysis*) → **FamilyForge** (family *interpretation*)

CallForge produces the annotated joint callset; FamilyForge consumes it (the annotated
`VCF.gz`, and optionally `cnv_calls.tsv` and `provenance.json`) and interprets it for one
family.

## Features

- **Genotype-quality gating** per sample (GQ / DP / allele-balance) before any model runs.
- **Pedigree Mendelian models** — de novo / dominant, autosomal-recessive homozygous,
  compound heterozygous (parental phase inferred when a parent is missing), and X-linked
  recessive — each checked against the affected/unaffected structure.
- **Monozygotic-twin handling** — MZ twins are treated as genetic replicates; a candidate
  must be concordant across the pair, and discordant genotypes are flagged as artifacts.
- **ACMG-style tiering (A/B/C)** — ClinVar, population frequency (dominant vs recessive
  thresholds), predicted consequence, in-silico damage (SIFT/PolyPhen/PhyloP,
  AlphaMissense/REVEL/CADD), conservation, gene disease-mechanism, and configurable
  gene-specific rules, with an explicit evidence trail.
- **CNV cross-check** — deletions overlapping panel genes (priority genes configurable) from
  a CallForge `cnv_calls.tsv`.
- **Genotype-based QC** — pairwise concordance (relatedness / MZ / sample-swap) and sex
  inference from chrX heterozygosity, with no external sites file.
- **Reports** — `clinician_report.{html,md}`, `publication_report.{html,md}`,
  `candidates.json`, plus a standalone letterhead **clinical report** engine.

## Installation

Requires Python ≥ 3.10. `bcftools` on `PATH` is recommended (multiallelic normalisation) but
optional.

```bash
git clone https://github.com/aduton1000/familyforge.git
cd familyforge
python3 -m venv .venv && source .venv/bin/activate
pip install -e .
familyforge --help
```

Or with conda/mamba (also brings `bcftools`):

```bash
mamba create -n familyforge -c conda-forge -c bioconda \
    python=3.11 cyvcf2 pandas numpy jinja2 pyyaml bcftools
mamba activate familyforge && pip install -e .
```

## Quickstart

```bash
familyforge \
  --vcf        /path/to/annotated.vcf.gz \
  --ped        /path/to/family.ped \
  --panel      /path/to/gene_panel.tsv \
  --cnv        /path/to/cnv_calls.tsv \
  --mz-twins   TWIN_1,TWIN_2 \
  --study-name "Family FAM01 exome" \
  --outdir     /path/to/results/reports
```

Then render the letterhead clinical report from the machine-readable output:

```bash
familyforge-clinical-report \
  --lab-config    /path/to/lab_branding.yaml \
  --report-config /path/to/report_meta.yaml \
  --candidates    /path/to/results/reports/candidates.json \
  --ped           /path/to/family.ped \
  --out           /path/to/results/reports/clinical_report.html
# render to PDF: chrome --headless=new --no-pdf-header-footer \
#   --print-to-pdf=clinical_report.pdf file:///.../clinical_report.html
```

Verify your install with the bundled self-test (plants a STAT3 de-novo, a DOCK8
homozygous stop, and a PGM3 compound-het; expects all three at Tier A):

```bash
python tests/selftest.py     # -> SELFTEST PASS
```

## Commands

### `familyforge` — run the full interpretation

| Flag | Req. | Purpose |
|---|---|---|
| `--vcf` | yes | Annotated joint VCF.gz (CallForge stage-11 output) |
| `--ped` | yes | 6-column PED pedigree |
| `--panel` | yes | Gene-panel TSV (`gene,tier,inheritance,mechanism,omim,notes`) |
| `--outdir` | yes | Output directory |
| `--cnv` | | CallForge `cnv_calls.tsv` for the CNV cross-check |
| `--config` | | Override YAML (merged over `config/default.yaml`) |
| `--study-name` | | Report title (default `Family study`) |
| `--mz-twins` | | Comma-separated MZ-twin sample IDs; repeatable |
| `--provenance` | | CallForge `provenance.json` for the methods section |
| `--extra-scores` | | TSV of finalist CADD/REVEL/AlphaMissense scores |
| `--skip-norm` | | Skip `bcftools norm` (VCF already biallelic-split) |

### `familyforge-clinical-report` — render the letterhead clinical report

| Flag | Purpose |
|---|---|
| `--lab-config` | Reusable lab-branding YAML (institute, logos, head-of-lab, confidentiality) |
| `--report-config` | Per-report YAML (accession, specimen, disease name, sample labels, QC rows, management) |
| `--candidates` | `candidates.json` from a `familyforge` run |
| `--ped` | The pedigree (sex/relationship fallback) |
| `--out` | Output HTML path (a `.md` twin is written alongside) |

## Inputs

- **VCF** — a biallelic-friendly, VEP+ClinVar+gnomAD+PhyloP annotated multi-sample VCF
  (FamilyForge normalises with `bcftools norm` when available).
- **Pedigree** — a 6-column PED. Unsequenced members (e.g. an absent father) are kept for
  relationships; a member is "sequenced" iff its ID appears in the VCF.
- **Panel** — a TSV with header `gene tier inheritance mechanism omim notes`.
- **Config** — thresholds and gene rules in `familyforge/config/default.yaml`; override with
  `--config`.

## Outputs

- `clinician_report.{html,md}` — concise, actionable priority findings + QC + caveats.
- `publication_report.{html,md}` — full methods, tiered tables (A/B/C), per-candidate
  evidence, provisional ACMG-code table, QC, panel, limitations.
- `candidates.json` — machine-readable ranked candidates + QC (input to the letterhead engine).
- `clinical_report.{html,md}` (+ PDF via headless Chrome) — institutional letterhead report.

## Documentation

A full user manual is in [`docs/manual/`](docs/manual/):
[`familyforge_manual.pdf`](docs/manual/familyforge_manual.pdf) /
[`.docx`](docs/manual/familyforge_manual.docx) (source
[`familyforge_manual.md`](docs/manual/familyforge_manual.md)), with an architecture
schematic and a tiering-decision diagram. It documents concepts, every input/output,
both commands, the configuration, the full interpretation logic, the report templates,
troubleshooting, and reproducibility.

## Limitations

Candidates are research outputs pending laboratory ACMG classification and orthogonal
confirmation (Sanger / MLPA / aCGH). Interpretation quality depends on the upstream
annotation (population-AF completeness, in-silico coverage). A missing parent makes some
phase/de-novo calls inferred, not confirmed. The primary pass is panel-restricted. The
letterhead engine is disease-agnostic; the clinician/publication templates currently carry
HIES-oriented narrative framing (the analysis and data are fully general) — see manual §11.

## Citation

See [`CITATION.cff`](CITATION.cff). Author: **Bright Adu**, Genomics & Bioinformatics Core,
Noguchi Memorial Institute for Medical Research, University of Ghana, Legon, Ghana.

## License

MIT — see [`LICENSE`](LICENSE).
