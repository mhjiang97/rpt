# RPT - Report Generation Tool

![License GPLv3](https://img.shields.io/badge/License-GPLv3-blue.svg)
![Python >= 3.11](https://img.shields.io/badge/Python-%3E%3D3.11-blue.svg)

RPT turns the structural variant (SV) output of the [smk_sv](https://github.com/jasonwong-lab/smk_sv.git) long-read pipeline into a patient-level PDF report.

For each sample, the report includes:

- **Patient and specimen details**, taken from a sample table
- **Results summary**, with the known gene rearrangement (if any) and the tumor SV burden
- **Biomarker findings**: an explanation of tumor SV burden
- **Genomic findings**: background on gene rearrangements and a description of the detected fusion
- **Therapies** relevant to the detected rearrangement (only shown when one is found)

## Prerequisites

- Python >= 3.11
- A [ReportLab](https://www.reportlab.com/) account, needed to install `rlextra` (ReportLab PLUS) from ReportLab's private package index
- smk_sv output for the sample (see [Input layout](#input-layout))

## Installation

1. Install RPT with [pipx](https://pipx.pypa.io/):

   ```bash
   git clone https://github.com/mhjiang97/rpt.git
   cd rpt

   pipx install .
   ```

2. Add `rlextra` to the same environment (you will be asked for your ReportLab credentials):

   ```bash
   pipx inject rpt "rlextra[pycairo]" --index-url https://www.reportlab.com/pypi/
   ```

Without a ReportLab PLUS licence, every page carries an "Evaluation copy of ReportLab PLUS" line at the bottom.

## Usage

```bash
rpt \
    --out-dir <output_directory> \
    --survivor-dir <smk_sv_run>/survivor/<sample> \
    --callers <caller_a> --callers <caller_b> ... \
    --sample-table <sample_table.csv>
```

The report is written to `<output_directory>/<sample>.pdf`. The sample name is taken from the last component of `--survivor-dir`.

### Options

| Option           | Description                                                                      | Default                                            |
| ---------------- | -------------------------------------------------------------------------------- | -------------------------------------------------- |
| `--out-dir`      | Directory for the PDF report (created if missing)                                | `tests/output`                                     |
| `--survivor-dir` | SURVIVOR output directory of one sample, *i.e.* `<smk_sv_run>/survivor/<sample>` | `tests/data/survivor/GYB`                          |
| `--callers`      | SV callers used in the smk_sv run; repeat the option once per caller             | `cutesv`, `sniffles`, `svim`, `severus`, `svision` |
| `--sample-table` | CSV file with patient and specimen details                                       | `tests/data/sample_info.csv`                       |

Run `rpt --help` to see all options.

### Input layout

RPT expects the directory structure produced by smk_sv:

```text
<smk_sv_run>/
├── survivor/<sample>/final/
│   ├── <sample>.BND.merged.vcf
│   ├── <sample>.DEL.merged.vcf
│   ├── <sample>.DUP.merged.vcf
│   ├── <sample>.INS.merged.vcf
│   └── <sample>.INV.merged.vcf
└── <caller>/<sample>/merged/filtered/<caller>.BND.tsv   # one per caller, with SnpEff annotation
```

- **Merged VCFs** — all five are used for the SV burden; the BND VCF supplies the translocation IDs.
- **Per-caller BND TSVs** — searched for SnpEff `gene_fusion` effects. `svision` is accepted in `--callers` but its TSV is not used for annotation.

### Sample table format

A CSV file with one row per sample and these columns:

| Column            | Shown in the report as                                        |
| ----------------- | ------------------------------------------------------------- |
| `sample_name`     | Patient name (must match the `--survivor-dir` directory name) |
| `diagnosis`       | Diagnosis / tumor type                                        |
| `birth_date`      | Date of birth                                                 |
| `gender`          | Sex                                                           |
| `specimen`        | Specimen site                                                 |
| `collection_date` | Date of collection                                            |

Example:

```csv
sample_name,diagnosis,birth_date,gender,specimen,collection_date
GYB,B-cell leukemia,01-Jan-1970,Male,Bone marrow,01-Jan-2022
```

Values are printed as given, so use whatever date format you want to appear in the report.

## How results are computed

- **Known gene rearrangement** — a translocation counts as a known rearrangement when its ID appears in the merged BND VCF and a caller's SnpEff annotation reports a `gene_fusion` covering both genes of a curated fusion (*e.g.* `KMT2A::AFF1`).
- **Tumor SV burden** — the total number of merged SV records (DEL, DUP, INS, INV and BND) divided by a genome size of 3,000 Mb, reported as variants/Mb.

## Customization

### Adding a known rearrangement

1. Add the fusion to `SVS_KNOWN` in [`src/rpt/__init__.py`](src/rpt/__init__.py), using `GENE1::GENE2` notation.
2. Add a matching entry to [`src/rpt/templates/descriptions.json`](src/rpt/templates/descriptions.json):

   ```json
   {
       "GENE1::GENE2": {
           "description": ["First paragraph.", "Second paragraph."],
           "therapies": [
               {"name": "Therapy class", "text": "Description of the therapy."}
           ]
       }
   }
   ```

   Each `description` item becomes a paragraph in the report. If `therapies` is empty, the Therapies section is left out.

### Changing the report layout

The report is built from the ReportLab RML template [`src/rpt/templates/rpt.prep`](src/rpt/templates/rpt.prep) using [preppy](https://pypi.org/project/preppy/). Colors are defined once in its `<docinit>` block and text styles in its `<stylesheet>`. See the [RML reference](tests/examples/rml-reference.pdf) for the available elements.

## Development

```bash
uv sync
uv pip install "rlextra[pycairo]" --index-url https://www.reportlab.com/pypi/

# Render the example report from the bundled test data
PYTHONPATH=src uv run python -m rpt --out-dir tests/output
```

This writes `tests/output/GYB.pdf` from the data in `tests/data`.

## License

RPT is licensed under the [GNU General Public License v3.0 or later](LICENSE).
