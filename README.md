# RPT - Report Generation Tool

![License GPLv3](https://img.shields.io/badge/License-GPLv3-blue.svg)

RPT is a command-line tool for generating PDF reports from SV VCF files. It processes SURVIVOR output and sample information to create comprehensive reports.

## Prerequisites

- Python >= 3.11
- A ReportLab account (required for rlextra installation)
- SURVIVOR output files (from [smk_sv](https://github.com/jasonwong-lab/smk_sv.git))

## Installation

1. Install the package using pipx:

   ```bash
   git clone https://github.com/mhjiang97/rpt.git
   cd rpt

   pipx install .
   ```

2. Activate the virtual environment:

   ```bash
   source ~/.local/pipx/venvs/rpt/bin/activate
   ```

3. Install rlextra (requires ReportLab account):

   ```bash
   pip install "rlextra[pycairo]" -i https://www.reportlab.com/pypi/
   ```

## Usage

After running [smk_sv](https://github.com/jasonwong-lab/smk_sv.git), generate a report with:

```bash
rpt \
    --out-dir <output_directory> \
    --survivor-dir <survivor_directory> \
    --callers <caller_a> --callers <caller_b> ... \
    --sample-table <sample_table.csv>
```

### Options

- `--out-dir`: Directory to save the generated report
- `--survivor-dir`: Path to SURVIVOR output directory containing `final/` subdirectory with merged VCF files and annotated TSV files
- `--callers`: List of SV callers used in the smk_sv pipeline (*e.g.*, `cutesv`, `sniffles`, `svim`, `severus`, `svision`)
- `--sample-table`: CSV file with sample information

### Sample Table Format

The CSV file must include the following columns:

- `sample_name`
- `diagnosis`
- `birth_date`
- `gender`
- `specimen`
- `collection_date`
