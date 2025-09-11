"""Generate a PDF report from a VCF."""

# %%
import json
import subprocess as sp
from datetime import datetime
from importlib.resources import files
from pathlib import Path

import polars as pl
import preppy
import typer
from rlextra.rml2pdf import rml2pdf
from typing_extensions import Annotated

from . import SIZE_GENOME, SVS_KNOWN


# %%
def find_known_rearrangement(snpeff_bnd):
    genes = (
        snpeff_bnd.filter(pl.col("EFFECT").str.contains("gene_fusion"))
        .select(["GENE"])
        .to_series()
        .unique()
        .str.split("&")
        .list.sort()
    )

    for lst in genes:
        for sv_known in SVS_KNOWN:
            if all(g in lst for g in sv_known.split("::")):
                return sv_known

    return "No known gene rearrangement detected"


def calculate_sv_burden(survivor_dir, patient_name):
    types_sv = ["DEL", "DUP", "INS", "INV", "BND"]
    vcfs = {
        sv: f"{survivor_dir}/final/{patient_name}.{sv}.merged.vcf" for sv in types_sv
    }
    cmds = {sv: f"grep -v '^#' {vcf} | wc -l" for sv, vcf in vcfs.items()}
    counts = {
        sv: int(sp.run(cmd, shell=True, capture_output=True, text=True).stdout.strip())
        for sv, cmd in cmds.items()
    }
    count_total = sum(counts.values())
    burden_sv = round(count_total / SIZE_GENOME, 2)

    return burden_sv


# %%
arg_out_dir = Annotated[
    Path, typer.Option(dir_okay=True, resolve_path=True, help="Output directory")
]
arg_survivor_dir = Annotated[
    Path, typer.Option(exists=True, help="Path to SURVIVOR directory")
]
arg_sample_table = Annotated[
    Path, typer.Option(exists=True, help="Path to sample table")
]


# %%
def generate_report(
    out_dir: arg_out_dir = "tests/output",
    survivor_dir: arg_survivor_dir = "tests/data/survivor/GYB",
    sample_table: arg_sample_table = "tests/data/sample_info.csv",
):
    out_dir.mkdir(parents=True, exist_ok=True)
    name_patient = survivor_dir.name
    path_snpeff_bnd = f"{survivor_dir}/final/{name_patient}.BND.snpeff.tsv"

    df_sample = pl.read_csv(sample_table)
    index_patient = df_sample["sample_name"].to_list().index(name_patient)

    diagnosis = df_sample["diagnosis"][index_patient]
    date_birth = df_sample["birth_date"][index_patient]
    sex = df_sample["gender"][index_patient]
    site_specimen = df_sample["specimen"][index_patient]
    date_collection = df_sample["collection_date"][index_patient]
    date_report = datetime.now().strftime("%d-%b-%Y")

    path_out = Path(f"{out_dir}/{name_patient}.pdf")

    path_template = files("rpt.templates").joinpath("rpt.prep")
    path_descriptions = files("rpt.templates").joinpath("descriptions.json")

    template = preppy.getModule(path_template)
    with open(path_descriptions) as f:
        descriptions = json.load(f)

    df_snpeff_bnd = pl.read_csv(path_snpeff_bnd, separator="\t")
    rearrangements_known = find_known_rearrangement(df_snpeff_bnd)

    if rearrangements_known == "No known gene rearrangement detected":
        description_rearrangement = "N/A"
    else:
        description_rearrangement = descriptions[rearrangements_known]

    burden_sv = calculate_sv_burden(survivor_dir, name_patient)

    text_rml = template.get(
        name_patient,
        sex,
        date_birth,
        diagnosis,
        site_specimen,
        date_collection,
        date_report,
        rearrangements_known,
        burden_sv,
        description_rearrangement,
    )

    rml2pdf.go(text_rml, outputFileName=path_out.as_posix())


# %%
