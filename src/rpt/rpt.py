"""Generate a PDF report from a VCF."""

# %%
import json
import subprocess as sp
from datetime import datetime
from importlib.resources import files
from pathlib import Path

import oxbow as ox
import polars as pl
import polars.selectors as cs
import preppy
import typer
from rlextra.rml2pdf import rml2pdf
from typing_extensions import Annotated

from . import NULLS, SCHEMA_OVERRIDES, SIZE_GENOME, SVS_KNOWN, TYPES_SV


# %%
def find_known_rearrangement(snpeff_bnd, translocation_IDS):
    genes = (
        snpeff_bnd.filter(pl.col("ID").is_in(translocation_IDS))
        .filter(pl.col("EFFECT").str.contains("gene_fusion"))
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
    vcfs = {
        sv: f"{survivor_dir}/final/{patient_name}.{sv}.merged.vcf" for sv in TYPES_SV
    }
    cmds = {sv: f"grep -v '^#' {vcf} | wc -l" for sv, vcf in vcfs.items()}
    counts = {
        sv: int(sp.run(cmd, shell=True, capture_output=True, text=True).stdout.strip())
        for sv, cmd in cmds.items()
    }
    count_total = sum(counts.values())
    burden_sv = round(count_total / SIZE_GENOME, 2)

    return burden_sv


def get_translocation_IDs(BND_vcf):
    df_bnd = ox.from_vcf(BND_vcf).pl()
    ids = [_ids[0] for _ids in df_bnd["id"].to_list()]

    return ids


def extract_annotation(caller_names, pipeline_work_dir, patient_name):
    dfs = []
    for caller in caller_names:
        if caller == "svision":
            continue

        path_tsv = f"{pipeline_work_dir}/{caller}/{patient_name}/merged/filtered/{caller}.BND.tsv"
        df = (
            pl.read_csv(
                path_tsv,
                separator="\t",
                null_values=NULLS,
                schema_overrides=SCHEMA_OVERRIDES,
            )
            .select(cs.starts_with("SnpEff"))
            .rename(lambda col_name: col_name.replace("SnpEff.", ""))
            .filter(pl.col("ID").is_not_null())
            .select(["ID", "EFFECT", "GENE"])
            .with_columns(pl.lit(f"{caller}").alias("Caller"))
        )
        dfs.append(df)

    df_merged = pl.concat(dfs, how="vertical_relaxed")

    return df_merged


def get_description(rearrangement):
    path_descriptions = files("rpt.templates").joinpath("descriptions.json")

    with open(path_descriptions) as f:
        descriptions = json.load(f)

    if rearrangement in descriptions:
        return descriptions[rearrangement]
    else:
        return "N/A"


def create_pdf_report(
    sample_table,
    patient_name,
    sv_burden,
    known_rearrangement,
    output_pdf,
):
    df_sample = pl.read_csv(sample_table)
    index_patient = df_sample["sample_name"].to_list().index(patient_name)
    diagnosis = df_sample["diagnosis"][index_patient]
    date_birth = df_sample["birth_date"][index_patient]
    sex = df_sample["gender"][index_patient]
    site_specimen = df_sample["specimen"][index_patient]
    date_collection = df_sample["collection_date"][index_patient]
    date_report = datetime.now().strftime("%d-%b-%Y")

    description_rearrangement = get_description(known_rearrangement)

    path_template = files("rpt.templates").joinpath("rpt.prep")
    template = preppy.getModule(path_template)
    text_rml = template.get(
        patient_name,
        sex,
        date_birth,
        diagnosis,
        site_specimen,
        date_collection,
        date_report,
        known_rearrangement,
        sv_burden,
        description_rearrangement,
    )

    rml2pdf.go(text_rml, outputFileName=output_pdf)


# %%
arg_out_dir = Annotated[
    Path, typer.Option(dir_okay=True, resolve_path=True, help="Output directory")
]
arg_survivor_dir = Annotated[
    Path, typer.Option(exists=True, help="Path to SURVIVOR directory")
]
arg_callers = Annotated[
    list[str], typer.Option(help="SV callers used in the smk_sv pipeline")
]
arg_sample_table = Annotated[
    Path, typer.Option(exists=True, help="Path to sample table")
]


# %%
def main(
    out_dir: arg_out_dir = "tests/output",
    survivor_dir: arg_survivor_dir = "tests/data/survivor/GYB",
    callers: arg_callers = ["cutesv", "sniffles", "svim", "severus", "svision"],
    sample_table: arg_sample_table = "tests/data/sample_info.csv",
):
    out_dir.mkdir(parents=True, exist_ok=True)
    name_patient = survivor_dir.name
    path_vcf_bnd = f"{survivor_dir}/final/{name_patient}.BND.merged.vcf"
    path_run = survivor_dir.parent.parent
    path_out = Path(f"{out_dir}/{name_patient}.pdf")

    ids = get_translocation_IDs(path_vcf_bnd)
    df_bnd = extract_annotation(callers, path_run, name_patient)
    rearrangement = find_known_rearrangement(df_bnd, ids)
    burden = calculate_sv_burden(survivor_dir, name_patient)

    create_pdf_report(
        sample_table=sample_table,
        patient_name=name_patient,
        sv_burden=burden,
        known_rearrangement=rearrangement,
        output_pdf=path_out.as_posix(),
    )


# %%
