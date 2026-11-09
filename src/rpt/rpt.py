"""Generate a PDF report from a VCF."""

# %%
import json
from datetime import datetime
from importlib.resources import files
from pathlib import Path
from typing import Annotated
from xml.sax.saxutils import escape

import oxbow as ox
import polars as pl
import polars.selectors as cs
import preppy
import typer
from rlextra.rml2pdf import rml2pdf

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

    return None


def calculate_sv_burden(survivor_dir, patient_name):
    vcfs = {
        sv: f"{survivor_dir}/final/{patient_name}.{sv}.merged.vcf" for sv in TYPES_SV
    }
    count_total = 0
    for vcf in vcfs.values():
        with open(vcf) as f:
            count_total += sum(1 for line in f if not line.startswith("#"))
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
    """Return the description paragraphs and therapies for a rearrangement."""
    path_descriptions = files("rpt.templates").joinpath("descriptions.json")

    with open(path_descriptions) as f:
        descriptions = json.load(f)

    return descriptions.get(rearrangement, {"description": [], "therapies": []})


def _quote_xml(value):
    return escape(str(value))


def create_pdf_report(
    sample_table,
    patient_name,
    sv_burden,
    known_rearrangement,
    output_pdf,
):
    df_sample = pl.read_csv(sample_table, infer_schema=False).filter(
        pl.col("sample_name") == patient_name
    )
    if df_sample.is_empty():
        raise ValueError(f"Sample {patient_name!r} not found in {sample_table}")
    sample = df_sample.row(0, named=True)

    knowledge = get_description(known_rearrangement)

    path_template = files("rpt.templates").joinpath("rpt.prep")
    template = preppy.getModule(str(path_template), savePyc=False)
    text_rml = template.get(
        patient_name=patient_name,
        sex=sample["gender"],
        birth_date=sample["birth_date"],
        tumor_type=sample["diagnosis"],
        specimen_site=sample["specimen"],
        collection_date=sample["collection_date"],
        report_date=datetime.now().strftime("%d-%b-%Y"),
        known_rearrangement=known_rearrangement,
        sv_burden=sv_burden,
        rearrangement_description=knowledge["description"],
        therapies=knowledge["therapies"],
        __quoteFunc__=_quote_xml,
    )

    rml2pdf.go(text_rml, outputFileName=output_pdf)


# %%
arg_out_dir = Annotated[
    Path, typer.Option(dir_okay=True, resolve_path=True, help="Output directory")
]
arg_survivor_dir = Annotated[
    Path, typer.Option(exists=True, help="Path to SURVIVOR directory")
]
DEFAULT_CALLERS = ("cutesv", "sniffles", "svim", "severus", "svision")
arg_callers = Annotated[
    list[str] | None,
    typer.Option(
        help="SV callers used in the smk_sv pipeline",
        show_default=", ".join(DEFAULT_CALLERS),
    ),
]
arg_sample_table = Annotated[
    Path, typer.Option(exists=True, help="Path to sample table")
]


# %%
def main(
    out_dir: arg_out_dir = "tests/output",
    survivor_dir: arg_survivor_dir = "tests/data/survivor/GYB",
    callers: arg_callers = None,
    sample_table: arg_sample_table = "tests/data/sample_info.csv",
):
    if not callers:
        callers = list(DEFAULT_CALLERS)
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
