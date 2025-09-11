"""RPT - Report Generation Tool for SV analysis."""

__version__ = "0.1.0"
__author__ = "Minghao Jiang"

import polars as pl

SVS_KNOWN = ("KMT2A::AFF1",)
SIZE_GENOME = 3000  # Mb
TYPES_SV = ["DEL", "DUP", "INS", "INV", "BND"]
SCHEMA_OVERRIDES = {
    "AnnotSV.GC_content_left": pl.String,
    "AnnotSV.GC_content_right": pl.String,
    "AnnotSV.SV_chrom": pl.String,
}
NULLS = ["", "NA", r"\."]
