#!/usr/bin/env python3
"""Build compact supplementary tables for performance-variation analyses."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Any
import xml.etree.ElementTree as ET
import zipfile

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
SIM_SCRIPT_DIR = SCRIPT_DIR.parent
PLOT_ROOT = SIM_SCRIPT_DIR.parents[1]
if str(PLOT_ROOT) not in sys.path:
    sys.path.insert(0, str(PLOT_ROOT))

from _shared.io import ensure_dir, write_tsv  # noqa: E402

OUTPUT_DIR = SCRIPT_DIR / "output"
PRIMARY_TRANSFORM = "logit_clip"

TABLE_PROTOCOL = "supplementary_table_S8A_protocol_associated_variation.tsv"
TABLE_DESCRIPTOR = "supplementary_table_S8B_peak_morphology_associations.tsv"
TABLE_SENSITIVITY = "supplementary_table_S8C_sensitivity_diagnostics.tsv"
WORKBOOK = "Supplementary_Table_S8_protocol_peak_morphology_performance_variation.xlsx"

SHEETS = [
    ("S8-notes", "notes"),
    ("S8A Protocol variation", TABLE_PROTOCOL),
    ("S8B Morphology assoc", TABLE_DESCRIPTOR),
    ("S8C Robustness diag", TABLE_SENSITIVITY),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--response-transform", default=PRIMARY_TRANSFORM)
    return parser.parse_args()


def read_tsv(input_dir: Path, name: str) -> pd.DataFrame:
    path = input_dir / name
    if not path.exists():
        raise FileNotFoundError(f"Missing required input: {path}")
    return pd.read_csv(path, sep="\t")


def excel_col_name(index: int) -> str:
    name = ""
    index += 1
    while index:
        index, rem = divmod(index - 1, 26)
        name = chr(65 + rem) + name
    return name


def xml_text_cell(row_idx: int, col_idx: int, value: Any) -> ET.Element:
    cell = ET.Element("c", {"r": f"{excel_col_name(col_idx)}{row_idx}", "t": "inlineStr"})
    inline = ET.SubElement(cell, "is")
    text = ET.SubElement(inline, "t")
    if isinstance(value, str) and (value.startswith(" ") or value.endswith(" ")):
        text.set("{http://www.w3.org/XML/1998/namespace}space", "preserve")
    text.text = "" if pd.isna(value) else str(value)
    return cell


def xml_numeric_cell(row_idx: int, col_idx: int, value: Any) -> ET.Element:
    cell = ET.Element("c", {"r": f"{excel_col_name(col_idx)}{row_idx}"})
    val = ET.SubElement(cell, "v")
    val.text = str(value)
    return cell


def worksheet_xml(df: pd.DataFrame) -> bytes:
    ns = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
    ET.register_namespace("", ns)
    sheet = ET.Element(f"{{{ns}}}worksheet")
    dimension = ET.SubElement(sheet, f"{{{ns}}}dimension")
    last_col = excel_col_name(max(0, len(df.columns) - 1))
    dimension.set("ref", f"A1:{last_col}{len(df) + 1}")
    sheet_data = ET.SubElement(sheet, f"{{{ns}}}sheetData")

    header = ET.SubElement(sheet_data, f"{{{ns}}}row", {"r": "1"})
    for col_idx, col in enumerate(df.columns):
        header.append(xml_text_cell(1, col_idx, col))

    for row_offset, (_, row) in enumerate(df.iterrows(), start=2):
        xml_row = ET.SubElement(sheet_data, f"{{{ns}}}row", {"r": str(row_offset)})
        for col_idx, value in enumerate(row.tolist()):
            if pd.isna(value):
                continue
            if isinstance(value, (int, float, np.integer, np.floating)) and np.isfinite(value):
                xml_row.append(xml_numeric_cell(row_offset, col_idx, value))
            else:
                xml_row.append(xml_text_cell(row_offset, col_idx, value))
    return ET.tostring(sheet, encoding="utf-8", xml_declaration=True)


def write_xlsx(path: Path, sheets: list[tuple[str, pd.DataFrame]]) -> None:
    content_types = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>
""" + "".join(
        f'  <Override PartName="/xl/worksheets/sheet{i}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>\n'
        for i in range(1, len(sheets) + 1)
    ) + "</Types>\n"
    root_rels = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
</Relationships>
"""
    workbook_xml = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
  <sheets>
""" + "".join(
        f'    <sheet name="{name}" sheetId="{i}" r:id="rId{i}"/>\n'
        for i, (name, _) in enumerate(sheets, start=1)
    ) + "  </sheets>\n</workbook>\n"
    workbook_rels = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
""" + "".join(
        f'  <Relationship Id="rId{i}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i}.xml"/>\n'
        for i in range(1, len(sheets) + 1)
    ) + "</Relationships>\n"

    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", content_types)
        zf.writestr("_rels/.rels", root_rels)
        zf.writestr("xl/workbook.xml", workbook_xml)
        zf.writestr("xl/_rels/workbook.xml.rels", workbook_rels)
        for idx, (_, df) in enumerate(sheets, start=1):
            zf.writestr(f"xl/worksheets/sheet{idx}.xml", worksheet_xml(df))


def primary_rows(df: pd.DataFrame, transform: str) -> pd.DataFrame:
    out = df[df["response_transform"].eq(transform)].copy()
    if "is_primary" in out.columns:
        out = out[out["is_primary"].astype(str).str.lower().isin(["true", "1"])]
    return out


def compact_flag(flags: list[Any]) -> str:
    clean = [str(flag) for flag in flags if pd.notna(flag) and str(flag) != ""]
    if not clean:
        return ""
    non_ok = [flag for flag in clean if flag != "ok"]
    return ";".join(sorted(set(non_ok or clean)))


def build_notes_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "item": "caption",
                "note": (
                    "Supplementary Table S8. Protocol- and peak-morphology-associated "
                    "variation in workflow performance."
                ),
            },
            {
                "item": "sheet_A",
                "note": (
                    "Sheet A summarizes protocol-associated variation in workflow "
                    "performance using instance-adjusted relative DE-APA performance."
                ),
            },
            {
                "item": "sheet_B",
                "note": (
                    "Sheet B summarizes workflow-specific associations between relative "
                    "DE-APA performance and individual peak descriptors."
                ),
            },
            {
                "item": "sheet_C",
                "note": (
                    "Sheet C summarizes descriptor/protocol variance-decomposition analyses "
                    "and robustness checks, including sample-collapsed, complete-case, "
                    "response-transform, and study leave-one-out sensitivity analyses."
                ),
            },
            {
                "item": "primary_inference",
                "note": (
                    "Primary P value and FDR columns in sheets A and B come from "
                    "sample-collapsed models fitted after averaging relative performance "
                    "within empirical_sample_id x workflow."
                ),
            },
            {
                "item": "study_proxy",
                "note": (
                    "Study groups are derived from the available source grouping field for leave-one-study-out "
                    "sensitivity because an explicit study identifier is not available in "
                    "the joined performance tables."
                ),
            },
            {
                "item": "r2_scale",
                "note": (
                    "R2 values were calculated on the clipped-logit instance-adjusted "
                    "relative-performance scale and should not be interpreted as variance "
                    "explained in raw F1, precision, or recall."
                ),
            },
            {
                "item": "descriptor_interpretation",
                "note": (
                    "Slopes were estimated on the clipped-logit instance-adjusted "
                    "relative-performance scale. Positive slopes indicate that a workflow "
                    "performed better relative to other workflows as the corresponding "
                    "descriptor increased. Descriptor rows should be interpreted as "
                    "individual or standalone associations, not independent contributions."
                ),
            },
        ]
    )


def rename_for_supplement(df: pd.DataFrame, rename_map: dict[str, str]) -> pd.DataFrame:
    """Rename compact-table columns for paper-facing supplementary outputs."""
    return df.rename(columns={old: new for old, new in rename_map.items() if old in df.columns})


def build_protocol_table(input_dir: Path, transform: str) -> pd.DataFrame:
    protocol = primary_rows(read_tsv(input_dir, "protocol_response_by_workflow.tsv"), transform)
    sample_protocol = primary_rows(read_tsv(input_dir, "sample_collapsed_protocol_response_by_workflow.tsv"), transform)
    estimates = primary_rows(read_tsv(input_dir, "protocol_response_estimates.tsv"), transform)

    idx = ["regime", "metric", "workflow"]
    high = (
        estimates.sort_values(idx + ["mean_relative_response"], ascending=[True, True, True, False])
        .groupby(idx, as_index=False)
        .first()
        .rename(
            columns={
                "protocol": "protocol_with_highest_relative_response",
                "mean_relative_response": "highest_mean_relative_response",
                "standard_error": "highest_protocol_standard_error",
                "n_instances": "highest_protocol_n_instances",
            }
        )
    )
    low = (
        estimates.sort_values(idx + ["mean_relative_response"], ascending=[True, True, True, True])
        .groupby(idx, as_index=False)
        .first()
        .rename(
            columns={
                "protocol": "protocol_with_lowest_relative_response",
                "mean_relative_response": "lowest_mean_relative_response",
                "standard_error": "lowest_protocol_standard_error",
                "n_instances": "lowest_protocol_n_instances",
            }
        )
    )
    keep_high = idx + [
        "protocol_with_highest_relative_response",
        "highest_mean_relative_response",
        "highest_protocol_standard_error",
        "highest_protocol_n_instances",
    ]
    keep_low = idx + [
        "protocol_with_lowest_relative_response",
        "lowest_mean_relative_response",
        "lowest_protocol_standard_error",
        "lowest_protocol_n_instances",
    ]
    out = protocol.merge(high[keep_high], on=idx, how="left").merge(low[keep_low], on=idx, how="left")
    out = out.rename(
        columns={
            "r2": "R2_protocol",
            "partial_r2": "partial_R2_protocol",
            "f_stat": "instance_level_F_stat",
            "p_value": "instance_level_p_value",
            "fdr_bh": "instance_level_FDR",
            "diagnostic_flag": "instance_level_diagnostic_flag",
        }
    )
    sample_keep = idx + ["n_instances", "r2", "partial_r2", "f_stat", "p_value", "fdr_bh", "diagnostic_flag"]
    sample_protocol = sample_protocol[sample_keep].rename(
        columns={
            "n_instances": "sample_collapsed_n_samples",
            "r2": "sample_collapsed_R2_protocol",
            "partial_r2": "sample_collapsed_partial_R2_protocol",
            "f_stat": "sample_collapsed_F_stat",
            "p_value": "p_value",
            "fdr_bh": "FDR",
            "diagnostic_flag": "sample_collapsed_diagnostic_flag",
        }
    )
    out = out.merge(sample_protocol, on=idx, how="left")
    columns = [
        "regime",
        "metric",
        "workflow",
        "n_instances",
        "n_protocols",
        "R2_protocol",
        "partial_R2_protocol",
        "sample_collapsed_n_samples",
        "sample_collapsed_R2_protocol",
        "sample_collapsed_partial_R2_protocol",
        "sample_collapsed_F_stat",
        "p_value",
        "FDR",
        "instance_level_F_stat",
        "instance_level_p_value",
        "instance_level_FDR",
        "sample_collapsed_diagnostic_flag",
        "instance_level_diagnostic_flag",
        "protocol_with_highest_relative_response",
        "highest_mean_relative_response",
        "highest_protocol_standard_error",
        "highest_protocol_n_instances",
        "protocol_with_lowest_relative_response",
        "lowest_mean_relative_response",
        "lowest_protocol_standard_error",
        "lowest_protocol_n_instances",
    ]
    out = out[columns].sort_values(["regime", "metric", "workflow"]).reset_index(drop=True)
    return rename_for_supplement(
        out,
        {
            "regime": "Regime",
            "metric": "Metric",
            "workflow": "Workflow",
            "n_instances": "Instance-level N instances",
            "n_protocols": "N protocols",
            "R2_protocol": "Instance-level R2",
            "partial_R2_protocol": "Instance-level partial R2",
            "sample_collapsed_n_samples": "Sample-level N samples",
            "sample_collapsed_R2_protocol": "Sample-level R2",
            "sample_collapsed_partial_R2_protocol": "Sample-level partial R2",
            "sample_collapsed_F_stat": "Sample-level F statistic",
            "p_value": "Sample-level P value",
            "FDR": "Sample-level FDR",
            "instance_level_F_stat": "Instance-level F statistic",
            "instance_level_p_value": "Instance-level P value",
            "instance_level_FDR": "Instance-level FDR",
            "sample_collapsed_diagnostic_flag": "Sample-level diagnostic flag",
            "instance_level_diagnostic_flag": "Instance-level diagnostic flag",
            "protocol_with_highest_relative_response": "Protocol with highest relative performance",
            "highest_mean_relative_response": "Highest mean relative performance",
            "highest_protocol_standard_error": "Highest protocol standard error",
            "highest_protocol_n_instances": "Highest protocol N instances",
            "protocol_with_lowest_relative_response": "Protocol with lowest relative performance",
            "lowest_mean_relative_response": "Lowest mean relative performance",
            "lowest_protocol_standard_error": "Lowest protocol standard error",
            "lowest_protocol_n_instances": "Lowest protocol N instances",
        },
    )


def build_descriptor_table(input_dir: Path, transform: str) -> pd.DataFrame:
    descriptor = primary_rows(read_tsv(input_dir, "descriptor_response_by_workflow.tsv"), transform)
    sample_descriptor = primary_rows(read_tsv(input_dir, "sample_collapsed_descriptor_response_by_workflow.tsv"), transform)
    idx = ["regime", "metric", "workflow", "descriptor"]
    out = descriptor.rename(
        columns={
            "n_instances": "instance_level_n_instances",
            "slope": "instance_level_slope",
            "standard_error": "instance_level_standard_error",
            "t_stat": "instance_level_t_stat",
            "p_value": "instance_level_p_value",
            "r2_descriptor": "R2_descriptor",
            "partial_r2": "partial_R2_descriptor",
            "fdr_bh": "instance_level_FDR",
            "diagnostic_flag": "instance_level_diagnostic_flag",
        }
    )
    sample_keep = idx + [
        "n_instances",
        "slope",
        "standard_error",
        "t_stat",
        "p_value",
        "r2_descriptor",
        "partial_r2",
        "fdr_bh",
        "diagnostic_flag",
    ]
    sample_descriptor = sample_descriptor[sample_keep].rename(
        columns={
            "n_instances": "sample_collapsed_n_samples",
            "slope": "slope",
            "standard_error": "standard_error",
            "t_stat": "t_stat",
            "r2_descriptor": "sample_collapsed_R2_descriptor",
            "partial_r2": "sample_collapsed_partial_R2_descriptor",
            "fdr_bh": "FDR",
            "diagnostic_flag": "sample_collapsed_diagnostic_flag",
        }
    )
    out = out.merge(sample_descriptor, on=idx, how="left")
    columns = [
        "regime",
        "metric",
        "workflow",
        "descriptor",
        "sample_collapsed_n_samples",
        "slope",
        "standard_error",
        "t_stat",
        "p_value",
        "FDR",
        "sample_collapsed_R2_descriptor",
        "sample_collapsed_partial_R2_descriptor",
        "instance_level_n_instances",
        "instance_level_slope",
        "instance_level_standard_error",
        "instance_level_t_stat",
        "instance_level_p_value",
        "instance_level_FDR",
        "R2_descriptor",
        "partial_R2_descriptor",
        "sample_collapsed_diagnostic_flag",
        "instance_level_diagnostic_flag",
    ]
    out = out[columns].sort_values(["regime", "metric", "workflow", "descriptor"]).reset_index(drop=True)
    return rename_for_supplement(
        out,
        {
            "regime": "Regime",
            "metric": "Metric",
            "workflow": "Workflow",
            "descriptor": "Descriptor",
            "sample_collapsed_n_samples": "Sample-level N samples",
            "slope": "Sample-level slope",
            "standard_error": "Sample-level standard error",
            "t_stat": "Sample-level t statistic",
            "p_value": "Sample-level P value",
            "FDR": "Sample-level FDR",
            "sample_collapsed_R2_descriptor": "Sample-level R2",
            "sample_collapsed_partial_R2_descriptor": "Sample-level partial R2",
            "instance_level_n_instances": "Instance-level N instances",
            "instance_level_slope": "Instance-level slope",
            "instance_level_standard_error": "Instance-level standard error",
            "instance_level_t_stat": "Instance-level t statistic",
            "instance_level_p_value": "Instance-level P value",
            "instance_level_FDR": "Instance-level FDR",
            "R2_descriptor": "Instance-level R2",
            "partial_R2_descriptor": "Instance-level partial R2",
            "sample_collapsed_diagnostic_flag": "Sample-level diagnostic flag",
            "instance_level_diagnostic_flag": "Instance-level diagnostic flag",
        },
    )


def median_descriptor_r2(df: pd.DataFrame, value_col: str, out_col: str, transform: str) -> pd.DataFrame:
    primary = primary_rows(df, transform)
    return (
        primary.groupby(["regime", "metric", "workflow"], as_index=False)[value_col]
        .median()
        .rename(columns={value_col: out_col})
    )


def sensitivity_consistency(
    primary: pd.DataFrame,
    sensitivity: pd.DataFrame,
    *,
    value_col: str,
    out_prefix: str,
    transform: str,
) -> pd.DataFrame:
    idx = ["regime", "metric", "workflow"]
    base = primary_rows(primary, transform)
    sens = primary_rows(sensitivity, transform)
    if "descriptor" in base.columns and "descriptor" in sens.columns:
        idx = [*idx, "descriptor"]
    merged = base[idx + [value_col, "diagnostic_flag"]].merge(
        sens[idx + [value_col, "diagnostic_flag"]],
        on=idx,
        how="left",
        suffixes=("_primary", "_sensitivity"),
    )
    merged[f"{out_prefix}_R2_spearman_by_workflow"] = np.nan
    summary_rows = []
    for keys, sub in merged.groupby(["regime", "metric", "workflow"]):
        valid = sub[f"{value_col}_primary"].notna() & sub[f"{value_col}_sensitivity"].notna()
        if valid.sum() >= 2:
            corr = sub.loc[valid, f"{value_col}_primary"].corr(
                sub.loc[valid, f"{value_col}_sensitivity"], method="spearman"
            )
        elif valid.sum() == 1:
            corr = 1.0
        else:
            corr = np.nan
        summary_rows.append(
            {
                "regime": keys[0],
                "metric": keys[1],
                "workflow": keys[2],
                f"{out_prefix}_R2_spearman_by_workflow": corr,
                f"{out_prefix}_diagnostic_flag": compact_flag(
                    sub["diagnostic_flag_primary"].tolist() + sub["diagnostic_flag_sensitivity"].tolist()
                ),
            }
        )
    return pd.DataFrame(summary_rows)


def build_sensitivity_table(input_dir: Path, transform: str) -> pd.DataFrame:
    explanation = primary_rows(read_tsv(input_dir, "descriptor_explains_protocol_response.tsv"), transform)
    protocol = read_tsv(input_dir, "protocol_response_by_workflow.tsv")
    descriptor = read_tsv(input_dir, "descriptor_response_by_workflow.tsv")
    sample_protocol = read_tsv(input_dir, "sample_collapsed_protocol_response_by_workflow.tsv")
    sample_descriptor = read_tsv(input_dir, "sample_collapsed_descriptor_response_by_workflow.tsv")
    complete_protocol = read_tsv(input_dir, "complete_case_protocol_response_by_workflow.tsv")
    complete_descriptor = read_tsv(input_dir, "complete_case_descriptor_response_by_workflow.tsv")
    study_protocol = read_tsv(input_dir, "study_loo_protocol_sensitivity_summary.tsv")
    study_descriptor = read_tsv(input_dir, "study_loo_descriptor_sensitivity_summary.tsv")
    coverage = read_tsv(input_dir, "morphology_coverage_diagnostics.tsv")
    consistency = read_tsv(input_dir, "transform_consistency_summary.tsv")

    idx = ["regime", "metric", "workflow"]
    out = explanation.rename(
        columns={
            "r2_protocol": "R2_protocol",
            "r2_all_descriptors": "R2_all_descriptors",
            "partial_r2_protocol_after_descriptors": "partial_R2_protocol_after_descriptors",
            "partial_r2_within_protocol_descriptors": "partial_R2_within_protocol_descriptors",
            "fdr_all_descriptors": "FDR_all_descriptors",
            "fdr_protocol_after_descriptors": "FDR_protocol_after_descriptors",
            "fdr_within_protocol_descriptors": "FDR_within_protocol_descriptors",
            "diagnostic_flag": "primary_diagnostic_flag",
        }
    )

    sample_p = primary_rows(sample_protocol, transform)[idx + ["r2", "fdr_bh", "diagnostic_flag"]].rename(
        columns={
            "r2": "sample_collapsed_protocol_R2",
            "fdr_bh": "sample_collapsed_protocol_FDR",
            "diagnostic_flag": "sample_collapsed_protocol_diagnostic_flag",
        }
    )
    complete_p = primary_rows(complete_protocol, transform)[idx + ["r2", "fdr_bh", "diagnostic_flag"]].rename(
        columns={
            "r2": "complete_case_protocol_R2",
            "fdr_bh": "complete_case_protocol_FDR",
            "diagnostic_flag": "complete_case_protocol_diagnostic_flag",
        }
    )
    sample_desc = median_descriptor_r2(sample_descriptor, "r2_descriptor", "sample_collapsed_descriptor_median_R2", transform)
    complete_desc = median_descriptor_r2(complete_descriptor, "r2_descriptor", "complete_case_descriptor_median_R2", transform)
    primary_desc = median_descriptor_r2(descriptor, "r2_descriptor", "primary_descriptor_median_R2", transform)
    study_p = primary_rows(study_protocol, transform)[
        idx + ["loo_R2_max_abs_delta", "all_loo_significant", "n_studies_tested", "n_studies_ok", "diagnostic_flag"]
    ].rename(
        columns={
            "loo_R2_max_abs_delta": "study_loo_protocol_R2_max_abs_delta",
            "all_loo_significant": "study_loo_protocol_all_significant",
            "n_studies_tested": "study_loo_protocol_n_studies_tested",
            "n_studies_ok": "study_loo_protocol_n_studies_ok",
            "diagnostic_flag": "study_loo_protocol_diagnostic_flag",
        }
    )
    study_d = primary_rows(study_descriptor, transform)
    study_desc = (
        study_d.groupby(idx, as_index=False)
        .agg(
            study_loo_descriptor_median_direction_agreement=(
                "loo_slope_direction_agreement_fraction",
                "median",
            ),
            study_loo_descriptor_median_R2_max_abs_delta=("loo_R2_max_abs_delta", "median"),
            study_loo_descriptor_min_n_studies_ok=("n_studies_ok", "min"),
            study_loo_descriptor_all_significant=("all_loo_significant", "all"),
            study_loo_descriptor_diagnostic_flag=("diagnostic_flag", lambda values: compact_flag(list(values))),
        )
    )

    sample_desc_cons = sensitivity_consistency(
        descriptor,
        sample_descriptor,
        value_col="r2_descriptor",
        out_prefix="sample_collapsed_descriptor",
        transform=transform,
    )
    complete_desc_cons = sensitivity_consistency(
        descriptor,
        complete_descriptor,
        value_col="r2_descriptor",
        out_prefix="complete_case_descriptor",
        transform=transform,
    )
    overall_coverage = coverage[
        coverage["diagnostic_type"].eq("overall_coverage") & coverage["response_transform"].eq(transform)
    ][
        [
            "regime",
            "metric",
            "n_empirical_samples",
            "k_min",
            "k_median",
            "k_max",
            "n_complete_case_instances",
            "frac_complete_case_instances",
            "n_duplicate_instance_workflow_rows",
        ]
    ].drop_duplicates()
    transform_desc = consistency[consistency["comparison_type"].eq("descriptor_response")][
        [
            "regime",
            "metric",
            "direction_agreement_fraction",
            "slope_spearman",
            "abs_slope_spearman",
            "r2_spearman",
        ]
    ].rename(
        columns={
            "direction_agreement_fraction": "transform_descriptor_direction_agreement",
            "slope_spearman": "transform_descriptor_slope_spearman",
            "abs_slope_spearman": "transform_descriptor_abs_slope_spearman",
            "r2_spearman": "transform_descriptor_R2_spearman",
        }
    )
    transform_proto = consistency[consistency["comparison_type"].eq("protocol_response")][
        ["regime", "metric", "r2_spearman"]
    ].rename(columns={"r2_spearman": "transform_protocol_R2_spearman"})

    for table in [
        sample_p,
        complete_p,
        primary_desc,
        sample_desc,
        complete_desc,
        study_p,
        study_desc,
        sample_desc_cons,
        complete_desc_cons,
    ]:
        out = out.merge(table, on=idx, how="left")
    out["sample_collapsed_protocol_R2_delta"] = out["sample_collapsed_protocol_R2"] - out["R2_protocol"]
    out["complete_case_protocol_R2_delta"] = out["complete_case_protocol_R2"] - out["R2_protocol"]
    out["sample_collapsed_descriptor_median_R2_delta"] = (
        out["sample_collapsed_descriptor_median_R2"] - out["primary_descriptor_median_R2"]
    )
    out["complete_case_descriptor_median_R2_delta"] = (
        out["complete_case_descriptor_median_R2"] - out["primary_descriptor_median_R2"]
    )
    out = out.merge(overall_coverage, on=["regime", "metric"], how="left")
    out = out.merge(transform_desc, on=["regime", "metric"], how="left")
    out = out.merge(transform_proto, on=["regime", "metric"], how="left")
    out["diagnostic_flag"] = out.apply(
        lambda row: compact_flag(
            [
                row.get("primary_diagnostic_flag"),
                row.get("sample_collapsed_protocol_diagnostic_flag"),
                row.get("complete_case_protocol_diagnostic_flag"),
                row.get("sample_collapsed_descriptor_diagnostic_flag"),
                row.get("complete_case_descriptor_diagnostic_flag"),
                row.get("study_loo_protocol_diagnostic_flag"),
                row.get("study_loo_descriptor_diagnostic_flag"),
            ]
        ),
        axis=1,
    )
    columns = [
        "regime",
        "metric",
        "workflow",
        "R2_protocol",
        "R2_all_descriptors",
        "partial_R2_protocol_after_descriptors",
        "partial_R2_within_protocol_descriptors",
        "FDR_all_descriptors",
        "FDR_protocol_after_descriptors",
        "FDR_within_protocol_descriptors",
        "sample_collapsed_protocol_R2",
        "sample_collapsed_protocol_R2_delta",
        "sample_collapsed_protocol_FDR",
        "primary_descriptor_median_R2",
        "sample_collapsed_descriptor_median_R2",
        "sample_collapsed_descriptor_median_R2_delta",
        "sample_collapsed_descriptor_R2_spearman_by_workflow",
        "complete_case_protocol_R2",
        "complete_case_protocol_R2_delta",
        "complete_case_protocol_FDR",
        "complete_case_descriptor_median_R2",
        "complete_case_descriptor_median_R2_delta",
        "complete_case_descriptor_R2_spearman_by_workflow",
        "study_loo_protocol_R2_max_abs_delta",
        "study_loo_protocol_all_significant",
        "study_loo_protocol_n_studies_tested",
        "study_loo_protocol_n_studies_ok",
        "study_loo_descriptor_median_direction_agreement",
        "study_loo_descriptor_median_R2_max_abs_delta",
        "study_loo_descriptor_min_n_studies_ok",
        "study_loo_descriptor_all_significant",
        "n_empirical_samples",
        "k_min",
        "k_median",
        "k_max",
        "n_complete_case_instances",
        "frac_complete_case_instances",
        "n_duplicate_instance_workflow_rows",
        "transform_descriptor_direction_agreement",
        "transform_descriptor_slope_spearman",
        "transform_descriptor_abs_slope_spearman",
        "transform_descriptor_R2_spearman",
        "transform_protocol_R2_spearman",
        "diagnostic_flag",
    ]
    out = out[columns].sort_values(["regime", "metric", "workflow"]).reset_index(drop=True)
    return rename_for_supplement(
        out,
        {
            "regime": "Regime",
            "metric": "Metric",
            "workflow": "Workflow",
            "R2_protocol": "Protocol-associated R2",
            "R2_all_descriptors": "All-descriptor R2",
            "partial_R2_protocol_after_descriptors": "Partial R2 protocol after descriptors",
            "partial_R2_within_protocol_descriptors": "Partial R2 within-protocol descriptors",
            "FDR_all_descriptors": "All-descriptor FDR",
            "FDR_protocol_after_descriptors": "Protocol-after-descriptors FDR",
            "FDR_within_protocol_descriptors": "Within-protocol descriptor FDR",
            "sample_collapsed_protocol_R2": "Sample-level protocol-associated R2",
            "sample_collapsed_protocol_R2_delta": "Sample-level protocol-associated R2 delta",
            "sample_collapsed_protocol_FDR": "Sample-level protocol-associated FDR",
            "primary_descriptor_median_R2": "Instance-level median descriptor R2",
            "sample_collapsed_descriptor_median_R2": "Sample-level median morphology R2",
            "sample_collapsed_descriptor_median_R2_delta": "Sample-level median morphology R2 delta",
            "sample_collapsed_descriptor_R2_spearman_by_workflow": "Sample-level morphology R2 Spearman",
            "complete_case_protocol_R2": "Complete-case protocol-associated R2",
            "complete_case_protocol_R2_delta": "Complete-case protocol-associated R2 delta",
            "complete_case_protocol_FDR": "Complete-case protocol-associated FDR",
            "complete_case_descriptor_median_R2": "Complete-case median morphology R2",
            "complete_case_descriptor_median_R2_delta": "Complete-case median morphology R2 delta",
            "complete_case_descriptor_R2_spearman_by_workflow": "Complete-case morphology R2 Spearman",
            "study_loo_protocol_R2_max_abs_delta": "Study LOO protocol-associated R2 max abs delta",
            "study_loo_protocol_all_significant": "Study LOO protocol all P<0.05",
            "study_loo_protocol_n_studies_tested": "Study LOO N studies tested",
            "study_loo_protocol_n_studies_ok": "Study LOO N studies ok",
            "study_loo_descriptor_median_direction_agreement": "Study LOO median morphology direction agreement",
            "study_loo_descriptor_median_R2_max_abs_delta": "Study LOO median morphology R2 max abs delta",
            "study_loo_descriptor_min_n_studies_ok": "Study LOO min N studies ok",
            "study_loo_descriptor_all_significant": "Study LOO morphology all P<0.05",
            "n_empirical_samples": "N empirical samples",
            "k_min": "Minimum workflows per instance",
            "k_median": "Median workflows per instance",
            "k_max": "Maximum workflows per instance",
            "n_complete_case_instances": "N complete-case instances",
            "frac_complete_case_instances": "Fraction complete-case instances",
            "n_duplicate_instance_workflow_rows": "N duplicate instance-workflow rows",
            "transform_descriptor_direction_agreement": "Transform morphology direction agreement",
            "transform_descriptor_slope_spearman": "Transform morphology slope Spearman",
            "transform_descriptor_abs_slope_spearman": "Transform morphology abs slope Spearman",
            "transform_descriptor_R2_spearman": "Transform morphology R2 Spearman",
            "transform_protocol_R2_spearman": "Transform protocol-associated R2 Spearman",
            "diagnostic_flag": "Diagnostic flag",
        },
    )


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir
    output_dir = ensure_dir(args.output_dir)

    protocol = build_protocol_table(input_dir, args.response_transform)
    descriptor = build_descriptor_table(input_dir, args.response_transform)
    sensitivity = build_sensitivity_table(input_dir, args.response_transform)
    notes = build_notes_table()

    write_tsv(protocol, output_dir / TABLE_PROTOCOL)
    write_tsv(descriptor, output_dir / TABLE_DESCRIPTOR)
    write_tsv(sensitivity, output_dir / TABLE_SENSITIVITY)
    write_xlsx(
        output_dir / WORKBOOK,
        [
            ("S8-notes", notes),
            ("S8A Protocol variation", protocol),
            ("S8B Morphology assoc", descriptor),
            ("S8C Robustness diag", sensitivity),
        ],
    )
    print(f"Wrote compact supplementary tables to: {output_dir}")


if __name__ == "__main__":
    main()
