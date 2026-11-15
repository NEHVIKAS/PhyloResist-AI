"""Prepare the public Zymoseptoria tritici European-panel phenotype intake.

Run from the repository root in the phyloresist environment:
    python notebooks/15_prepare_ztritici_phenotype_intake.py

This reads the source-preserved Table S1 workbook in data/reference and writes
derived isolate and EC50 phenotype tables. It does not download sequencing data,
call variants, run GWAS, or create resistance labels. The source workbook remains
unchanged.
"""

from __future__ import annotations

import re
import sys
from datetime import datetime
from pathlib import Path
from zipfile import ZipFile
from xml.etree import ElementTree as ET

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "reference" / "zymoseptoria_tritici_european_panel_table_s1.xlsx"
OUTPUT = ROOT / "data" / "processed" / "z_tritici_european_panel_2025"
REPORTS = ROOT / "results" / "reports"
STUDY_ID = "z_tritici_european_dmi_panel_2025"
PATHOGEN = "Zymoseptoria tritici"
EC50_COLUMNS = {
    "EC50_CCZ": "cyproconazole",
    "EC50_EPZ": "epoxiconazole",
    "EC50_PLZ": "prochloraz",
    "EC50_TCZ": "tebuconazole",
    "EC50_PTZ": "prothioconazole",
}
MISSING_VALUES = {"", "NA", "Nd", "nd", "No"}
EXCEL_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"


def excel_column_index(reference: str) -> int:
    """Convert an Excel cell reference such as AB12 to a zero-based column."""
    match = re.match(r"[A-Z]+", reference.upper())
    if not match:
        raise ValueError(f"Invalid Excel cell reference: {reference!r}")
    index = 0
    for letter in match.group(0):
        index = index * 26 + ord(letter) - ord("A") + 1
    return index - 1


def read_first_worksheet(path: Path) -> pd.DataFrame:
    """Read an .xlsx worksheet using the standard library, including blank cells."""
    if not path.is_file():
        raise FileNotFoundError(f"Source workbook not found: {path}")
    namespace = f"{{{EXCEL_NS}}}"
    with ZipFile(path) as workbook:
        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in workbook.namelist():
            root = ET.fromstring(workbook.read("xl/sharedStrings.xml"))
            shared_strings = [
                "".join(node.text or "" for node in item.iter(namespace + "t"))
                for item in root.findall(namespace + "si")
            ]
        sheet = ET.fromstring(workbook.read("xl/worksheets/sheet1.xml"))
        row_values: list[dict[int, str]] = []
        for row in sheet.findall(f".//{namespace}row"):
            values: dict[int, str] = {}
            for cell in row.findall(namespace + "c"):
                value_node = cell.find(namespace + "v")
                value = value_node.text if value_node is not None else ""
                if cell.attrib.get("t") == "s" and value:
                    value = shared_strings[int(value)]
                values[excel_column_index(cell.attrib.get("r", ""))] = value
            if values:
                row_values.append(values)
    if not row_values:
        raise ValueError("The source workbook has no data rows.")
    headers = [row_values[0].get(index, "").strip() for index in range(max(row_values[0]) + 1)]
    if "Isolate_name" not in headers:
        raise ValueError(f"Expected Isolate_name in Table S1 headers, found: {headers}")
    records = [
        {headers[index]: row.get(index, "") for index in range(len(headers))}
        for row in row_values[1:]
    ]
    return pd.DataFrame(records)


def present(value: object) -> bool:
    """Return whether a table value represents an observed measurement."""
    return str(value).strip() not in MISSING_VALUES


def main() -> int:
    source = read_first_worksheet(SOURCE)
    required = {"Isolate_name", "State", "Year", "Latitude", "Longitude", "BioProject_accession", *EC50_COLUMNS}
    missing_columns = sorted(required - set(source.columns))
    if missing_columns:
        raise ValueError(f"Table S1 is missing expected columns: {missing_columns}")

    source["Isolate_name"] = source["Isolate_name"].astype(str).str.strip()
    if source["Isolate_name"].eq("").any() or source["Isolate_name"].duplicated().any():
        raise ValueError("Table S1 has blank or duplicate isolate names; intake stopped.")

    OUTPUT.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)
    manifest_path = OUTPUT / "isolate_manifest.csv"
    phenotype_path = OUTPUT / "phenotype_measurements.csv"
    coverage_path = OUTPUT / "phenotype_coverage.csv"
    report_path = REPORTS / "ztritici_european_phenotype_intake.md"

    manifest = pd.DataFrame({
        "study_id": STUDY_ID,
        "pathogen": PATHOGEN,
        "scientific_name": PATHOGEN,
        "isolate_id": source["Isolate_name"],
        "alternate_ids": "",
        "collection_year": source["Year"],
        "collection_country": source["State"],
        "collection_region": "",
        "host": "wheat (study context)",
        "lineage_or_population": source["Continent"],
        "source_bioproject": source["BioProject_accession"],
        "latitude": source["Latitude"],
        "longitude": source["Longitude"],
        "included_in_snp_gwas": source["Strains_in_SNPs_GWAS"],
        "source_reference": "Puccetti et al. 2025/2026, Zenodo Table S1",
        "notes": "Source-preserved metadata; no sequencing-run crosswalk included in this intake.",
    })
    manifest.to_csv(manifest_path, index=False)

    long_frames: list[pd.DataFrame] = []
    coverage_rows: list[dict[str, object]] = []
    for source_column, fungicide in EC50_COLUMNS.items():
        values = source[source_column].astype(str).str.strip()
        keep = values.map(present)
        numeric = pd.to_numeric(values.where(keep), errors="coerce")
        invalid = int(keep.sum() - numeric.notna().sum())
        if invalid:
            raise ValueError(f"{source_column} has {invalid} non-numeric non-missing values.")
        coverage_rows.append({
            "fungicide": fungicide,
            "source_column": source_column,
            "n_ec50": int(keep.sum()),
            "missing_ec50": int((~keep).sum()),
            "minimum_mg_per_l": float(numeric.min()),
            "maximum_mg_per_l": float(numeric.max()),
        })
        long_frames.append(pd.DataFrame({
            "study_id": STUDY_ID,
            "isolate_id": source.loc[keep, "Isolate_name"],
            "fungicide": fungicide,
            "trait_name": "EC50",
            "value": values.loc[keep],
            "unit": "mg/L",
            "measurement_type": "study-reported mean EC50",
            "replicate_id": "",
            "replicate_count": "2 biological × 2 technical replicates",
            "assay_protocol": "Liquid assay; 6 days; non-linear curve fit. See Puccetti et al. methods.",
            "assay_batch": "",
            "standard_error": "",
            "degrees_of_freedom": "",
            "source_lower_CL_u": "",
            "source_pper_CL": "",
            "resistance_threshold": "",
            "threshold_source": "",
            "source_reference": "Puccetti et al. 2025/2026, Zenodo Table S1",
            "notes": "Separate pathogen/study module. Values remain in source mg/L units; no transformation or resistance label applied.",
        }))
    phenotypes = pd.concat(long_frames, ignore_index=True)
    phenotypes.to_csv(phenotype_path, index=False)
    coverage = pd.DataFrame(coverage_rows)
    coverage.to_csv(coverage_path, index=False)

    lines = [
        "# Zymoseptoria tritici European-panel phenotype intake",
        "",
        f"Generated: {datetime.now().astimezone().isoformat(timespec='seconds')}",
        "",
        "## What was prepared",
        "",
        f"- Source Table S1 isolates: {len(source):,}",
        f"- Canonical manifest rows: {len(manifest):,}",
        f"- EC50 phenotype rows: {len(phenotypes):,}",
        f"- Country/location labels represented: {source['State'].replace('', pd.NA).nunique()}",
        f"- Distinct BioProject accessions: {source['BioProject_accession'].map(present).sum()} entries; {source.loc[source['BioProject_accession'].map(present), 'BioProject_accession'].nunique()} unique IDs",
        "",
        "## EC50 coverage",
        "",
        "| Fungicide | Source column | EC50 values | Missing | Minimum (mg/L) | Maximum (mg/L) |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for row in coverage_rows:
        lines.append(
            f"| {row['fungicide']} | {row['source_column']} | {row['n_ec50']} | {row['missing_ec50']} | "
            f"{row['minimum_mg_per_l']:.6g} | {row['maximum_mg_per_l']:.6g} |"
        )
    lines.extend([
        "",
        "## Current boundary",
        "",
        "This is a phenotype and metadata intake only. The source table supplies BioProject-level accessions, not a ready-to-use genotype matrix/VCF or isolate-to-run crosswalk. Do not run genotype prediction or GWAS from these outputs. Obtain the final filtered VCF/genotype matrix, sequencing-run crosswalk, exact IPO323 reference build/version, and reuse terms first.",
        "",
        "The EC50 values are reported in mg/L and come from a different pathogen, assay, and unit than the current *Fusarium graminearum* study. Do not pool them or compare raw values without an explicit harmonization plan.",
        "",
        "The source workbook was read only and was not modified.",
    ])
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("Zymoseptoria tritici phenotype intake completed.")
    print(f"Manifest: {manifest_path}")
    print(f"Phenotype measurements: {phenotype_path}")
    print(f"Coverage table: {coverage_path}")
    print(f"Report: {report_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"Zymoseptoria tritici phenotype intake failed: {error}", file=sys.stderr)
        raise
