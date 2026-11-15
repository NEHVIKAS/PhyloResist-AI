"""Validate canonical study-intake files before analysis.

Run after notebooks/13_prepare_study_intake.py, or point to another canonical
intake folder containing a study_config.json:
    python notebooks/14_validate_study_intake.py
    python notebooks/14_validate_study_intake.py --intake-dir data/processed/other_study

This checks the derived intake tables against each other and against the source
genotype matrix. It writes a CSV of individual checks plus a readable report.
It does not modify source data or assign resistance labels.
"""

from __future__ import annotations

import json
import argparse
import math
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INTAKE = ROOT / "data" / "processed" / "fg_graminearum_2025"
REPORTS = ROOT / "results" / "reports"

MANIFEST_REQUIRED = {
    "study_id", "pathogen", "scientific_name", "isolate_id", "alternate_ids",
    "collection_year", "collection_region", "host", "lineage_or_population",
}
PHENOTYPE_REQUIRED = {
    "study_id", "isolate_id", "fungicide", "trait_name", "value", "unit",
    "measurement_type", "replicate_count", "standard_error", "degrees_of_freedom",
    "source_lower_CL_u", "source_pper_CL", "resistance_threshold",
}
def add(checks: list[dict[str, str]], name: str, status: str, details: str) -> None:
    checks.append({"check": name, "status": status, "details": details})


def read_required(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Intake file not found: {path}")
    return pd.read_csv(path, dtype=str, keep_default_na=True, low_memory=False)


def normalize_id(value: Any) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def validate(intake: Path) -> tuple[list[dict[str, str]], dict[str, Any]]:
    checks: list[dict[str, str]] = []
    config_path = intake / "study_config.json"
    if not config_path.is_file():
        raise FileNotFoundError(f"Study configuration not found: {config_path}. Run the study intake preparation step first.")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    required_config = {"study_id", "pathogen", "fungicides", "phenotype", "genotype_manifest", "phenotype_table", "isolate_manifest"}
    missing_config = sorted(required_config - set(config))
    if missing_config:
        raise ValueError(f"Study configuration is missing required fields: {missing_config}")
    manifest_path = intake / config["isolate_manifest"]
    phenotype_path = intake / config["phenotype_table"]
    genotype_manifest_path = intake / config["genotype_manifest"]
    expected_fungicides = {str(value).strip().casefold() for value in config["fungicides"]}
    manifest = read_required(manifest_path)
    phenotypes = read_required(phenotype_path)
    if not genotype_manifest_path.is_file():
        raise FileNotFoundError(f"Genotype input manifest not found: {genotype_manifest_path}")
    genotype_info = json.loads(genotype_manifest_path.read_text(encoding="utf-8"))

    missing_manifest_columns = sorted(MANIFEST_REQUIRED - set(manifest.columns))
    missing_phenotype_columns = sorted(PHENOTYPE_REQUIRED - set(phenotypes.columns))
    add(checks, "Manifest required columns", "FAIL" if missing_manifest_columns else "PASS", f"Missing={missing_manifest_columns}")
    add(checks, "Phenotype required columns", "FAIL" if missing_phenotype_columns else "PASS", f"Missing={missing_phenotype_columns}")
    if missing_manifest_columns or missing_phenotype_columns:
        return checks, {"manifest_rows": len(manifest), "phenotype_rows": len(phenotypes), "status": "FAIL"}

    manifest["isolate_id"] = manifest["isolate_id"].map(normalize_id)
    phenotypes["isolate_id"] = phenotypes["isolate_id"].map(normalize_id)
    duplicate_manifest = int(manifest["isolate_id"].duplicated().sum())
    duplicate_records = int(phenotypes.duplicated(subset=["study_id", "isolate_id", "fungicide", "trait_name"]).sum())
    add(checks, "Unique manifest isolate IDs", "FAIL" if duplicate_manifest else "PASS", f"Duplicate rows={duplicate_manifest}")
    add(checks, "Unique phenotype isolate/fungicide/trait records", "FAIL" if duplicate_records else "PASS", f"Duplicate rows={duplicate_records}")

    manifest_ids = set(manifest["isolate_id"]) - {""}
    phenotype_ids = set(phenotypes["isolate_id"]) - {""}
    only_manifest = sorted(manifest_ids - phenotype_ids)
    only_phenotype = sorted(phenotype_ids - manifest_ids)
    add(
        checks,
        "Manifest/phenotype isolate matching",
        "FAIL" if only_manifest or only_phenotype else "PASS",
        f"Manifest-only={len(only_manifest)}; phenotype-only={len(only_phenotype)}; shared={len(manifest_ids & phenotype_ids)}",
    )

    fungicides = set(phenotypes["fungicide"].dropna().str.strip().str.casefold())
    add(checks, "Configured fungicide coverage", "FAIL" if fungicides != expected_fungicides else "PASS", f"Expected={sorted(expected_fungicides)}; observed={sorted(fungicides)}")
    study_ids = set(phenotypes["study_id"].dropna().str.strip()) | set(manifest["study_id"].dropna().str.strip())
    add(checks, "Configured study identifier", "FAIL" if study_ids != {str(config["study_id"])} else "PASS", f"Expected={config['study_id']!r}; observed={sorted(study_ids)}")
    pathogens = set(manifest["scientific_name"].dropna().str.strip())
    add(checks, "Configured pathogen identity", "FAIL" if pathogens != {str(config["pathogen"])} else "PASS", f"Expected={config['pathogen']!r}; observed={sorted(pathogens)}")
    expected_trait = str(config["phenotype"].get("trait_name", "")).strip()
    expected_unit = str(config["phenotype"].get("unit", "")).strip()
    traits = set(phenotypes["trait_name"].dropna().str.strip())
    units = set(phenotypes["unit"].dropna().str.strip())
    add(checks, "Configured phenotype trait", "FAIL" if traits != {expected_trait} else "PASS", f"Expected={expected_trait!r}; observed={sorted(traits)}")

    values = pd.to_numeric(phenotypes["value"], errors="coerce")
    non_numeric = int(phenotypes["value"].notna().sum() - values.notna().sum())
    missing_values = int(phenotypes["value"].isna().sum())
    non_finite = int((~values.dropna().map(math.isfinite)).sum())
    add(checks, "Phenotype values finite and numeric", "FAIL" if non_numeric or missing_values or non_finite else "PASS", f"Missing={missing_values}; non-numeric={non_numeric}; non-finite={non_finite}")
    add(checks, "Configured phenotype units", "FAIL" if units != {expected_unit} or not expected_unit else "PASS", f"Expected={expected_unit!r}; observed={sorted(units)}")
    threshold_count = int(phenotypes["resistance_threshold"].notna().sum())
    if threshold_count == 0:
        add(checks, "Resistance threshold supplied", "INFO", "No threshold supplied; report measured EC50 values and do not assign susceptible/resistant labels.")
    else:
        add(checks, "Resistance threshold source", "WARN", "Threshold fields exist; confirm a validated pathogen/fungicide/assay-specific source before classification.")

    source_file = genotype_info.get("source_file", "")
    source_path = (ROOT / source_file).resolve() if source_file else None
    source_checks = "not checked"
    if not source_path or not source_path.is_file():
        add(checks, "Genotype source file available", "WARN", f"Could not locate source file from manifest: {source_file!r}")
    else:
        genotype = pd.read_csv(source_path, dtype=str, keep_default_na=True, low_memory=False)
        source_ids = [normalize_id(value) for value in genotype.columns[3:]]
        source_variants = int(genotype.shape[0])
        source_codes = sorted(set(genotype.iloc[:, 3:].stack().dropna().str.strip().str.upper().tolist()))
        expected_codes = sorted(genotype_info.get("observed_genotype_codes", []))
        count_match = len(source_ids) == genotype_info.get("isolate_count") and source_variants == genotype_info.get("variant_count")
        code_match = source_codes == expected_codes
        id_match = set(source_ids) == manifest_ids
        add(checks, "Genotype source dimensions match manifest", "PASS" if count_match else "FAIL", f"Source={source_variants} variants × {len(source_ids)} isolates; manifest={genotype_info.get('variant_count')} × {genotype_info.get('isolate_count')}")
        add(checks, "Genotype source isolate IDs match manifest", "PASS" if id_match else "FAIL", f"Source-only={len(set(source_ids) - manifest_ids)}; manifest-only={len(manifest_ids - set(source_ids))}")
        add(checks, "Genotype coding matches manifest", "PASS" if code_match else "FAIL", f"Source={source_codes}; manifest={expected_codes}")
        source_checks = "checked"

    build = str(genotype_info.get("genome_build", "unknown"))
    if "not specified" in build.casefold() or build.casefold() == "unknown":
        add(checks, "Genome build identified", "WARN", "Genome build is unknown; resolve before cross-study coordinate harmonization.")
    else:
        add(checks, "Genome build identified", "PASS", build)
    accession_ids = genotype_info.get("sequence_accessions_observed_in_variant_ids", [])
    reference_hypothesis = genotype_info.get("coordinate_reference_hypothesis", "")
    if accession_ids:
        add(checks, "Coordinate accession evidence", "INFO", f"Observed={accession_ids}. {reference_hypothesis}")

    failed = [row for row in checks if row["status"] == "FAIL"]
    summary = {
        "manifest_rows": len(manifest),
        "phenotype_rows": len(phenotypes),
        "isolate_count": len(manifest_ids),
        "fungicides": sorted(fungicides),
        "study_id": config["study_id"],
        "pathogen": config["pathogen"],
        "failed_checks": len(failed),
        "source_matrix_checked": source_checks,
        "status": "FAIL" if failed else "PASS WITH WARNINGS" if any(row["status"] == "WARN" for row in checks) else "PASS",
    }
    return checks, summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--intake-dir",
        type=Path,
        default=DEFAULT_INTAKE,
        help="Folder containing study_config.json and the canonical study files.",
    )
    args = parser.parse_args()
    intake = args.intake_dir.resolve()
    checks, summary = validate(intake)
    REPORTS.mkdir(parents=True, exist_ok=True)
    CHECKS_PATH = REPORTS / "study_intake_validation_checks.csv"
    REPORT_PATH = REPORTS / "study_intake_validation.md"
    pd.DataFrame(checks).to_csv(CHECKS_PATH, index=False)
    lines = [
        "# Study intake validation",
        "",
        f"Generated: {datetime.now().astimezone().isoformat(timespec='seconds')}",
        f"Study: `{summary.get('study_id', 'unknown')}` ({summary.get('pathogen', 'unknown')})",
        f"Status: **{summary['status']}**",
        "",
        "## Summary",
        "",
        f"- Isolates: {summary.get('isolate_count', 'not available')}",
        f"- Phenotype rows: {summary.get('phenotype_rows', 'not available')}",
        f"- Fungicides: {', '.join(summary.get('fungicides', []))}",
        f"- Failed checks: {summary.get('failed_checks', 0)}",
        f"- Original genotype matrix checked: {summary.get('source_matrix_checked', 'not checked')}",
        "",
        "## Checks",
        "",
        "| Status | Check | Details |",
        "|---|---|---|",
    ]
    for row in checks:
        details = row["details"].replace("|", "\\|")
        lines.append(f"| {row['status']} | {row['check']} | {details} |")
    lines.extend([
        "",
        "This is an intake-format/QC check, not an external validation or resistance classifier.",
        "A missing threshold means no susceptible/resistant label should be assigned.",
    ])
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Study intake validation: {summary['status']}")
    print(f"Report: {REPORT_PATH}")
    print(f"Check table: {CHECKS_PATH}")
    return 1 if summary["status"] == "FAIL" else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"Study intake validation failed: {error}", file=sys.stderr)
        raise
