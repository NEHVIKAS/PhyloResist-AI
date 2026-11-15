"""Beginner-friendly integrity checks for the PhyloResist AI starter dataset.

Run from the repository root with:
    python notebooks/03_data_integrity.py

This script reads the three CSV files in data/raw and writes a human-readable
report plus a machine-readable check summary under results/reports. It does not
change raw files. The study author confirmed `trait1` is an isolate-level
least-squares mean EC50 in mM across three experimental replications. Values
remain on their supplied scale; no log transformation is applied.
"""

from __future__ import annotations

import csv
import math
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd


# Resolve paths from this script, so it also works when launched from elsewhere.
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = PROJECT_ROOT / "data" / "raw"
REPORT_DIR = PROJECT_ROOT / "results" / "reports"
GENOTYPE_FILE = RAW_DIR / "genotype.csv"
PHENOTYPE_FILES = {
    "propiconazole": RAW_DIR / "phenotype_propiconazole.csv",
    "tebuconazole": RAW_DIR / "phenotype_tebuconazole.csv",
}
METADATA_COLUMNS = 3  # rs#, chrom, pos in the current genotype file.
EXPECTED_ALLELES = {"A", "C", "G", "T"}


def normalize_id(value: Any) -> str:
    """Normalize IDs without changing meaningful leading zeros or case."""
    if pd.isna(value):
        return ""
    value = str(value).strip()
    # CSV readers can represent integer-like IDs as e.g. 11157.0.
    if re.fullmatch(r"\d+\.0+", value):
        value = value.split(".", 1)[0]
    return value


def add_check(checks: list[dict[str, str]], name: str, status: str, details: str) -> None:
    checks.append({"check": name, "status": status, "details": details})


def duplicate_values(values: list[str]) -> list[str]:
    seen: set[str] = set()
    duplicates: set[str] = set()
    for value in values:
        if value in seen:
            duplicates.add(value)
        seen.add(value)
    return sorted(duplicates)


def main() -> int:
    checks: list[dict[str, str]] = []
    lines: list[str] = []
    lines.append("PHYLORESIST AI — DATA INTEGRITY REPORT")
    lines.append("=" * 72)
    lines.append(f"Generated: {datetime.now().astimezone().isoformat(timespec='seconds')}")
    lines.append(f"Project: {PROJECT_ROOT}")
    lines.append("")
    lines.append("Scientific interpretation safeguard")
    lines.append("-" * 72)
    lines.append("The study author confirmed `trait1` is an isolate-level least-squares mean EC50 in mM.")
    lines.append("Each isolate estimate summarizes three experimental replications; values remain on the supplied scale.")
    lines.append("No log transformation is applied. Review the final peer-reviewed methods for full curve-fit exclusions.")
    lines.append("")

    required = [GENOTYPE_FILE, *PHENOTYPE_FILES.values()]
    missing_files = [path for path in required if not path.is_file()]
    if missing_files:
        for path in missing_files:
            add_check(checks, f"Input file exists: {path.name}", "FAIL", f"Not found: {path}")
        lines.extend(["INPUT FILE ERROR", "-" * 72])
        lines.extend(f"Missing: {path}" for path in missing_files)
        REPORT_DIR.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(checks).to_csv(REPORT_DIR / "integrity_checks.csv", index=False)
        (REPORT_DIR / "integrity_report.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
        print("Required input files are missing. See results/reports/integrity_report.txt")
        return 1

    # Read the header separately so duplicate raw CSV column labels are detectable.
    with GENOTYPE_FILE.open("r", newline="", encoding="utf-8-sig") as handle:
        genotype_header = next(csv.reader(handle))
    genotype = pd.read_csv(GENOTYPE_FILE, dtype=str, keep_default_na=True, low_memory=False)
    genotype_metadata = list(genotype.columns[:METADATA_COLUMNS])
    genotype_isolates = list(genotype.columns[METADATA_COLUMNS:])
    genotype_ids = [normalize_id(value) for value in genotype_isolates]
    duplicate_headers = duplicate_values([str(value).strip() for value in genotype_header])
    duplicate_genotype_ids = duplicate_values(genotype_ids)

    lines.extend(["GENOTYPE DATA", "-" * 72])
    lines.append(f"Shape: {genotype.shape[0]:,} variant rows × {genotype.shape[1]:,} columns")
    lines.append(f"Metadata columns (first three): {genotype_metadata}")
    lines.append(f"Isolate columns identified: {len(genotype_isolates):,}")
    add_check(checks, "Genotype metadata columns", "INFO", f"First three columns: {genotype_metadata}")
    add_check(checks, "Genotype isolate columns", "INFO", f"Identified {len(genotype_isolates)} columns after metadata")
    add_check(checks, "Duplicate genotype header labels", "FAIL" if duplicate_headers else "PASS", f"Count={len(duplicate_headers)}; examples={duplicate_headers[:10]}")
    add_check(checks, "Duplicate genotype isolate IDs", "FAIL" if duplicate_genotype_ids else "PASS", f"Count={len(duplicate_genotype_ids)}; examples={duplicate_genotype_ids[:10]}")
    lines.append(f"Duplicate header labels: {len(duplicate_headers)}")
    lines.append(f"Duplicate isolate IDs: {len(duplicate_genotype_ids)}")

    # Missing values, alleles, and metadata validation.
    missing_by_column = genotype.isna().sum()
    missing_cells = int(missing_by_column.sum())
    add_check(checks, "Genotype missing cells", "PASS" if not missing_cells else "WARN", f"{missing_cells} missing cells")
    lines.append(f"Missing genotype cells: {missing_cells:,}")
    allele_cells = genotype[genotype_isolates].apply(lambda col: col.astype("string").str.strip().str.upper())
    observed_alleles = set(allele_cells.stack().dropna().tolist())
    unexpected_alleles = sorted(observed_alleles - EXPECTED_ALLELES)
    add_check(checks, "Allele coding", "PASS" if not unexpected_alleles else "WARN", f"Observed={sorted(observed_alleles)}; unexpected={unexpected_alleles[:25]}")
    lines.append(f"Observed genotype codes: {sorted(observed_alleles)}")
    lines.append(f"Unexpected genotype codes (outside A/C/G/T): {unexpected_alleles[:25]}")

    metadata_missing = {column: int(genotype[column].isna().sum()) for column in genotype_metadata}
    metadata_duplicates = int(genotype.duplicated(subset=genotype_metadata).sum()) if genotype_metadata else 0
    add_check(checks, "Variant metadata missing values", "PASS" if not any(metadata_missing.values()) else "WARN", str(metadata_missing))
    add_check(checks, "Duplicate variant metadata rows", "PASS" if not metadata_duplicates else "WARN", f"{metadata_duplicates} repeated rows across {genotype_metadata}")
    lines.append(f"Missing metadata by column: {metadata_missing}")
    lines.append(f"Duplicate variant metadata combinations: {metadata_duplicates:,}")
    if "pos" in genotype.columns:
        numeric_pos = pd.to_numeric(genotype["pos"], errors="coerce")
        invalid_pos = int(numeric_pos.isna().sum())
        nonpositive_pos = int((numeric_pos.dropna() <= 0).sum())
        add_check(checks, "Variant position values", "PASS" if not invalid_pos and not nonpositive_pos else "WARN", f"Non-numeric/missing={invalid_pos}; non-positive={nonpositive_pos}")
        lines.append(f"Position check: non-numeric/missing={invalid_pos}; non-positive={nonpositive_pos}")

    phenotype_data: dict[str, tuple[pd.DataFrame, str, list[str]]] = {}
    for fungicide, path in PHENOTYPE_FILES.items():
        phenotype = pd.read_csv(path, dtype=str, keep_default_na=True, low_memory=False)
        id_candidates = [column for column in phenotype.columns if any(term in str(column).strip().casefold() for term in ("isolate", "sample", "strain", "phenotype", "id"))]
        # Prefer an explicit ID-like label, then a phenotype-labeled field, then column 1.
        preferred = [column for column in id_candidates if re.search(r"(^|[^a-z])(id|isolate|sample|strain)([^a-z]|$)", str(column).strip().casefold())]
        id_column = (preferred or id_candidates or list(phenotype.columns[:1]))[0]
        phenotype_ids = [normalize_id(value) for value in phenotype[id_column].tolist()]
        phenotype_duplicates = duplicate_values(phenotype_ids)
        genotype_set = {value for value in genotype_ids if value}
        phenotype_set = {value for value in phenotype_ids if value}
        only_genotype = sorted(genotype_set - phenotype_set)
        only_phenotype = sorted(phenotype_set - genotype_set)
        missing_cells = int(phenotype.isna().sum().sum())
        lines.extend(["", f"PHENOTYPE DATA — {fungicide.upper()}", "-" * 72])
        lines.append(f"Shape: {phenotype.shape[0]:,} rows × {phenotype.shape[1]:,} columns")
        lines.append(f"Isolate ID column detected: {id_column!r}; candidate labels: {id_candidates}")
        lines.append(f"Duplicate phenotype IDs: {len(phenotype_duplicates)}; examples={phenotype_duplicates[:10]}")
        lines.append(f"Missing phenotype cells: {missing_cells:,}")
        lines.append(f"IDs in genotype but absent from phenotype: {len(only_genotype)}; examples={only_genotype[:15]}")
        lines.append(f"IDs in phenotype but absent from genotype: {len(only_phenotype)}; examples={only_phenotype[:15]}")
        add_check(checks, f"{fungicide}: phenotype ID column", "INFO", f"Detected {id_column!r} from columns {list(phenotype.columns)}")
        add_check(checks, f"{fungicide}: duplicate phenotype IDs", "FAIL" if phenotype_duplicates else "PASS", f"Count={len(phenotype_duplicates)}; examples={phenotype_duplicates[:10]}")
        add_check(checks, f"{fungicide}: genotype/phenotype ID match", "PASS" if not only_genotype and not only_phenotype else "WARN", f"Genotype-only={len(only_genotype)}; phenotype-only={len(only_phenotype)}; shared={len(genotype_set & phenotype_set)}")
        add_check(checks, f"{fungicide}: missing values", "PASS" if not missing_cells else "WARN", f"{missing_cells} missing cells")
        for column in phenotype.columns:
            if column == id_column:
                continue
            numeric = pd.to_numeric(phenotype[column], errors="coerce")
            source_nonmissing = int(phenotype[column].notna().sum())
            invalid_numeric = source_nonmissing - int(numeric.notna().sum())
            finite = numeric.dropna()
            if len(finite):
                stats = {
                    "n": int(finite.size), "mean": float(finite.mean()), "std": float(finite.std()),
                    "min": float(finite.min()), "q25": float(finite.quantile(0.25)),
                    "median": float(finite.median()), "q75": float(finite.quantile(0.75)), "max": float(finite.max()),
                    "non_finite": int((~finite.map(math.isfinite)).sum()),
                }
                stat_text = ", ".join(f"{key}={value:.6g}" if isinstance(value, float) else f"{key}={value}" for key, value in stats.items())
            else:
                stats = {"n": 0, "mean": float("nan"), "std": float("nan"), "min": float("nan"), "q25": float("nan"), "median": float("nan"), "q75": float("nan"), "max": float("nan"), "non_finite": 0}
                stat_text = "no numeric observations"
            lines.append(f"Numeric field {column!r}: {stat_text}; non-numeric values={invalid_numeric}")
            add_check(checks, f"{fungicide}: numeric distribution {column}", "WARN" if invalid_numeric or not len(finite) else "INFO", f"{stat_text}; non-numeric={invalid_numeric}")
        phenotype_data[fungicide] = (phenotype, id_column, phenotype_ids)

    # Save an isolate-level reconciliation table to make mismatches actionable.
    all_ids = sorted({value for value in genotype_ids if value} | {value for _, _, ids in phenotype_data.values() for value in ids if value})
    id_report = pd.DataFrame({"isolate_id": all_ids, "in_genotype": [value in set(genotype_ids) for value in all_ids]})
    for fungicide, (_, _, ids) in phenotype_data.items():
        id_report[f"in_{fungicide}"] = [value in set(ids) for value in all_ids]
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    checks_path = REPORT_DIR / "integrity_checks.csv"
    id_path = REPORT_DIR / "isolate_id_matching.csv"
    report_path = REPORT_DIR / "integrity_report.txt"
    pd.DataFrame(checks).to_csv(checks_path, index=False)
    id_report.to_csv(id_path, index=False)
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("Integrity checks completed.")
    print(f"Text report: {report_path}")
    print(f"Check summary CSV: {checks_path}")
    print(f"Isolate matching CSV: {id_path}")
    print("Phenotype meaning and units were confirmed by the study author; consult final paper methods for full estimation details.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"Data integrity check failed: {error}", file=sys.stderr)
        raise
