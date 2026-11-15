"""Map the first study to the shared intake format without changing source files.

Run from the repository root in the phyloresist environment:
    python notebooks/13_prepare_study_intake.py

Outputs are written under data/processed/fg_graminearum_2025/ and
results/reports/. The original genotype and phenotype CSVs are read-only inputs.
"""

from __future__ import annotations

import json
import math
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
REFERENCE = ROOT / "data" / "reference"
OUTPUT = ROOT / "data" / "processed" / "fg_graminearum_2025"
REPORTS = ROOT / "results" / "reports"
STUDY_ID = "fg_graminearum_fungicide_study_2025"
PATHOGEN = "Fusarium graminearum"
GENOTYPE_PATH = RAW / "genotype.csv"
METADATA_PATH = REFERENCE / "isolate_metadata.csv"
PHENOTYPE_PATHS = {
    "propiconazole": RAW / "phenotype_propiconazole.csv",
    "tebuconazole": RAW / "phenotype_tebuconazole.csv",
}


def normalize_id(value: Any) -> str:
    """Normalize numeric-looking IDs while preserving leading zeroes."""
    if pd.isna(value):
        return ""
    value = str(value).strip()
    if re.fullmatch(r"\d+\.0+", value):
        return value.split(".", 1)[0]
    return value


def find_id_column(frame: pd.DataFrame) -> str:
    """Find the phenotype ID column despite source label/case differences."""
    preferred = [
        col for col in frame.columns
        if re.search(r"(^|[^a-z])(id|isolate|sample|strain)([^a-z]|$)", str(col).strip().casefold())
    ]
    candidates = preferred or [
        col for col in frame.columns
        if any(term in str(col).strip().casefold() for term in ("phenotype", "isolate", "sample", "strain", "id"))
    ]
    if not candidates:
        raise ValueError(f"Could not identify phenotype isolate-ID column in {list(frame.columns)}")
    return candidates[0]


def read_csv(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Required source file not found: {path}")
    return pd.read_csv(path, dtype=str, keep_default_na=True, low_memory=False)


def main() -> int:
    genotype = read_csv(GENOTYPE_PATH)
    metadata = read_csv(METADATA_PATH)
    if genotype.shape[1] < 4:
        raise ValueError("Genotype matrix must have three variant metadata columns plus isolate columns.")

    variant_columns = list(genotype.columns[:3])
    genotype_ids = [normalize_id(value) for value in genotype.columns[3:]]
    if len(set(genotype_ids)) != len(genotype_ids):
        raise ValueError("Duplicate genotype isolate IDs found; intake stopped.")
    if not {"isolate_id", "alternate_id"}.issubset(metadata.columns):
        raise ValueError("isolate_metadata.csv must contain isolate_id and alternate_id.")
    metadata["isolate_id"] = metadata["isolate_id"].map(normalize_id)
    if metadata["isolate_id"].duplicated().any():
        raise ValueError("Duplicate isolate IDs found in source metadata; intake stopped.")
    genotype_set = set(genotype_ids)
    metadata_set = set(metadata["isolate_id"])
    missing_metadata = sorted(genotype_set - metadata_set)
    extra_metadata = sorted(metadata_set - genotype_set)
    if missing_metadata or extra_metadata:
        raise ValueError(
            f"Metadata/genotype IDs do not match exactly. Missing metadata={missing_metadata[:10]}; "
            f"metadata-only={extra_metadata[:10]}"
        )

    # Canonical isolate manifest: one row per isolate, retaining source metadata.
    manifest_columns = {
        "alternate_id": "alternate_ids",
        "collection_state": "collection_region",
        "collection_year": "collection_year",
        "host": "host",
        "population_label": "lineage_or_population",
        "source_table": "source_reference",
    }
    manifest = metadata.rename(columns=manifest_columns).copy()
    manifest.insert(0, "study_id", STUDY_ID)
    manifest.insert(1, "pathogen", PATHOGEN)
    manifest.insert(2, "scientific_name", PATHOGEN)
    keep_manifest = [
        "study_id", "pathogen", "scientific_name", "isolate_id", "alternate_ids",
        "collection_year", "collection_country", "collection_region", "host",
        "lineage_or_population", "source_reference", "notes",
    ]
    for column in keep_manifest:
        if column not in manifest.columns:
            manifest[column] = ""
    manifest = manifest[keep_manifest].sort_values("isolate_id", kind="stable")

    # Canonical long phenotype table; author-confirmed values are LS means in mM.
    # Values are kept as supplied: no log transform, EC50 conversion, or threshold.
    phenotype_frames: list[pd.DataFrame] = []
    numeric_check_rows: list[dict[str, Any]] = []
    for fungicide, path in PHENOTYPE_PATHS.items():
        source = read_csv(path)
        id_column = find_id_column(source)
        source = source.rename(columns={id_column: "isolate_id"})
        source["isolate_id"] = source["isolate_id"].map(normalize_id)
        if source["isolate_id"].duplicated().any():
            raise ValueError(f"Duplicate phenotype IDs in {path.name}; intake stopped.")
        phenotype_ids = set(source["isolate_id"])
        if phenotype_ids != genotype_set:
            raise ValueError(
                f"{fungicide} phenotype IDs do not match genotype IDs: "
                f"genotype-only={len(genotype_set - phenotype_ids)}, phenotype-only={len(phenotype_ids - genotype_set)}"
            )
        if "trait1" not in source.columns:
            raise ValueError(f"trait1 was not found in {path.name}.")
        for column in ("trait1", "SE", "df", "lower.CL u", "pper.CL"):
            if column in source.columns:
                numeric = pd.to_numeric(source[column], errors="coerce")
                bad = int(source[column].notna().sum() - numeric.notna().sum())
                missing = int(source[column].isna().sum())
                non_finite = int((~numeric.dropna().map(math.isfinite)).sum())
                numeric_check_rows.append({
                    "fungicide": fungicide,
                    "source_column": column,
                    "n": int(numeric.notna().sum()),
                    "missing": missing,
                    "non_numeric": bad,
                    "non_finite": non_finite,
                    "minimum": float(numeric.min()) if numeric.notna().any() else None,
                    "maximum": float(numeric.max()) if numeric.notna().any() else None,
                })
                if column == "trait1" and (missing or bad or non_finite):
                    raise ValueError(
                        f"{fungicide} trait1 must contain finite numeric values for every isolate "
                        f"(missing={missing}, non-numeric={bad}, non-finite={non_finite})."
                    )

        canonical = pd.DataFrame({
            "study_id": STUDY_ID,
            "isolate_id": source["isolate_id"],
            "fungicide": fungicide,
            "trait_name": "trait1",
            "value": source["trait1"],
            "unit": "mM",
            "measurement_type": "author-reported isolate least-squares mean EC50",
            "replicate_id": "",
            "replicate_count": "3",
            "assay_protocol": "See study methods; model-based isolate LS mean across three experimental replications",
            "assay_batch": "",
            "standard_error": source["SE"] if "SE" in source else "",
            "degrees_of_freedom": source["df"] if "df" in source else "",
            # Keep unusual source headings intact rather than guessing their meaning.
            "source_lower_CL_u": source["lower.CL u"] if "lower.CL u" in source else "",
            "source_pper_CL": source["pper.CL"] if "pper.CL" in source else "",
            "resistance_threshold": "",
            "threshold_source": "",
            "source_reference": "Current study files; EC50 units and LS-mean definition confirmed by study author",
            "notes": "Trait values retained on supplied scale; no transformation applied.",
        })
        phenotype_frames.append(canonical)

    phenotypes = pd.concat(phenotype_frames, ignore_index=True)
    # Validate genotype coding and report the source matrix schema without expanding it.
    genotype_values = genotype.iloc[:, 3:].apply(lambda col: col.astype("string").str.strip().str.upper())
    observed_codes = sorted(set(genotype_values.stack().dropna().tolist()))
    unexpected_codes = sorted(set(observed_codes) - {"A", "C", "G", "T"})
    missing_genotypes = int(genotype_values.isna().sum().sum())
    sequence_accessions = sorted({
        str(value).rsplit("_", 1)[0].removeprefix("S")
        for value in genotype[variant_columns[0]].dropna()
        if "_" in str(value)
    })

    if unexpected_codes:
        raise ValueError(f"Unexpected genotype coding found: {unexpected_codes[:20]}")
    if missing_genotypes:
        raise ValueError(f"Missing genotype calls found ({missing_genotypes}); review before intake.")
    if len(variant_columns) != 3:
        raise ValueError(f"Expected three leading variant metadata columns, got {variant_columns}.")

    OUTPUT.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)
    manifest_path = OUTPUT / "isolate_manifest.csv"
    phenotype_path = OUTPUT / "phenotype_measurements.csv"
    genotype_manifest_path = OUTPUT / "genotype_input_manifest.json"
    intake_report_path = REPORTS / "study_intake_qc.md"
    numeric_checks_path = OUTPUT / "phenotype_numeric_checks.csv"
    study_config_path = OUTPUT / "study_config.json"
    manifest.to_csv(manifest_path, index=False)
    phenotypes.to_csv(phenotype_path, index=False)
    pd.DataFrame(numeric_check_rows).to_csv(numeric_checks_path, index=False)

    genotype_manifest = {
        "study_id": STUDY_ID,
        "pathogen": PATHOGEN,
        "source_file": "data/raw/genotype.csv",
        "representation": "variant rows by isolate columns",
        "variant_metadata_columns": variant_columns,
        "isolate_count": len(genotype_ids),
        "variant_count": int(genotype.shape[0]),
        "observed_genotype_codes": observed_codes,
        "missing_genotype_cells": missing_genotypes,
        "genome_build": "not specified by the study data repository; resolve before coordinate-based external harmonization",
        "coordinate_fields": {"chromosome": variant_columns[1], "position": variant_columns[2]},
        "sequence_accessions_observed_in_variant_ids": sequence_accessions,
        "coordinate_reference_hypothesis": "The observed HG970330-HG970335 accession family is consistent with the PH-1 RRes v4.0 reference, but the fungicide study repository does not explicitly confirm the build.",
        "allele_orientation": "source coding retained; reference/alternate orientation not inferred",
        "source_files_modified": False,
    }
    genotype_manifest_path.write_text(json.dumps(genotype_manifest, indent=2) + "\n", encoding="utf-8")

    study_config = {
        "study_id": STUDY_ID,
        "pathogen": PATHOGEN,
        "fungicides": sorted(PHENOTYPE_PATHS),
        "phenotype": {
            "trait_name": "trait1",
            "biological_interpretation": "isolate least-squares mean EC50",
            "unit": "mM",
            "replicate_count": 3,
            "transformation": "none; supplied scale retained",
            "resistance_threshold": None,
        },
        "genotype_manifest": "genotype_input_manifest.json",
        "phenotype_table": "phenotype_measurements.csv",
        "isolate_manifest": "isolate_manifest.csv",
        "interpretation_note": "Research intake only; this configuration does not define a susceptible/resistant cutoff.",
    }
    study_config_path.write_text(json.dumps(study_config, indent=2) + "\n", encoding="utf-8")

    lines = [
        "# Study intake and QC summary",
        "",
        f"Generated: {datetime.now().astimezone().isoformat(timespec='seconds')}",
        f"Study: `{STUDY_ID}` (*{PATHOGEN}*)",
        "",
        "## Intake results",
        "",
        f"- Isolate manifest rows: {len(manifest)}",
        f"- Genotype isolates: {len(genotype_ids)}; variants: {len(genotype):,}",
        f"- Genotype coding: {', '.join(observed_codes)}; missing calls: {missing_genotypes}",
        f"- Phenotype rows: {len(phenotypes)} (two fungicides × {len(manifest)} isolates)",
        "- Isolate IDs match exactly across genotype, metadata, and both phenotype files.",
        "- Phenotype `trait1` is represented as author-confirmed isolate least-squares mean EC50 in mM across three experimental replications.",
        "- Trait values are retained on their supplied scale; no log transform or resistance cutoff is applied.",
        "- Source fields `lower.CL u` and `pper.CL` are carried under source-labeled columns; their semantics are not inferred from the headings.",
        "- Genome build and reference/alternate allele orientation remain unspecified and must be resolved before coordinate-based external harmonization.",
        f"- Variant identifiers contain sequence accessions: {', '.join(sequence_accessions)}. These are consistent with the PH-1 RRes v4.0 accession family, but the study repository does not explicitly name the build.",
        "",
        "## Derived files",
        "",
        f"- `{manifest_path.relative_to(ROOT).as_posix()}`",
        f"- `{phenotype_path.relative_to(ROOT).as_posix()}`",
        f"- `{genotype_manifest_path.relative_to(ROOT).as_posix()}` (describes the original matrix without expanding or duplicating it)",
        f"- `{study_config_path.relative_to(ROOT).as_posix()}` (study-specific expectations used by the validator)",
        f"- `{numeric_checks_path.relative_to(ROOT).as_posix()}`",
        "",
        "Original source CSVs were read only and were not modified.",
    ]
    intake_report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("Study intake mapping completed.")
    print(f"Isolate manifest: {manifest_path}")
    print(f"Phenotype measurements: {phenotype_path}")
    print(f"Genotype input manifest: {genotype_manifest_path}")
    print(f"QC report: {intake_report_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"Study intake failed: {error}", file=sys.stderr)
        raise
