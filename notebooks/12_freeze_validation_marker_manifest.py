"""Freeze the marker set for a future external validation analysis.

Run from the repository root in the phyloresist environment:
    python notebooks/12_freeze_validation_marker_manifest.py

Three cross-fungicide markers from the current shortlist are primary tests.
The 2024 North Dakota tebuconazole lead is listed separately as a published
benchmark. Current-panel allele summaries conservatively exclude one confirmed
overlap and one unresolved alias conflict; these are discovery descriptions,
not external validation results.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "data" / "raw"
REVIEW_DIR = ROOT / "results" / "gwas" / "mrmlm" / "review"
OUTPUT_CSV = REVIEW_DIR / "external_validation_marker_manifest.csv"
OUTPUT_MD = REVIEW_DIR / "external_validation_marker_manifest.md"
POSSIBLE_OVERLAP_IDS = {"23262", "23269"}
PUBLISHED_LEAD = {
    "marker_id": "SHG970334.1_297186",
    "gene_id": "FGRAPHPH1_01G16325 (local annotation: FGRAMPH1_01G16325)",
    "fungicide": "tebuconazole",
    "role": "published_external_benchmark",
}
FUNGICIDE_FILES = {
    "propiconazole": RAW_DIR / "phenotype_propiconazole.csv",
    "tebuconazole": RAW_DIR / "phenotype_tebuconazole.csv",
}
FINAL_RESULT_FILES = {
    fungicide: ROOT / "results" / "gwas" / "mrmlm" / fungicide / "1_Final result.csv"
    for fungicide in FUNGICIDE_FILES
}


def normalize_id(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    return re.sub(r"\.0+$", "", text)


def phenotype_table(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, dtype=str, keep_default_na=True, low_memory=False)
    id_candidates = [
        col for col in frame.columns
        if any(token in str(col).strip().casefold() for token in ("phenotype", "isolate", "sample", "strain", "id"))
    ]
    trait = next((col for col in frame.columns if str(col).strip().casefold() == "trait1"), None)
    if not id_candidates or trait is None:
        raise ValueError(f"Could not identify isolate ID and trait1 in {path.name}.")
    result = pd.DataFrame(
        {
            "isolate_id": frame[id_candidates[0]].map(normalize_id),
            "trait1_mM": pd.to_numeric(frame[trait], errors="coerce"),
        }
    )
    if result["isolate_id"].duplicated().any() or result["trait1_mM"].isna().any():
        raise ValueError(f"Duplicate IDs or missing trait values in {path.name}.")
    return result


def marker_methods(fungicide: str, marker: str) -> str:
    result = pd.read_csv(FINAL_RESULT_FILES[fungicide], dtype=str, keep_default_na=False)
    hits = result.loc[result["RS#"].astype(str).str.strip().eq(marker), "Method"].astype(str).str.strip()
    return ";".join(sorted(set(hits)))


def main() -> int:
    evidence_path = REVIEW_DIR / "candidate_evidence_matrix.csv"
    genotype_path = RAW_DIR / "genotype.csv"
    required = [evidence_path, genotype_path, *FUNGICIDE_FILES.values(), *FINAL_RESULT_FILES.values()]
    missing = [path for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Required inputs missing; run step 06 first:\n" + "\n".join(map(str, missing)))

    evidence = pd.read_csv(evidence_path, dtype=str, keep_default_na=False)
    genotype = pd.read_csv(genotype_path, dtype=str, keep_default_na=True, low_memory=False)
    marker_rows = genotype.set_index("rs#", drop=False)
    phenotype = {name: phenotype_table(path) for name, path in FUNGICIDE_FILES.items()}

    # Identify exact markers reported within the default gene window for both traits.
    primary_markers: dict[str, tuple[str, str]] = {}
    for _, gene in evidence.iterrows():
        if str(gene["exact_marker_shared_across_fungicides"]).strip().casefold() != "true":
            continue
        reported: dict[str, set[str]] = {name: set() for name in FUNGICIDE_FILES}
        for item in str(gene["reported_markers_within_default_window"]).split(";"):
            if "::" not in item:
                continue
            fungicide, marker = item.split("::", 1)
            if fungicide in reported:
                reported[fungicide].add(marker)
        shared = sorted(reported["propiconazole"] & reported["tebuconazole"])
        if not shared:
            raise ValueError(f"Flagged shared gene {gene['gene_id']} has no exact shared marker in its records.")
        primary_markers[shared[0]] = (gene["gene_id"], gene["product_description"])

    if len(primary_markers) != 3:
        raise ValueError(f"Expected three primary shared markers; found {len(primary_markers)}: {list(primary_markers)}")

    planned = [
        {
            "marker_id": marker,
            "gene_id": gene_id,
            "gene_description": description,
            "fungicide": "propiconazole;tebuconazole",
            "role": "primary_shared_marker_candidate",
        }
        for marker, (gene_id, description) in primary_markers.items()
    ]
    planned.append(PUBLISHED_LEAD | {"gene_description": "bHLH-related annotation in the external paper", "role": "published_external_benchmark"})

    overlap_excluded_genotypes = genotype.iloc[:, 3:].drop(columns=list(POSSIBLE_OVERLAP_IDS), errors="ignore")
    if len(overlap_excluded_genotypes.columns) != 150:
        raise ValueError("Expected 150 isolates after excluding the two known overlap candidates.")

    rows: list[dict[str, object]] = []
    for item in planned:
        marker = item["marker_id"]
        if marker not in marker_rows.index:
            raise ValueError(f"Frozen marker {marker} is absent from genotype.csv.")
        marker_record = marker_rows.loc[marker]
        if isinstance(marker_record, pd.DataFrame):
            raise ValueError(f"Duplicate genotype rows for frozen marker {marker}.")
        raw_calls = marker_record.iloc[3:].copy()
        raw_calls.index = [normalize_id(value) for value in raw_calls.index]
        calls = raw_calls.rename("allele").rename_axis("isolate_id").reset_index()
        calls["isolate_id"] = calls["isolate_id"].map(normalize_id)
        calls["allele"] = calls["allele"].astype("string").str.strip().str.upper()
        calls = calls.loc[~calls["isolate_id"].isin(POSSIBLE_OVERLAP_IDS)]
        allele_counts = calls["allele"].value_counts(dropna=True).to_dict()
        if len(allele_counts) != 2:
            raise ValueError(f"Expected two observed alleles at {marker}; found {allele_counts}.")
        alleles = sorted(allele_counts)

        for fungicide in item["fungicide"].split(";"):
            if fungicide not in phenotype:
                continue
            joined = calls.merge(phenotype[fungicide], on="isolate_id", how="inner", validate="one_to_one")
            means = joined.groupby("allele", observed=True)["trait1_mM"].mean()
            rows.append(
                {
                    **item,
                    "fungicide_trait": fungicide,
                    "marker_chromosome_in_current_file": marker_record["chrom"],
                    "marker_position_in_current_file": marker_record["pos"],
                    "current_marker_reported_by_mrmlm": bool(marker_methods(fungicide, marker)),
                    "mrmlm_methods_if_reported": marker_methods(fungicide, marker),
                    "analysis_n_after_overlap_candidate_exclusion": len(joined),
                    "overlap_candidate_ids_excluded": ";".join(sorted(POSSIBLE_OVERLAP_IDS)),
                    "allele_1": alleles[0],
                    "allele_1_n": int(allele_counts[alleles[0]]),
                    "allele_1_mean_trait1_mM_discovery_only": float(means[alleles[0]]),
                    "allele_2": alleles[1],
                    "allele_2_n": int(allele_counts[alleles[1]]),
                    "allele_2_mean_trait1_mM_discovery_only": float(means[alleles[1]]),
                    "mean_difference_allele1_minus_allele2_mM_discovery_only": float(means[alleles[0]] - means[alleles[1]]),
                    "external_allele_harmonization": "pending reference-build and REF/ALT crosswalk",
                    "external_test_predefined": True,
                }
            )

    REVIEW_DIR.mkdir(parents=True, exist_ok=True)
    manifest = pd.DataFrame(rows)
    manifest.to_csv(OUTPUT_CSV, index=False)

    lines = [
        "# Frozen External Validation Marker Manifest",
        "",
        "This manifest was generated before receipt of the North Dakota isolate-level genotype and phenotype files. It records the planned tebuconazole tests without selecting or changing markers based on external phenotypes.",
        "",
        "## Predefined marker set",
        "",
        "- **Three primary candidate markers:** exact markers reported for both fungicides within the current study's default gene windows, one each for FGRAMPH1_01G08383, FGRAMPH1_01G21231, and FGRAMPH1_01G25701.",
        "- **One benchmark marker:** the North Dakota paper's tebuconazole lead at chromosome 3 position 297186, represented as `SHG970334.1_297186` in the current genotype file. This is a prior published signal, not one of the three primary candidates.",
        "- **Trait tested externally:** tebuconazole only. The North Dakota panel does not include propiconazole.",
        "",
        "Current-panel allele counts and mean trait1 values below are included only to freeze discovery directions. They conservatively exclude ID 23262, which matches the supplement by alias, county, year, and chemotype, and ID 23269, whose alias resembles a supplement entry but whose county conflicts. The latter exclusion is provisional until the external isolate crosswalk resolves it. No external allele direction is assumed; genome build and REF/ALT harmonization remain required.",
        "",
        "| Role | Gene | Marker | Fungicide | Allele 1 (n; mean mM) | Allele 2 (n; mean mM) | Difference mM |",
        "|---|---|---|---|---|---|---:|",
    ]
    for row in rows:
        if row["fungicide_trait"] != "tebuconazole":
            continue
        lines.append(
            f"| {row['role']} | {row['gene_id']} | `{row['marker_id']}` | tebuconazole | "
            f"{row['allele_1']} (n={row['allele_1_n']}; {row['allele_1_mean_trait1_mM_discovery_only']:.6f}) | "
            f"{row['allele_2']} (n={row['allele_2_n']}; {row['allele_2_mean_trait1_mM_discovery_only']:.6f}) | "
            f"{row['mean_difference_allele1_minus_allele2_mM_discovery_only']:.6f} |"
        )
    lines.extend(
        [
            "",
            "Do not call the future test external validation until duplicate isolates are removed, exact marker coordinates and alleles are harmonized, and independent genotypes with matched phenotype values are analyzed. No association p-values are selected or reported by this manifest.",
            "",
            "Full marker-level fields are in `external_validation_marker_manifest.csv`.",
            "",
        ]
    )
    OUTPUT_MD.write_text("\n".join(lines), encoding="utf-8")
    print("External validation marker set frozen.")
    print(f"Manifest: {OUTPUT_CSV}")
    print(f"Guide: {OUTPUT_MD}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"Validation marker freeze failed: {error}", file=sys.stderr)
        raise
