"""Compare a published external tebuconazole lead coordinate with this panel.

Run from the repository root in the phyloresist environment:
    python notebooks/11_external_lead_marker_check.py

This is a descriptive locus check in the current panel. It is not external
validation, and it does not replace a structure-corrected GWAS.
"""

from __future__ import annotations

import runpy
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
REFERENCE = ROOT / "data" / "reference"
GWAS_RESULTS = ROOT / "results" / "gwas" / "mrmlm" / "tebuconazole" / "1_Final result.csv"
OUTPUT_DIR = ROOT / "results" / "reports"
REPORT_PATH = OUTPUT_DIR / "north_dakota_lead_marker_crosscheck.md"
CSV_PATH = OUTPUT_DIR / "north_dakota_lead_marker_crosscheck.csv"

# North Dakota Phytopathology 2024 paper: SFG3_297186, annotated there to
# FGRAPHPH1_01G16325. The current panel uses scaffold-position marker IDs.
CURRENT_MARKER = "SHG970334.1_297186"
PUBLISHED_CHROMOSOME = "3"
PUBLISHED_POSITION = 297186
PUBLISHED_GENE_ID = "FGRAPHPH1_01G16325"
CURRENT_GENE_ID = "FGRAMPH1_01G16325"
POSSIBLE_OVERLAP_IDS = {"23262", "23269"}


def normalize_id(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip()
    return text[:-2] if text.endswith(".0") else text


def find_id_column(frame: pd.DataFrame) -> str:
    candidates = [
        col for col in frame.columns
        if any(token in str(col).casefold() for token in ("phenotype", "isolate", "sample", "strain", "id"))
    ]
    if not candidates:
        raise ValueError(f"Could not identify isolate ID column in {list(frame.columns)}")
    return candidates[0]


def main() -> int:
    genotype_path = RAW / "genotype.csv"
    phenotype_path = RAW / "phenotype_tebuconazole.csv"
    metadata_path = REFERENCE / "isolate_metadata.csv"
    required = [genotype_path, phenotype_path, metadata_path, GWAS_RESULTS]
    missing = [path for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing required input(s):\n" + "\n".join(map(str, missing)))

    genotype = pd.read_csv(genotype_path, dtype=str, keep_default_na=True, low_memory=False)
    if "rs#" not in genotype.columns:
        raise ValueError("Genotype marker ID column 'rs#' is missing.")
    marker_rows = genotype.loc[genotype["rs#"].astype(str).str.strip().eq(CURRENT_MARKER)]
    if len(marker_rows) != 1:
        raise ValueError(f"Expected one genotype row for {CURRENT_MARKER}; found {len(marker_rows)}.")

    calls = marker_rows.iloc[0, 3:].copy()
    calls.index = [normalize_id(value) for value in calls.index]
    calls.index.name = "isolate_id"
    calls = calls.astype("string").str.strip().str.upper()

    phenotype = pd.read_csv(phenotype_path, dtype=str, keep_default_na=True, low_memory=False)
    id_col = find_id_column(phenotype)
    trait_col = next((c for c in phenotype.columns if str(c).strip().casefold() == "trait1"), None)
    if trait_col is None:
        raise ValueError("trait1 column is missing from tebuconazole phenotype data.")
    phenotype_frame = pd.DataFrame(
        {
            "isolate_id": phenotype[id_col].map(normalize_id),
            "trait1_mM": pd.to_numeric(phenotype[trait_col], errors="coerce"),
        }
    )
    if phenotype_frame["isolate_id"].duplicated().any() or phenotype_frame["trait1_mM"].isna().any():
        raise ValueError("Duplicate isolate IDs or missing trait1 values found in phenotype data.")
    data = calls.rename("allele").reset_index().merge(phenotype_frame, on="isolate_id", validate="one_to_one")
    if len(data) != 152:
        raise ValueError(f"Expected 152 joined isolate records, found {len(data)}.")

    metadata = pd.read_csv(metadata_path, dtype=str, keep_default_na=False)
    metadata["isolate_id"] = metadata["isolate_id"].map(normalize_id)
    data = data.merge(metadata[["isolate_id", "publication", "alternate_id"]], on="isolate_id", how="left", validate="one_to_one")

    # Reuse the workbook reader that step 06 already documents and uses.
    helper_namespace = runpy.run_path(
        str(ROOT / "notebooks" / "06_post_gwas_candidate_review.py"),
        run_name="candidate_review_import",
    )
    annotations = helper_namespace["read_annotation_workbook"](REFERENCE / "full_annotation.xlsx")
    gene = annotations.loc[annotations["gene_id"].astype(str).eq(CURRENT_GENE_ID)]
    if len(gene) != 1:
        raise ValueError(f"Expected one annotation row for {CURRENT_GENE_ID}; found {len(gene)}.")
    gene_row = gene.iloc[0]
    marker_fields = CURRENT_MARKER.rsplit("_", 1)
    marker_scaffold = helper_namespace["normalize_scaffold"](marker_fields[0])
    marker_position = int(marker_fields[1])
    if marker_scaffold != str(gene_row["scaffold"]):
        raise ValueError(
            f"Marker scaffold {marker_scaffold} does not match gene scaffold {gene_row['scaffold']}."
        )
    if marker_position != PUBLISHED_POSITION:
        raise ValueError("Published and current marker coordinates are not equal.")
    marker_inside_gene = int(gene_row["gene_start"]) <= marker_position <= int(gene_row["gene_stop"])

    gwas = pd.read_csv(GWAS_RESULTS, dtype=str, keep_default_na=False)
    marker_reported = gwas["RS#"].astype(str).str.strip().eq(CURRENT_MARKER).any()

    output_rows: list[dict[str, object]] = []
    for group_label, excluded_ids in [
        ("all_current_isolates", set()),
        ("exclude_overlap_candidates", POSSIBLE_OVERLAP_IDS),
    ]:
        subset = data.loc[~data["isolate_id"].isin(excluded_ids)].copy()
        allele_groups = sorted(subset["allele"].dropna().unique().tolist())
        if len(allele_groups) != 2:
            raise ValueError(f"Expected two observed alleles in {group_label}; found {allele_groups}.")
        stats = subset.groupby("allele", observed=True)["trait1_mM"].agg(["count", "mean", "median", "std"])
        first, second = allele_groups
        test = ttest_ind(
            subset.loc[subset["allele"].eq(first), "trait1_mM"],
            subset.loc[subset["allele"].eq(second), "trait1_mM"],
            equal_var=False,
        )
        output_rows.append(
            {
                "analysis_set": group_label,
                "n_isolates": len(subset),
                "allele_1": first,
                "allele_1_n": int(stats.loc[first, "count"]),
                "allele_1_mean_trait1_mM": float(stats.loc[first, "mean"]),
                "allele_1_median_trait1_mM": float(stats.loc[first, "median"]),
                "allele_2": second,
                "allele_2_n": int(stats.loc[second, "count"]),
                "allele_2_mean_trait1_mM": float(stats.loc[second, "mean"]),
                "allele_2_median_trait1_mM": float(stats.loc[second, "median"]),
                "mean_difference_allele1_minus_allele2_mM": float(stats.loc[first, "mean"] - stats.loc[second, "mean"]),
                "descriptive_welch_p_unadjusted": float(test.pvalue),
                "reported_in_current_tebuconazole_mrmlm_results": bool(marker_reported),
            }
        )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    result = pd.DataFrame(output_rows)
    result.to_csv(CSV_PATH, index=False)

    lines = [
        "# North Dakota Tebuconazole Lead Marker Cross-Check",
        "",
        "The 2024 North Dakota GWAS reports marker `SFG3_297186` for tebuconazole and annotates it to `FGRAPHPH1_01G16325`. The current genotype file contains `SHG970334.1_297186`; its normalized scaffold and position match the current annotation workbook's `FGRAMPH1_01G16325` interval. The spelling of the gene ID differs between the published table and the local annotation, so this is recorded as a coordinate-supported candidate mapping, not an assumed gene-ID alias.",
        "",
        f"The variant is {'inside' if marker_inside_gene else 'outside'} the local annotated gene interval ({gene_row['scaffold']}:{int(gene_row['gene_start'])}-{int(gene_row['gene_stop'])}). Exact REF/ALT harmonization against the external assembly remains pending.",
        "",
        f"The marker {'was' if marker_reported else 'was not'} listed among the current panel's tebuconazole mrMLM final-result markers.",
        "",
        "## Descriptive allele summaries in the current panel",
        "",
        "| Set | n | Allele 1 (n; mean mM) | Allele 2 (n; mean mM) | Difference (mM) | Unadjusted Welch p |",
        "|---|---:|---|---|---:|---:|",
    ]
    for row in output_rows:
        lines.append(
            f"| {row['analysis_set']} | {row['n_isolates']} | "
            f"{row['allele_1']} ({row['allele_1_n']}; {row['allele_1_mean_trait1_mM']:.6f}) | "
            f"{row['allele_2']} ({row['allele_2_n']}; {row['allele_2_mean_trait1_mM']:.6f}) | "
            f"{row['mean_difference_allele1_minus_allele2_mM']:.6f} | "
            f"{row['descriptive_welch_p_unadjusted']:.4f} |"
        )
    lines.extend(
        [
            "",
            "The sensitivity check conservatively excludes IDs 23262 and 23269. Supplementary Table S1 confirms 23262 / Fg08-007 against Fg08_7 by county, year, and chemotype. The 23269 / Fg08-018 alias resembles Fg08_18, but the current county (Steele) conflicts with Table S1 (Griggs). Obtain the external crosswalk to resolve the second alias and check for any additional overlaps.",
            "",
            "These are unadjusted summaries in the current 152-isolate panel, not a replication test. They do not correct for population structure, relatedness, selection, or testing multiple markers. Do not compare effect directions until the published reference allele is mapped to the current C/T calls. The external paper's published association remains the evidence for its original cohort; this report only checks whether that coordinate is represented and how it behaves descriptively here.",
            "",
            "No independent data were available for this check. The next validation step still requires the North Dakota isolate crosswalk and genotype/phenotype files.",
            "",
            "Source: [North Dakota GWAS](https://pubmed.ncbi.nlm.nih.gov/38079375/). Full values: `north_dakota_lead_marker_crosscheck.csv`.",
            "",
        ]
    )
    REPORT_PATH.write_text("\n".join(lines), encoding="utf-8")
    print("External lead-marker cross-check completed.")
    print(f"Report: {REPORT_PATH}")
    print(f"CSV: {CSV_PATH}")
    print("Interpret descriptively; this is not independent validation.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"External lead-marker check failed: {error}", file=sys.stderr)
        raise
