"""Create a concise, evidence-separated candidate follow-up shortlist.

Run from the repository root in the phyloresist environment:
    python notebooks/10_candidate_shortlist.py

The script reshapes the existing candidate evidence matrix from step 06. It
does not recalculate GWAS statistics or treat proximity/LD as causality.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "results" / "gwas" / "mrmlm" / "review" / "candidate_evidence_matrix.csv"
OUT_DIR = ROOT / "results" / "gwas" / "mrmlm" / "review"
CSV_OUT = OUT_DIR / "candidate_shortlist.csv"
MD_OUT = OUT_DIR / "candidate_shortlist.md"


def truthy(value: object) -> bool:
    return str(value).strip().casefold() in {"true", "1", "yes"}


def number(value: object, default: float = 0.0) -> float:
    parsed = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    return default if pd.isna(parsed) else float(parsed)


def classify(row: pd.Series) -> tuple[str, str, str]:
    shared = truthy(row["exact_marker_shared_across_fungicides"])
    distance = number(row["distance_qtn_to_gene_interval_bp"], default=float("inf"))
    inside_count = int(number(row["genotyped_variants_inside_gene_with_pairwise_r2"]))
    spans_gene = bool(str(row["high_ld_marker_spanning_gene_interval"]).strip())

    if shared:
        return (
            "A — cross-fungicide follow-up",
            "The same reported marker occurs for both fungicides within the study's default gene window.",
            "Check the exact marker, local haplotypes, and replication in an independent isolate panel.",
        )
    if spans_gene:
        return (
            "B — LD-region follow-up",
            "A high-LD marker interval spans the annotated gene, but the nearest reported QTN is distant.",
            "Recheck coordinates and fine-map the region with denser genotypes before elevating the gene.",
        )
    if distance <= 3_000:
        local = (
            f" {inside_count} genotyped variants inside the annotated interval were available for local LD review."
            if inside_count
            else " No genotyped variants fell inside the annotated interval in this panel."
        )
        return (
            "B — local positional follow-up",
            f"The nearest reported QTN is within 3 kb of or overlaps the annotated gene interval.{local}",
            "Inspect the associated marker and local allele pattern; treat this as positional evidence only.",
        )
    return (
        "C — lower-priority / evidence gap",
        "The nearest reported QTN is outside the default 3 kb window and the current LD scan does not span the gene.",
        "Keep as a paper-mentioned candidate; revisit only if the final manuscript or denser data support it.",
    )


def main() -> int:
    if not INPUT.is_file():
        raise FileNotFoundError(f"Candidate evidence matrix not found: {INPUT}. Run notebook 06 first.")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    evidence = pd.read_csv(INPUT, keep_default_na=False)
    required = {
        "gene_id", "product_description", "nearest_qtn_fungicide", "nearest_qtn_marker",
        "nearest_qtn_methods", "distance_qtn_to_gene_interval_bp", "reported_markers_within_default_window",
        "methods_within_default_window", "exact_marker_shared_across_fungicides",
        "genotyped_variants_inside_gene_with_pairwise_r2", "maximum_qtn_r2_with_variant_inside_gene",
        "high_ld_marker_spanning_gene_interval", "spanning_marker_r2", "study_candidate_rationale",
    }
    missing = sorted(required - set(evidence.columns))
    if missing:
        raise ValueError(f"Candidate evidence matrix is missing required columns: {missing}")

    rows = []
    for _, source in evidence.iterrows():
        tier, evidence_summary, next_action = classify(source)
        rows.append(
            {
                "follow_up_tier": tier,
                "gene_id": source["gene_id"],
                "annotated_product": source["product_description"],
                "study_rationale": source["study_candidate_rationale"],
                "evidence_summary": evidence_summary,
                "nearest_qtn_fungicide": source["nearest_qtn_fungicide"],
                "nearest_qtn_marker": source["nearest_qtn_marker"],
                "nearest_qtn_methods": source["nearest_qtn_methods"],
                "distance_qtn_to_gene_interval_bp": source["distance_qtn_to_gene_interval_bp"],
                "reported_markers_within_default_window": source["reported_markers_within_default_window"],
                "methods_within_default_window": source["methods_within_default_window"],
                "same_marker_reported_for_both_fungicides": source["exact_marker_shared_across_fungicides"],
                "genotyped_variants_inside_gene": source["genotyped_variants_inside_gene_with_pairwise_r2"],
                "maximum_qtn_r2_with_variant_inside_gene": source["maximum_qtn_r2_with_variant_inside_gene"],
                "high_ld_marker_spanning_gene": source["high_ld_marker_spanning_gene_interval"],
                "spanning_marker_r2": source["spanning_marker_r2"],
                "suggested_next_check": next_action,
                "causal_status": "Unvalidated candidate; association, proximity, and LD do not establish causality.",
            }
        )

    tier_order = {
        "A — cross-fungicide follow-up": 0,
        "B — LD-region follow-up": 1,
        "B — local positional follow-up": 2,
        "C — lower-priority / evidence gap": 3,
    }
    shortlist = pd.DataFrame(rows)
    shortlist["_tier_order"] = shortlist["follow_up_tier"].map(tier_order)
    shortlist = shortlist.sort_values(
        ["_tier_order", "gene_id"], kind="stable"
    ).drop(columns="_tier_order")
    shortlist.to_csv(CSV_OUT, index=False)

    report = [
        "# Candidate Follow-up Shortlist",
        "",
        "This is a prioritization of candidates already discussed in the study materials, using the generated mrMLM, annotation, and exploratory local-LD evidence. Tiers indicate follow-up order, not confidence that a gene causes fungicide sensitivity.",
        "",
        "## Evidence source and status",
        "",
        "The manuscript record located for this review is the 2025 bioRxiv preprint (DOI: 10.1101/2025.09.25.678624) and its linked data repository. The attached supplement and current repository outputs were used for the local review. Reconcile marker counts, coordinates, and gene descriptions against the peer-reviewed version and its final supplement when available.",
        "",
        "The same genotype panel and 152 isolates underlie both fungicide analyses. A shared marker is therefore cross-trait evidence in this panel, not an independent replication. Markers inside a gene interval may be in LD with one another; the reported maximum r² is descriptive and is not a fine-mapping posterior probability.",
        "",
        "## Follow-up tiers",
        "",
        "- **Tier A:** same exact reported marker in both fungicide analyses within the default gene window.",
        "- **Tier B:** a localized gene-near QTN or a high-LD marker interval spanning a study-named gene; the distinction is shown per row.",
        "- **Tier C:** candidate named in study materials, but current default-window and local-LD evidence is weaker or incomplete.",
        "",
        "## Shortlist",
        "",
        "| Tier | Gene | Product | Evidence in this panel | Next check |",
        "|---|---|---|---|---|",
    ]
    for _, row in shortlist.iterrows():
        cells = [
            row["follow_up_tier"], row["gene_id"], row["annotated_product"],
            row["evidence_summary"], row["suggested_next_check"],
        ]
        report.append("| " + " | ".join(str(value).replace("|", "/").replace("\n", " ") for value in cells) + " |")
    report.extend(
        [
            "",
            "## Interpretation limits",
            "",
            "The genotype panel has 152 isolates and a single predominant source population. The geographic holdout analysis did not outperform the training-mean baseline on pooled RMSE, and it is not independent validation. Do not call these genes validated resistance markers or present this shortlist as a surveillance-ready assay.",
            "",
            "For the full statistics, inspect `candidate_evidence_matrix.csv`, `candidate_gene_ld_check.csv`, and `candidate_gene_ld_review.md`. The generated table is `candidate_shortlist.csv`.",
            "",
        ]
    )
    MD_OUT.write_text("\n".join(report), encoding="utf-8")
    print("Candidate shortlist completed.")
    print(f"CSV: {CSV_OUT}")
    print(f"Report: {MD_OUT}")
    print("Candidate ranking is for follow-up only; it does not imply causality or validation.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"Candidate shortlist failed: {error}", file=sys.stderr)
        raise
