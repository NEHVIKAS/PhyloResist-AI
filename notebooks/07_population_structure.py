"""Exploratory population-structure diagnostics for the starter genotype panel.

Run from the repository root with the phyloresist environment active:
    python notebooks/07_population_structure.py

The script treats each variant as a haploid biallelic marker (one of two
nucleotide calls per isolate), filters on minor-allele frequency, computes
standardized PCA scores, and summarizes pairwise allele sharing. The pairwise
similarity is descriptive only; it is not a formal kinship estimate. PCA does
not assign biological population labels. Interpret it alongside isolate
provenance, lineage, and collection metadata when those data are available.
"""

from __future__ import annotations

import math
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Any
from xml.sax.saxutils import escape

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
GENOTYPE_PATH = ROOT / "data" / "raw" / "genotype.csv"
ISOLATE_METADATA_PATH = ROOT / "data" / "reference" / "isolate_metadata.csv"
OUTPUT_DIR = ROOT / "results" / "population_structure"
METADATA_COLUMN_COUNT = 3
MINIMUM_MAF = 0.05
MAXIMUM_PCS = 10


def top_eigenpairs_symmetric(matrix: np.ndarray, count: int) -> tuple[np.ndarray, np.ndarray]:
    """Estimate leading eigenpairs using power iteration and rank-one deflation.

    This avoids vendor BLAS/LAPACK calls, which have failed in some Windows
    conda builds used by this project. The matrix is a small symmetric Gram
    matrix, and only its leading components are needed for exploratory PCA.
    """
    residual = matrix.copy()
    original = matrix
    rng = np.random.default_rng(2026)
    values: list[float] = []
    vectors: list[np.ndarray] = []

    for _ in range(count):
        vector = rng.standard_normal(matrix.shape[0])
        norm = math.sqrt(float(np.sum(vector * vector)))
        vector /= norm
        for _iteration in range(2_000):
            updated = np.sum(residual * vector[np.newaxis, :], axis=1)
            for previous in vectors:
                updated -= float(np.sum(updated * previous)) * previous
            updated_norm = math.sqrt(float(np.sum(updated * updated)))
            if updated_norm <= 1e-12:
                break
            updated /= updated_norm
            # Eigenvectors are sign-invariant; compare both possible signs.
            difference = min(
                math.sqrt(float(np.sum((updated - vector) ** 2))),
                math.sqrt(float(np.sum((updated + vector) ** 2))),
            )
            vector = updated
            if difference < 1e-9:
                break

        product = np.sum(original * vector[np.newaxis, :], axis=1)
        eigenvalue = float(np.sum(vector * product))
        if eigenvalue <= 1e-10:
            break
        values.append(eigenvalue)
        vectors.append(vector.copy())
        residual -= eigenvalue * np.outer(vector, vector)

    return np.asarray(values), np.column_stack(vectors) if vectors else np.empty((matrix.shape[0], 0))


def write_pca_svg(
    scores: np.ndarray,
    isolate_ids: list[str],
    explained_variance_ratio: np.ndarray,
    output_path: Path,
    metadata: pd.DataFrame,
) -> None:
    """Write a state-colored PC1/PC2 scatter plot using standard-library SVG."""
    if scores.shape[1] < 2:
        return
    width, height = 1180, 820
    left, top, plot_width, plot_height = 110, 90, 740, 540
    x_values, y_values = scores[:, 0], scores[:, 1]
    x_min, x_max = float(np.min(x_values)), float(np.max(x_values))
    y_min, y_max = float(np.min(y_values)), float(np.max(y_values))
    x_padding = (x_max - x_min) * 0.08 or 0.5
    y_padding = (y_max - y_min) * 0.08 or 0.5
    x_min, x_max = x_min - x_padding, x_max + x_padding
    y_min, y_max = y_min - y_padding, y_max + y_padding

    def x_pixel(value: float) -> float:
        return left + (value - x_min) / (x_max - x_min) * plot_width

    def y_pixel(value: float) -> float:
        return top + plot_height - (value - y_min) / (y_max - y_min) * plot_height

    state_values = metadata["collection_state"].fillna("").astype(str).tolist()
    population_values = metadata["population_label"].fillna("").astype(str).tolist()
    categories = sorted(set(state_values))
    palette = ["#356aa0", "#d97732", "#4c956c", "#9b5de5", "#d1495b", "#00798c", "#e0a458", "#6c757d"]
    colors = {category: palette[index % len(palette)] for index, category in enumerate(categories)}
    group_counts = pd.Series(state_values).value_counts().to_dict()

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">',
        '<title id="title">Exploratory genotype PCA, PC1 versus PC2</title>',
        f'<desc id="desc">Scatter plot of {len(isolate_ids)} isolates, colored by source state. PC1 explains {explained_variance_ratio[0] * 100:.2f}% and PC2 explains {explained_variance_ratio[1] * 100:.2f}% of variance. Hover for isolate, state, and population label.</desc>',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<text x="570" y="42" text-anchor="middle" font-family="Arial, sans-serif" font-size="24" font-weight="600" fill="#203040">Genotype PCA — exploratory, colored by collection state</text>',
    ]
    tick_count = 6
    for tick in range(tick_count):
        fraction = tick / (tick_count - 1)
        xv = x_min + fraction * (x_max - x_min)
        yv = y_min + fraction * (y_max - y_min)
        xp = x_pixel(xv)
        yp = y_pixel(yv)
        parts.append(f'<line x1="{xp:.2f}" y1="{top}" x2="{xp:.2f}" y2="{top + plot_height}" stroke="#e4e9ef"/>')
        parts.append(f'<line x1="{left}" y1="{yp:.2f}" x2="{left + plot_width}" y2="{yp:.2f}" stroke="#e4e9ef"/>')
        parts.append(f'<text x="{xp:.2f}" y="{top + plot_height + 24}" text-anchor="middle" font-family="Arial, sans-serif" font-size="13" fill="#44515e">{xv:.2f}</text>')
        parts.append(f'<text x="{left - 12}" y="{yp + 4:.2f}" text-anchor="end" font-family="Arial, sans-serif" font-size="13" fill="#44515e">{yv:.2f}</text>')
    parts.extend([
        f'<line x1="{left}" y1="{top + plot_height}" x2="{left + plot_width}" y2="{top + plot_height}" stroke="#526170" stroke-width="1.4"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_height}" stroke="#526170" stroke-width="1.4"/>',
    ])
    for isolate_id, state, population, x_value, y_value in zip(isolate_ids, state_values, population_values, x_values, y_values):
        tooltip = escape(f"Isolate {isolate_id} | State {state} | Population {population}")
        parts.append(
            f'<circle cx="{x_pixel(float(x_value)):.2f}" cy="{y_pixel(float(y_value)):.2f}" r="4.2" '
            f'fill="{colors[state]}" fill-opacity="0.82" stroke="#ffffff" stroke-width="0.8">'
            f'<title>{tooltip}</title></circle>'
        )
    legend_x, legend_y = 895, 120
    parts.append(f'<text x="{legend_x}" y="{legend_y}" font-family="Arial, sans-serif" font-size="17" font-weight="600" fill="#203040">Collection state</text>')
    for index, category in enumerate(categories):
        y = legend_y + 34 + index * 32
        parts.append(f'<circle cx="{legend_x + 7}" cy="{y - 5}" r="6" fill="{colors[category]}"/>')
        parts.append(f'<text x="{legend_x + 22}" y="{y}" font-family="Arial, sans-serif" font-size="14" fill="#44515e">{escape(category or "Unspecified")} (n={group_counts[category]})</text>')
    parts.extend([
        f'<text x="{left + plot_width / 2}" y="{height - 92}" text-anchor="middle" font-family="Arial, sans-serif" font-size="16" fill="#203040">PC1 ({explained_variance_ratio[0] * 100:.2f}% variance)</text>',
        f'<text x="28" y="{top + plot_height / 2}" text-anchor="middle" transform="rotate(-90 28 {top + plot_height / 2})" font-family="Arial, sans-serif" font-size="16" fill="#203040">PC2 ({explained_variance_ratio[1] * 100:.2f}% variance)</text>',
        '<text x="570" y="795" text-anchor="middle" font-family="Arial, sans-serif" font-size="13" fill="#526170">State is a geographic source label, not an independent population. Hover for isolate and source population label.</text>',
        '</svg>',
    ])
    output_path.write_text("\n".join(parts), encoding="utf-8")


def normalize_id(value: Any) -> str:
    """Keep isolate IDs as strings while reconciling integer-like CSV values."""
    text = str(value).strip()
    if re.fullmatch(r"\d+\.0+", text):
        text = text.split(".", 1)[0]
    return text


def main() -> int:
    if not GENOTYPE_PATH.is_file():
        raise FileNotFoundError(f"Genotype file not found: {GENOTYPE_PATH}")

    # Read all columns as text to preserve isolate IDs and nucleotide calls.
    genotype = pd.read_csv(GENOTYPE_PATH, dtype=str, keep_default_na=False, low_memory=False)
    if genotype.shape[1] <= METADATA_COLUMN_COUNT:
        raise ValueError("Expected three variant metadata columns followed by isolate columns.")

    isolate_ids = [normalize_id(value) for value in genotype.columns[METADATA_COLUMN_COUNT:]]
    if any(not value for value in isolate_ids) or len(isolate_ids) != len(set(isolate_ids)):
        raise ValueError("Isolate column IDs are blank or duplicated.")

    if not ISOLATE_METADATA_PATH.is_file():
        raise FileNotFoundError(
            f"Source-backed isolate metadata not found: {ISOLATE_METADATA_PATH}. "
            "Run notebooks/08_prepare_isolate_metadata.py first."
        )
    isolate_metadata = pd.read_csv(ISOLATE_METADATA_PATH, dtype=str, keep_default_na=False)
    if "isolate_id" not in isolate_metadata.columns or isolate_metadata["isolate_id"].duplicated().any():
        raise ValueError("The isolate metadata file must contain unique isolate_id values.")
    isolate_metadata["isolate_id"] = isolate_metadata["isolate_id"].map(normalize_id)
    if set(isolate_metadata["isolate_id"]) != set(isolate_ids):
        raise ValueError("Metadata isolate IDs do not exactly match the genotype columns.")
    isolate_metadata = isolate_metadata.set_index("isolate_id").loc[isolate_ids].reset_index()

    raw_calls = genotype.iloc[:, METADATA_COLUMN_COUNT:].copy()
    raw_calls.columns = isolate_ids
    raw_calls = raw_calls.apply(lambda column: column.astype("string").str.strip().str.upper())
    raw_calls = raw_calls.replace({"": pd.NA, "NA": pd.NA, "NAN": pd.NA})

    metadata = genotype.iloc[:, :METADATA_COLUMN_COUNT].copy()
    marker_name = str(genotype.columns[0])
    marker_ids = metadata.iloc[:, 0].astype(str).str.strip()

    # Convert only biallelic polymorphic variants to 0/1. The labels remain
    # nominal A/C/G/T alleles; the minor-allele indicator is a coding device.
    encoded_rows: list[np.ndarray] = []
    marker_records: list[dict[str, Any]] = []
    excluded_multiallelic = 0
    excluded_monomorphic = 0
    excluded_low_maf = 0
    excluded_no_calls = 0
    observed_unexpected: set[str] = set()

    for row_index, calls in raw_calls.iterrows():
        values = calls.to_numpy(dtype=object)
        valid_values = [str(value) for value in values if not pd.isna(value)]
        alleles = sorted(set(valid_values))
        observed_unexpected.update(set(alleles) - {"A", "C", "G", "T"})
        if not alleles:
            excluded_no_calls += 1
            continue
        if len(alleles) == 1:
            excluded_monomorphic += 1
            continue
        if len(alleles) != 2:
            excluded_multiallelic += 1
            continue

        counts = pd.Series(valid_values).value_counts()
        minor_allele = str(counts.idxmin())
        called_count = int(counts.sum())
        maf = float(counts.min() / called_count)
        if maf < MINIMUM_MAF:
            excluded_low_maf += 1
            continue

        # A missing call stays NaN here and is mean-imputed only for PCA.
        coded = np.array(
            [np.nan if pd.isna(value) else float(str(value) == minor_allele) for value in values],
            dtype=float,
        )
        encoded_rows.append(coded)
        marker_records.append(
            {
                "source_row_index": int(row_index),
                "marker_id": marker_ids.iloc[row_index],
                "chromosome": metadata.iloc[row_index, 1],
                "position": metadata.iloc[row_index, 2],
                "alleles": ";".join(alleles),
                "minor_allele": minor_allele,
                "maf": maf,
                "called_isolates": called_count,
            }
        )

    if observed_unexpected:
        raise ValueError(f"Unexpected genotype codes found: {sorted(observed_unexpected)[:20]}")
    if not encoded_rows:
        raise ValueError("No biallelic variants passed the configured MAF filter.")

    # Matrix orientation is isolates × variants. Impute any missing calls with
    # the variant allele frequency, then standardize each marker before PCA.
    allele_matrix = np.vstack(encoded_rows).T
    pairwise_allele_matrix = allele_matrix.copy()
    allele_frequencies = np.nanmean(allele_matrix, axis=0)
    missing_mask = np.isnan(allele_matrix)
    allele_matrix[missing_mask] = np.broadcast_to(allele_frequencies, allele_matrix.shape)[missing_mask]
    centered = allele_matrix - allele_frequencies
    scale = np.sqrt(allele_frequencies * (1.0 - allele_frequencies))
    if np.any(scale == 0):
        raise ValueError("A monomorphic marker passed the MAF filter unexpectedly.")
    standardized = centered / scale

    component_count = min(MAXIMUM_PCS, len(isolate_ids) - 1, standardized.shape[1])
    # Build the smaller isolate-by-isolate Gram matrix in bounded chunks.
    # np.einsum avoids the Windows BLAS path that failed in the local env.
    gram = np.zeros((len(isolate_ids), len(isolate_ids)), dtype=float)
    for start in range(0, standardized.shape[1], 256):
        block = standardized[:, start : start + 256]
        gram += np.einsum("ik,jk->ij", block, block, optimize=False)
    eigenvalues, eigenvectors = top_eigenpairs_symmetric(gram, component_count)
    scores = eigenvectors * np.sqrt(eigenvalues)
    total_variance = float(np.trace(gram))
    if total_variance <= 0:
        raise ValueError("PCA has zero total variance after filtering and standardization.")
    explained_variance_ratio = eigenvalues / total_variance
    pc_columns = [f"PC{index + 1}" for index in range(component_count)]

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    scores_frame = pd.DataFrame(scores, columns=pc_columns)
    scores_frame.insert(0, "isolate_id", isolate_ids)
    scores_frame.to_csv(OUTPUT_DIR / "pca_scores.csv", index=False)
    scores_with_metadata = scores_frame.merge(isolate_metadata, on="isolate_id", how="left", validate="one_to_one")
    scores_with_metadata.to_csv(OUTPUT_DIR / "pca_scores_with_metadata.csv", index=False)
    write_pca_svg(scores, isolate_ids, explained_variance_ratio, OUTPUT_DIR / "pca_pc1_pc2.svg", isolate_metadata)

    variance = pd.DataFrame(
        {
            "principal_component": pc_columns,
            "variance_explained_fraction": explained_variance_ratio,
            "variance_explained_percent": explained_variance_ratio * 100,
            "cumulative_variance_percent": np.cumsum(explained_variance_ratio) * 100,
        }
    )
    variance.to_csv(OUTPUT_DIR / "pca_variance.csv", index=False)

    # Pairwise allele sharing is the proportion of jointly called filtered
    # markers with matching allele states. It is not a kinship coefficient.
    pair_rows: list[dict[str, Any]] = []
    nearest: dict[str, tuple[str, float, int]] = {}
    for left_index in range(len(isolate_ids)):
        left = pairwise_allele_matrix[left_index]
        for right_index in range(left_index + 1, len(isolate_ids)):
            right = pairwise_allele_matrix[right_index]
            valid = np.isfinite(left) & np.isfinite(right)
            shared_count = int(valid.sum())
            if shared_count == 0:
                continue
            similarity = float(np.mean(left[valid] == right[valid]))
            left_id, right_id = isolate_ids[left_index], isolate_ids[right_index]
            pair_rows.append(
                {
                    "isolate_id_1": left_id,
                    "isolate_id_2": right_id,
                    "shared_variant_count": shared_count,
                    "allele_sharing_fraction": similarity,
                }
            )
            for sample_id, partner_id in ((left_id, right_id), (right_id, left_id)):
                current = nearest.get(sample_id)
                if current is None or similarity > current[1]:
                    nearest[sample_id] = (partner_id, similarity, shared_count)

    pair_frame = pd.DataFrame(pair_rows)
    pair_frame.to_csv(OUTPUT_DIR / "pairwise_allele_sharing.csv", index=False)
    nearest_rows = [
        {
            "isolate_id": sample_id,
            "nearest_genetic_neighbor_id": values[0],
            "allele_sharing_fraction": values[1],
            "shared_variant_count": values[2],
        }
        for sample_id, values in nearest.items()
    ]
    pd.DataFrame(nearest_rows).to_csv(OUTPUT_DIR / "nearest_genetic_neighbors.csv", index=False)
    pd.DataFrame(marker_records).to_csv(OUTPUT_DIR / "variants_used_for_pca.csv", index=False)

    lines = [
        "PHYLORESIST AI — EXPLORATORY POPULATION STRUCTURE REPORT",
        "=" * 72,
        f"Generated: {datetime.now().astimezone().isoformat(timespec='seconds')}",
        f"Genotype input: {GENOTYPE_PATH}",
        "",
        "INPUT AND FILTER SUMMARY",
        "-" * 72,
        f"Isolates: {len(isolate_ids)}",
        f"Variants in source: {len(genotype)}",
        f"Variants used in PCA (biallelic, MAF >= {MINIMUM_MAF:.2f}): {len(marker_records)}",
        f"Excluded monomorphic variants: {excluded_monomorphic}",
        f"Excluded multiallelic variants: {excluded_multiallelic}",
        f"Excluded variants below MAF threshold: {excluded_low_maf}",
        f"Excluded variants with no calls: {excluded_no_calls}",
        f"Missing genotype calls imputed for PCA: {int(missing_mask.sum())}",
        "",
        "PCA VARIANCE EXPLAINED",
        "-" * 72,
    ]
    lines.extend(
        f"{row.principal_component}: {row.variance_explained_percent:.2f}% "
        f"(cumulative {row.cumulative_variance_percent:.2f}%)"
        for row in variance.itertuples(index=False)
    )
    lines.extend(
        [
            "",
            "PAIRWISE ALLELE SHARING",
            "-" * 72,
            f"Isolate pairs summarized: {len(pair_frame)}",
            f"Pairwise allele-sharing fraction: mean={pair_frame['allele_sharing_fraction'].mean():.4f}, "
            f"min={pair_frame['allele_sharing_fraction'].min():.4f}, "
            f"max={pair_frame['allele_sharing_fraction'].max():.4f}",
            "",
            "INTERPRETATION LIMITS",
            "-" * 72,
            "PCA is an exploratory visualization/diagnostic, not a formal population assignment.",
            "Pairwise allele sharing is descriptive and is not a formal kinship coefficient.",
            "Points are colored by source collection state. State groups are geographic, not independent populations.",
            "The source population labels are highly imbalanced; do not use them to claim population-held-out validation.",
            "Some states have very small sample counts. State-held-out scores need per-state counts and cautious interpretation.",
            "These results do not validate prediction in an independent population.",
            "",
            "OUTPUT FILES",
            "-" * 72,
            "pca_pc1_pc2.svg; pca_scores.csv; pca_scores_with_metadata.csv; pca_variance.csv; pairwise_allele_sharing.csv;",
            "nearest_genetic_neighbors.csv; variants_used_for_pca.csv",
            "",
        ]
    )
    (OUTPUT_DIR / "population_structure_report.txt").write_text("\n".join(lines), encoding="utf-8")

    print("Exploratory population-structure diagnostics completed.")
    print(f"Isolates: {len(isolate_ids)}; variants used for PCA: {len(marker_records)}")
    print(f"Report: {OUTPUT_DIR / 'population_structure_report.txt'}")
    print(f"PCA scores: {OUTPUT_DIR / 'pca_scores.csv'}")
    print(f"PCA scores with metadata: {OUTPUT_DIR / 'pca_scores_with_metadata.csv'}")
    print(f"PCA plot: {OUTPUT_DIR / 'pca_pc1_pc2.svg'}")
    print("Pairwise allele sharing is descriptive, not a kinship estimate.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"Population-structure diagnostics failed: {error}", file=sys.stderr)
        raise
