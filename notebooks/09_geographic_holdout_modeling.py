"""Exploratory leave-one-state-out prediction for the starter dataset.

Run from the repository root in the phyloresist environment:
    python notebooks/09_geographic_holdout_modeling.py

Each outer test fold holds out every isolate from one source collection state.
Ridge alpha is selected with group-aware inner CV using only the remaining
states. This tests geographic transfer within the supplied study panel; it is
not independent external validation or population-held-out validation.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    from sklearn.model_selection import GroupKFold
except ImportError as exc:
    raise SystemExit(
        "This script needs pandas, numpy, and scikit-learn in the phyloresist "
        "environment. Install with: conda install pandas numpy scikit-learn"
    ) from exc


ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "data" / "raw"
METADATA_PATH = ROOT / "data" / "reference" / "isolate_metadata.csv"
OUTPUT_ROOT = ROOT / "results" / "geographic_holdout_validation"
GENOTYPE_PATH = RAW_DIR / "genotype.csv"
PHENOTYPES = {
    "propiconazole": RAW_DIR / "phenotype_propiconazole.csv",
    "tebuconazole": RAW_DIR / "phenotype_tebuconazole.csv",
}
ALPHA_GRID = [0.1, 1.0, 10.0, 100.0]
META_COLUMN_COUNT = 3


def normalize_id(value: Any) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    return text.split(".", 1)[0] if re.fullmatch(r"\d+\.0+", text) else text


def load_genotypes(path: Path) -> pd.DataFrame:
    """Return isolates as rows and categorical variant calls as columns."""
    if not path.is_file():
        raise FileNotFoundError(f"Genotype file not found: {path}")
    variants = pd.read_csv(path, dtype=str, keep_default_na=True, low_memory=False)
    if variants.shape[1] <= META_COLUMN_COUNT:
        raise ValueError("Expected three variant metadata columns plus isolate columns.")
    isolate_ids = [normalize_id(value) for value in variants.columns[META_COLUMN_COUNT:]]
    if any(not value for value in isolate_ids) or len(isolate_ids) != len(set(isolate_ids)):
        raise ValueError("Genotype isolate IDs are blank or duplicated.")

    features: list[str] = []
    for row_number, row in variants.iterrows():
        marker = normalize_id(row.iloc[0]) or "marker"
        chromosome = normalize_id(row.iloc[1]) or "unknown_chr"
        position = normalize_id(row.iloc[2]) or "unknown_pos"
        features.append(f"variant_{row_number + 1}_{chromosome}_{position}_{marker}")
    matrix = variants.iloc[:, META_COLUMN_COUNT:].copy()
    matrix.columns = isolate_ids
    matrix.index = features
    matrix = matrix.transpose()
    matrix.index.name = "isolate_id"
    return matrix.apply(lambda column: column.astype("string").str.strip().str.upper())


def load_phenotype(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Phenotype file not found: {path}")
    frame = pd.read_csv(path, dtype=str, keep_default_na=True, low_memory=False)
    id_candidates = [
        column for column in frame.columns
        if re.search(r"(^|[^a-z])(id|isolate|sample|strain)([^a-z]|$)", str(column).strip().casefold())
    ]
    phenotype_labels = [column for column in frame.columns if "phenotype" in str(column).casefold()]
    id_column = (id_candidates or phenotype_labels or list(frame.columns))[0]
    trait_column = next((column for column in frame.columns if str(column).strip().casefold() == "trait1"), None)
    if trait_column is None:
        raise ValueError(f"Column trait1 is missing from {path.name}.")
    result = pd.DataFrame(
        {
            "isolate_id": frame[id_column].map(normalize_id),
            "trait1": pd.to_numeric(frame[trait_column], errors="coerce"),
        }
    )
    if result["isolate_id"].eq("").any() or result["isolate_id"].duplicated().any():
        raise ValueError(f"Blank or duplicate isolate IDs in {path.name}.")
    if result["trait1"].isna().any() or not np.isfinite(result["trait1"].to_numpy(dtype=float)).all():
        raise ValueError(f"Missing or non-finite trait1 in {path.name}.")
    return result.set_index("isolate_id")


def encode_biallelic_genotypes(genotypes: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    """Encode one nominal allele indicator per biallelic variant.

    A/C/G/T labels have no numeric order. The 0/1 indicator is only a feature
    coding choice; swapping the two alleles leaves a linear model with an
    intercept and symmetric regularization unchanged.
    """
    values = genotypes.to_numpy(dtype=object)
    encoded = np.full(values.shape, np.nan, dtype=float)
    dropped: list[str] = []
    for column_index, feature in enumerate(genotypes.columns):
        calls = values[:, column_index]
        called = [str(value).strip().upper() for value in calls if not pd.isna(value) and str(value).strip()]
        alleles = sorted(set(called))
        if len(alleles) > 2:
            raise ValueError(f"Variant {feature} has more than two observed alleles: {alleles}")
        if len(alleles) == 0:
            dropped.append(str(feature))
            continue
        if len(alleles) == 1:
            encoded[:, column_index] = 0.0
            continue
        minor_code = alleles[1]
        for row_index, value in enumerate(calls):
            if pd.isna(value) or not str(value).strip():
                continue
            call = str(value).strip().upper()
            if call not in alleles:
                raise ValueError(f"Unexpected call {call!r} at variant {feature}")
            encoded[row_index, column_index] = float(call == minor_code)
    keep = [index for index, feature in enumerate(genotypes.columns) if str(feature) not in dropped]
    return encoded[:, keep], [str(genotypes.columns[index]) for index in keep]


def prepare_ridge_fold(
    x_train: np.ndarray,
    y_train: np.ndarray,
    x_valid: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Build dual Ridge matrices without the local Windows BLAS backend."""
    with np.errstate(invalid="ignore"):
        means = np.nanmean(x_train, axis=0)
    means = np.where(np.isfinite(means), means, 0.0)
    train_filled = np.where(np.isfinite(x_train), x_train, means)
    valid_filled = np.where(np.isfinite(x_valid), x_valid, means)
    scales = np.sqrt(np.mean((train_filled - means) ** 2, axis=0))
    scales = np.where(scales > 1e-12, scales, 1.0)
    train_scaled = (train_filled - means) / scales
    valid_scaled = (valid_filled - means) / scales

    # einsum and row-wise elementwise sums avoid safe_sparse_dot/BLAS, which
    # raised a Windows native error in the supplied phyloresist conda env.
    gram = np.einsum("ik,jk->ij", train_scaled, train_scaled, optimize=False)
    cross_gram = np.einsum("ik,jk->ij", valid_scaled, train_scaled, optimize=False)
    y_mean = float(np.mean(y_train))
    return gram, cross_gram, y_train - y_mean, y_mean


def solve_dual_ridge(gram: np.ndarray, centered_target: np.ndarray, alpha: float) -> np.ndarray:
    """Solve (X X' + alpha I)a = centered_y with conjugate gradients."""
    system = gram.copy()
    diagonal = np.diag_indices_from(system)
    system[diagonal] += alpha
    solution = np.zeros_like(centered_target, dtype=float)
    residual = centered_target.copy()
    direction = residual.copy()
    residual_norm_sq = float(np.sum(residual * residual))
    initial_norm = float(np.sqrt(residual_norm_sq))
    if initial_norm <= 1e-12:
        return solution

    tolerance = 1e-9 * max(1.0, initial_norm)
    for _iteration in range(max(100, min(3000, 20 * len(centered_target)))):
        product = np.sum(system * direction[np.newaxis, :], axis=1)
        denominator = float(np.sum(direction * product))
        if denominator <= 1e-20:
            break
        step = residual_norm_sq / denominator
        solution += step * direction
        residual -= step * product
        new_norm_sq = float(np.sum(residual * residual))
        if np.sqrt(new_norm_sq) <= tolerance:
            return solution
        direction = residual + (new_norm_sq / residual_norm_sq) * direction
        residual_norm_sq = new_norm_sq

    remaining = float(np.sqrt(np.sum(residual * residual)))
    if remaining > 1e-6 * max(1.0, initial_norm):
        raise RuntimeError(f"Dual Ridge solver did not converge (residual={remaining:.3g}).")
    return solution


def predict_dual_ridge(
    x_train: np.ndarray,
    y_train: np.ndarray,
    x_valid: np.ndarray,
    alpha: float,
) -> np.ndarray:
    gram, cross_gram, centered_target, y_mean = prepare_ridge_fold(x_train, y_train, x_valid)
    dual_weights = solve_dual_ridge(gram, centered_target, alpha)
    return y_mean + np.sum(cross_gram * dual_weights[np.newaxis, :], axis=1)


def rmse(observed: np.ndarray, predicted: np.ndarray) -> float:
    return float(np.sqrt(np.mean((observed - predicted) ** 2)))


def mae(observed: np.ndarray, predicted: np.ndarray) -> float:
    return float(np.mean(np.abs(observed - predicted)))


def r2_or_blank(observed: np.ndarray, predicted: np.ndarray) -> float | str:
    if len(observed) < 2:
        return ""
    total = float(np.sum((observed - np.mean(observed)) ** 2))
    return "" if total <= 1e-12 else float(1.0 - np.sum((observed - predicted) ** 2) / total)


def run_trait(
    fungicide: str,
    phenotype_path: Path,
    genotypes: pd.DataFrame,
    metadata: pd.DataFrame,
) -> dict[str, Any]:
    phenotype = load_phenotype(phenotype_path)
    genotype_ids = set(genotypes.index)
    phenotype_ids = set(phenotype.index)
    shared = sorted(genotype_ids & phenotype_ids)
    if genotype_ids != phenotype_ids:
        raise ValueError(
            f"Genotype/phenotype isolate IDs differ for {fungicide}: "
            f"genotype-only={sorted(genotype_ids - phenotype_ids)[:10]}, "
            f"phenotype-only={sorted(phenotype_ids - genotype_ids)[:10]}"
        )
    if set(metadata.index) != genotype_ids:
        raise ValueError("Metadata IDs do not exactly match the genotype isolate IDs.")

    x = genotypes.loc[shared].reset_index(drop=True)
    y = phenotype.loc[shared, "trait1"].to_numpy(dtype=float)
    states = metadata.loc[shared, "collection_state"].astype(str).to_numpy()
    if np.any(pd.Series(states).str.strip().eq("")):
        raise ValueError("Some isolates have no source collection state; cannot define geographic groups.")
    unique_states = sorted(set(states))
    if len(unique_states) < 3:
        raise ValueError("At least three collection states are needed for grouped validation.")

    x_encoded, encoded_features = encode_biallelic_genotypes(x)

    state_metrics: list[dict[str, Any]] = []
    prediction_rows: list[dict[str, Any]] = []
    split_count = len(unique_states)

    for split_number, held_state in enumerate(unique_states, start=1):
        test_idx = np.flatnonzero(states == held_state)
        train_idx = np.flatnonzero(states != held_state)
        train_states = states[train_idx]
        inner_group_count = len(set(train_states))
        inner_splits = min(3, inner_group_count)
        inner = GroupKFold(n_splits=inner_splits)
        print(
            f"[{fungicide}] holding out {held_state} "
            f"({len(test_idx)} isolates; state fold {split_number}/{split_count})...",
            flush=True,
        )
        alpha_scores: dict[float, list[float]] = {alpha: [] for alpha in ALPHA_GRID}
        inner_x = x_encoded[train_idx]
        inner_y = y[train_idx]
        for inner_train, inner_valid in inner.split(inner_x, inner_y, groups=train_states):
            for alpha in ALPHA_GRID:
                inner_prediction = predict_dual_ridge(
                    inner_x[inner_train], inner_y[inner_train], inner_x[inner_valid], alpha
                )
                alpha_scores[alpha].append(rmse(inner_y[inner_valid], inner_prediction))
        selected_alpha = min(ALPHA_GRID, key=lambda alpha: float(np.mean(alpha_scores[alpha])))
        predicted = predict_dual_ridge(x_encoded[train_idx], y[train_idx], x_encoded[test_idx], selected_alpha)
        train_mean = float(np.mean(y[train_idx]))
        baseline = np.full(len(test_idx), train_mean, dtype=float)
        state_metrics.append(
            {
                "held_out_state": held_state,
                "n_train": len(train_idx),
                "n_test": len(test_idx),
                "n_training_states": inner_group_count,
                "ridge_rmse_trait1_mM": rmse(y[test_idx], predicted),
                "training_mean_rmse_trait1_mM": rmse(y[test_idx], baseline),
                "ridge_r2": r2_or_blank(y[test_idx], predicted),
                "best_alpha_inner_state_cv": float(selected_alpha),
            }
        )
        for index, observed, estimate, base_estimate in zip(test_idx, y[test_idx], predicted, baseline):
            sample_id = shared[index]
            prediction_rows.append(
                {
                    "fungicide": fungicide,
                    "isolate_id": sample_id,
                    "collection_state": held_state,
                    "population_label": metadata.loc[sample_id, "population_label"],
                    "observed_trait1_mM": float(observed),
                    "predicted_trait1_mM": float(estimate),
                    "training_mean_baseline_mM": float(base_estimate),
                    "best_alpha_inner_state_cv": float(selected_alpha),
                }
            )

    pred_frame = pd.DataFrame(prediction_rows)
    state_frame = pd.DataFrame(state_metrics)
    observed_all = pred_frame["observed_trait1_mM"].to_numpy(dtype=float)
    predicted_all = pred_frame["predicted_trait1_mM"].to_numpy(dtype=float)
    baseline_all = pred_frame["training_mean_baseline_mM"].to_numpy(dtype=float)
    summary = {
        "fungicide": fungicide,
        "n_isolates": len(shared),
        "n_variants": len(encoded_features),
        "n_states": len(unique_states),
        "grouping_field": "collection_state",
        "outer_scheme": "Leave one source state out",
        "inner_scheme": "GroupKFold by source state (up to 3 folds)",
        "ridge_pooled_rmse_trait1_mM": rmse(observed_all, predicted_all),
        "training_mean_pooled_rmse_trait1_mM": rmse(observed_all, baseline_all),
        "ridge_pooled_mae_trait1_mM": mae(observed_all, predicted_all),
        "ridge_pooled_r2": r2_or_blank(observed_all, predicted_all),
        "mean_state_rmse_trait1_mM": float(state_frame["ridge_rmse_trait1_mM"].mean()),
        "mean_state_baseline_rmse_trait1_mM": float(state_frame["training_mean_rmse_trait1_mM"].mean()),
        "trait_interpretation": "trait1 is isolate least-squares mean EC50 in mM; original scale, no log transform",
        "interpretation": "Exploratory within-study geographic holdout; not independent external validation or population-held-out validation.",
    }

    out_dir = OUTPUT_ROOT / fungicide
    out_dir.mkdir(parents=True, exist_ok=True)
    state_frame.to_csv(out_dir / "state_fold_metrics.csv", index=False)
    pred_frame.to_csv(out_dir / "held_out_state_predictions.csv", index=False)
    pd.DataFrame([summary]).to_csv(out_dir / "summary.csv", index=False)
    (out_dir / "README.txt").write_text(
        "Exploratory within-study geographic holdout prediction\n"
        "======================================================\n"
        "Each outer fold holds out all isolates from one collection state.\n"
        "Ridge alpha is tuned with inner group-aware folds using only the training states.\n"
        "The response is trait1 on its original mM scale; no log transformation is used.\n"
        "This is an internal geographic-transfer check, not external validation, and\n"
        "collection states are not equivalent to independent pathogen populations.\n"
        "Some held-out states have very small sample counts; inspect state_fold_metrics.csv.\n",
        encoding="utf-8",
    )
    return summary


def main() -> int:
    if not METADATA_PATH.is_file():
        raise FileNotFoundError(f"Isolate metadata not found: {METADATA_PATH}. Run notebook 08 first.")
    genotypes = load_genotypes(GENOTYPE_PATH)
    metadata = pd.read_csv(METADATA_PATH, dtype=str, keep_default_na=False)
    if metadata["isolate_id"].duplicated().any():
        raise ValueError("Duplicate isolate_id values found in isolate_metadata.csv.")
    metadata["isolate_id"] = metadata["isolate_id"].map(normalize_id)
    metadata = metadata.set_index("isolate_id")
    if "collection_state" not in metadata or "population_label" not in metadata:
        raise ValueError("Metadata must contain collection_state and population_label.")

    summaries = [run_trait(name, path, genotypes, metadata) for name, path in PHENOTYPES.items()]
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(summaries).to_csv(OUTPUT_ROOT / "geographic_holdout_summary.csv", index=False)

    state_counts = metadata["collection_state"].value_counts().sort_index()
    population_counts = metadata["population_label"].value_counts().to_dict()
    summary_by_name = {row["fungicide"]: row for row in summaries}
    lines = [
        "# Geographic Holdout Validation",
        "",
        "Exploratory leave-one-state-out validation was run for the two fungicides. Each outer test fold contains all isolates from one collection state; inner alpha selection also holds out state groups.",
        "",
        "## Group coverage",
        "",
        *[f"- {state}: n={count}" for state, count in state_counts.items()],
        "",
        f"Source population labels: {population_counts}",
        "",
        "Population labels are highly imbalanced (the panel is predominantly NA1), so geographic state is used as the only viable documented group field here. Geographic holdout is an internal transfer scenario, not an independent population validation. Tiny held-out state groups make individual fold metrics unstable.",
        "",
        "Both traits are analyzed as supplied: least-squares mean EC50 in mM, without log transformation. See each fungicide subfolder for state-fold metrics and isolate predictions.",
        "",
        "## Pooled results",
        "",
    ]
    for fungicide in PHENOTYPES:
        result = summary_by_name[fungicide]
        lines.append(
            f"- {fungicide}: Ridge RMSE {result['ridge_pooled_rmse_trait1_mM']:.6f} mM; "
            f"training-mean baseline RMSE {result['training_mean_pooled_rmse_trait1_mM']:.6f} mM; "
            f"Ridge R² {result['ridge_pooled_r2']:.3f}."
        )
    lines.extend(
        [
            "",
            "In this run, pooled Ridge RMSE is slightly higher than the training-mean baseline for both fungicides, so the genotype model did not improve pooled prediction error. The selected alpha was 100 in every outer fold (the largest value tested); this suggests regularization tuning reached the current grid boundary and should be revisited only if a later analysis is justified. These exploratory results do not support a surveillance performance claim.",
            "",
        ]
    )
    (OUTPUT_ROOT / "geographic_holdout_report.md").write_text("\n".join(lines), encoding="utf-8")

    print("Geographic holdout modeling completed.")
    print(f"Summary: {OUTPUT_ROOT / 'geographic_holdout_summary.csv'}")
    print(f"Report: {OUTPUT_ROOT / 'geographic_holdout_report.md'}")
    print("Interpret as exploratory geographic transfer within this study, not external validation.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"Geographic holdout modeling failed: {error}", file=sys.stderr)
        raise
