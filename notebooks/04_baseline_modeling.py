"""Exploratory baseline prediction for the PhyloResist AI starter dataset.

Copy this file to the project's notebooks/ folder and run it from the project
root with the `phyloresist` conda environment active:

    python notebooks/04_baseline_modeling.py

The script evaluates whether genotype features predict the study's isolate-level
fungicide-sensitivity phenotype. It uses repeated outer cross-validation and
selects Ridge regularization only inside each training fold. This is an
exploratory baseline for a small dataset, not evidence of deployment readiness.

Phenotype provenance: Dhakal & Toomajian, bioRxiv 2025.09.25.678624. The paper
reports dose-response EC50 estimates and least-squares means across replicates.
The CSV's `trait1` is the isolate-level least-squares mean EC50 across three
replicates, stored in millimolar (mM), as confirmed by the study author. Some
replicate estimates were omitted when other dose-response parameters indicated
a poor curve fit. The script models trait1 as supplied and does not transform it.
"""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    from sklearn.compose import ColumnTransformer
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import Ridge
    from sklearn.metrics import mean_squared_error, r2_score
    from sklearn.model_selection import GridSearchCV, KFold, RepeatedKFold
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import OneHotEncoder, StandardScaler
except ImportError as exc:
    raise SystemExit(
        "This script needs pandas, numpy, and scikit-learn in the phyloresist "
        "environment. Install with: conda install pandas numpy scikit-learn"
    ) from exc


ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "data" / "raw"
OUTPUT_ROOT = ROOT / "results" / "baseline_modeling"
GENOTYPE_PATH = RAW_DIR / "genotype.csv"
PHENOTYPES = {
    "propiconazole": RAW_DIR / "phenotype_propiconazole.csv",
    "tebuconazole": RAW_DIR / "phenotype_tebuconazole.csv",
}
META_COLUMN_COUNT = 3
ALPHA_GRID = [0.1, 1.0, 10.0, 100.0]
OUTER_SPLITS = 5
OUTER_REPEATS = 2
RANDOM_SEED = 2026


def normalize_id(value: Any) -> str:
    """Make integer-like IDs comparable without losing leading zeros."""
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if re.fullmatch(r"\d+\.0+", text):
        text = text.split(".", 1)[0]
    return text


def load_genotypes(path: Path) -> pd.DataFrame:
    """Return isolates as rows and categorical variant calls as columns."""
    if not path.is_file():
        raise FileNotFoundError(f"Genotype file not found: {path}")
    variants = pd.read_csv(path, dtype=str, low_memory=False)
    if variants.shape[1] <= META_COLUMN_COUNT:
        raise ValueError("Genotype file has no isolate columns after the first three metadata columns.")

    isolate_columns = list(variants.columns[META_COLUMN_COUNT:])
    isolate_ids = [normalize_id(value) for value in isolate_columns]
    if any(not value for value in isolate_ids):
        raise ValueError("At least one genotype isolate column has an empty ID.")
    if len(isolate_ids) != len(set(isolate_ids)):
        raise ValueError("Duplicate genotype isolate IDs were found; resolve them before modelling.")

    # Make stable feature names from the variant metadata. The row index suffix
    # keeps names unique even if the source marker labels are repeated.
    features = []
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
    # Nucleotide calls remain categorical; no arbitrary A/C/G/T numeric order.
    matrix = matrix.apply(lambda column: column.astype("string").str.strip().str.upper())
    return matrix


def load_phenotype(path: Path) -> pd.DataFrame:
    """Read the isolate ID and trait1 only; other columns describe uncertainty."""
    if not path.is_file():
        raise FileNotFoundError(f"Phenotype file not found: {path}")
    phenotype = pd.read_csv(path, low_memory=False)
    if phenotype.shape[1] < 2:
        raise ValueError(f"Expected an isolate ID and trait1 in {path.name}.")

    id_candidates = [
        column for column in phenotype.columns
        if any(token in str(column).strip().casefold() for token in ("phenotype", "isolate", "sample", "strain", "id"))
    ]
    id_column = id_candidates[0] if id_candidates else phenotype.columns[0]
    if "trait1" not in phenotype.columns:
        raise ValueError(f"Column 'trait1' was not found in {path.name}; found {list(phenotype.columns)}")

    result = pd.DataFrame({
        "isolate_id": phenotype[id_column].map(normalize_id),
        "trait1": pd.to_numeric(phenotype["trait1"], errors="coerce"),
    })
    if result["isolate_id"].eq("").any():
        raise ValueError(f"Empty phenotype isolate IDs found in {path.name}.")
    if result["isolate_id"].duplicated().any():
        raise ValueError(f"Duplicate phenotype isolate IDs found in {path.name}.")
    if result["trait1"].isna().any() or not np.isfinite(result["trait1"]).all():
        raise ValueError(f"Missing or non-finite trait1 values found in {path.name}.")
    return result.set_index("isolate_id")


def make_pipeline() -> Pipeline:
    """Encode each variant as categorical alleles and fit regularized regression."""
    encode = ColumnTransformer(
        transformers=[
            (
                "alleles",
                Pipeline([
                    ("fill", SimpleImputer(strategy="most_frequent")),
                    ("one_hot", OneHotEncoder(handle_unknown="ignore")),
                ]),
                slice(0, None),
            )
        ],
        remainder="drop",
        verbose_feature_names_out=False,
    )
    return Pipeline([
        ("encode", encode),
        ("scale", StandardScaler(with_mean=False)),
        ("ridge", Ridge()),
    ])


def run_one_fungicide(name: str, path: Path, genotypes: pd.DataFrame) -> dict[str, Any]:
    phenotype = load_phenotype(path)
    shared = sorted(set(genotypes.index) & set(phenotype.index))
    genotype_only = sorted(set(genotypes.index) - set(phenotype.index))
    phenotype_only = sorted(set(phenotype.index) - set(genotypes.index))
    if len(shared) < 20:
        raise ValueError(f"Only {len(shared)} matching isolates for {name}; too few for this baseline.")

    # Keep only matched isolates and preserve a deterministic order.
    x = genotypes.loc[shared].reset_index(drop=True)
    y = phenotype.loc[shared, "trait1"].to_numpy(dtype=float)
    folds = min(OUTER_SPLITS, len(shared))
    if folds < 3:
        raise ValueError("At least three outer folds are required.")

    outer = RepeatedKFold(n_splits=folds, n_repeats=OUTER_REPEATS, random_state=RANDOM_SEED)
    predictions: list[dict[str, Any]] = []
    fold_metrics: list[dict[str, Any]] = []

    for split_number, (train_idx, test_idx) in enumerate(outer.split(x), start=1):
        total_splits = folds * OUTER_REPEATS
        print(f"[{name}] fitting outer fold {split_number}/{total_splits}...", flush=True)
        # The inner tuning loop sees only the outer training isolates.
        inner_splits = min(3, len(train_idx))
        inner = KFold(n_splits=inner_splits, shuffle=True, random_state=RANDOM_SEED + split_number)
        search = GridSearchCV(
            estimator=make_pipeline(),
            param_grid={"ridge__alpha": ALPHA_GRID},
            scoring="neg_root_mean_squared_error",
            cv=inner,
            # One worker avoids expensive process spawning on Windows/Miniforge.
            n_jobs=1,
            refit=True,
        )
        search.fit(x.iloc[train_idx], y[train_idx])
        predicted = search.predict(x.iloc[test_idx])
        train_mean = float(np.mean(y[train_idx]))
        baseline = np.full(len(test_idx), train_mean, dtype=float)
        fold_rmse = float(np.sqrt(mean_squared_error(y[test_idx], predicted)))
        baseline_rmse = float(np.sqrt(mean_squared_error(y[test_idx], baseline)))
        fold_metrics.append({
            "outer_split": split_number,
            "n_train": len(train_idx),
            "n_test": len(test_idx),
            "rmse_ridge_trait1_units": fold_rmse,
            "rmse_training_mean_trait1_units": baseline_rmse,
            "r2_ridge": float(r2_score(y[test_idx], predicted)),
            "best_alpha": float(search.best_params_["ridge__alpha"]),
        })
        repeat_number = (split_number - 1) // folds + 1
        fold_number = (split_number - 1) % folds + 1
        for index, actual, estimate, base_estimate in zip(test_idx, y[test_idx], predicted, baseline):
            predictions.append({
                "fungicide": name,
                "repeat": repeat_number,
                "fold": fold_number,
                "isolate_id": shared[index],
                "observed_trait1": float(actual),
                "predicted_trait1": float(estimate),
                "training_mean_baseline": float(base_estimate),
                "best_alpha": float(search.best_params_["ridge__alpha"]),
            })

    pred_df = pd.DataFrame(predictions)
    metric_df = pd.DataFrame(fold_metrics)
    out_dir = OUTPUT_ROOT / name
    out_dir.mkdir(parents=True, exist_ok=True)
    pred_df.to_csv(out_dir / "cross_validated_predictions.csv", index=False)
    metric_df.to_csv(out_dir / "fold_metrics.csv", index=False)

    observed = pred_df["observed_trait1"].to_numpy()
    predicted_all = pred_df["predicted_trait1"].to_numpy()
    baseline_all = pred_df["training_mean_baseline"].to_numpy()
    summary = {
        "fungicide": name,
        "n_isolates_matched": len(shared),
        "n_variants": genotypes.shape[1],
        "genotype_only_ids": len(genotype_only),
        "phenotype_only_ids": len(phenotype_only),
        "trait1_mean_original_scale": float(np.mean(y)),
        "trait1_sd_original_scale": float(np.std(y, ddof=1)),
        "trait1_unit_note": "mM, confirmed by the study author",
        "outer_splits": folds,
        "outer_repeats": OUTER_REPEATS,
        "ridge_cv_rmse_trait1_units": float(np.sqrt(mean_squared_error(observed, predicted_all))),
        "training_mean_cv_rmse_trait1_units": float(np.sqrt(mean_squared_error(observed, baseline_all))),
        "mean_fold_r2_ridge": float(metric_df["r2_ridge"].mean()),
        "sd_fold_r2_ridge": float(metric_df["r2_ridge"].std(ddof=1)),
        "mean_selected_alpha": float(metric_df["best_alpha"].mean()),
        "interpretation": "Exploratory internal cross-validation only; not independent external validation.",
    }
    pd.DataFrame([summary]).to_csv(out_dir / "summary.csv", index=False)
    (out_dir / "README.txt").write_text(
        "Exploratory genotype-to-phenotype baseline\n"
        "===========================================\n"
        f"Trait: {name} sensitivity, represented by the study's isolate-level EC50 least-squares mean.\n"
        "The original trait1 values are least-squares mean EC50 estimates across three replicates,\n"
        "in mM, confirmed by the study author. No log transformation is applied. Some replicate\n"
        "EC50 estimates were excluded before calculating means when other fit parameters indicated\n"
        "a poor dose-response fit.\n"
        "Ridge regression uses categorical nucleotide encoding. Alpha is selected within each outer\n"
        "training fold. Results are repeated internal cross-validation estimates, not external validation.\n"
        "A good score here does not establish transportability to other populations, years, or pathogens.\n",
        encoding="utf-8",
    )
    return summary


def main() -> int:
    if not GENOTYPE_PATH.is_file():
        raise FileNotFoundError(f"Could not find {GENOTYPE_PATH}. Run this script from the project copy.")
    genotypes = load_genotypes(GENOTYPE_PATH)
    summaries = []
    for fungicide, path in PHENOTYPES.items():
        summaries.append(run_one_fungicide(fungicide, path, genotypes))
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(summaries).to_csv(OUTPUT_ROOT / "baseline_summary.csv", index=False)

    print("Exploratory baseline modeling finished.")
    print(f"Summary: {OUTPUT_ROOT / 'baseline_summary.csv'}")
    print("Check each fungicide folder for fold metrics and isolate-level predictions.")
    print("Trait1 values were used on their original scale; no log transformation was applied.")
    print("Treat results as exploratory internal validation, not deployment performance.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"Baseline modeling failed: {error}", file=sys.stderr)
        raise
