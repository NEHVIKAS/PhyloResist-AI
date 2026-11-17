"""Quality-control utilities for replicate-level EC50 measurements."""

import numpy as np
import pandas as pd


def summarize_replicates(
    data: pd.DataFrame,
    isolate_col: str = "isolate_id",
    fungicide_col: str = "fungicide",
    replicate_col: str = "replicate",
    value_col: str = "ec50",
    expected_replicates: int = 3,
    max_cv_percent: float | None = None,
) -> pd.DataFrame:
    """
    Summarize raw EC50 replicates by isolate and fungicide.

    Invalid and missing observations are reported, not silently removed.
    Only valid positive finite values contribute to descriptive statistics.

    Groups with fewer or more observations than expected require review.
    Additional observations are flagged, not automatically rejected.
    """

    if not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a pandas DataFrame.")

    if data.empty:
        raise ValueError("data must contain at least one observation.")

    required_columns = [
        isolate_col,
        fungicide_col,
        replicate_col,
        value_col,
    ]

    if len(set(required_columns)) != len(required_columns):
        raise ValueError("Column names must be distinct.")

    missing_columns = [
        column for column in required_columns
        if column not in data.columns
    ]

    if missing_columns:
        raise ValueError(
            f"Missing required columns: {missing_columns}"
        )

    if (
        isinstance(expected_replicates, bool)
        or not isinstance(expected_replicates, int)
        or expected_replicates < 1
    ):
        raise ValueError(
            "expected_replicates must be a positive integer."
        )

    if max_cv_percent is not None:
        if (
            isinstance(max_cv_percent, bool)
            or not isinstance(max_cv_percent, (int, float))
            or not np.isfinite(max_cv_percent)
            or max_cv_percent < 0
        ):
            raise ValueError(
                "max_cv_percent must be None or a "
                "non-negative finite number."
            )

    working = data[required_columns].copy()

    grouping_columns = [isolate_col, fungicide_col]

    if working[grouping_columns].isna().any().any():
        raise ValueError(
            "Isolate and fungicide identifiers cannot be missing."
        )

    if working[replicate_col].isna().any():
        raise ValueError(
            "Replicate identifiers cannot be missing."
        )

    duplicate_mask = working.duplicated(
        subset=[
            isolate_col,
            fungicide_col,
            replicate_col,
        ],
        keep=False,
    )

    if duplicate_mask.any():
        duplicates = working.loc[
            duplicate_mask,
            [
                isolate_col,
                fungicide_col,
                replicate_col,
            ],
        ]

        raise ValueError(
            "Duplicate replicate identifiers found within "
            "isolate/fungicide groups:\n"
            f"{duplicates.to_string(index=False)}"
        )

    numeric_values = pd.to_numeric(
        working[value_col],
        errors="coerce",
    )

    missing_mask = working[value_col].isna()

    invalid_mask = (
        numeric_values.isna()
        & ~missing_mask
    )

    non_finite_mask = (
        numeric_values.notna()
        & ~np.isfinite(numeric_values)
    )

    non_positive_mask = (
        numeric_values.notna()
        & np.isfinite(numeric_values)
        & (numeric_values <= 0)
    )

    valid_mask = (
        numeric_values.notna()
        & np.isfinite(numeric_values)
        & (numeric_values > 0)
    )

    working["_numeric_ec50"] = numeric_values
    working["_is_missing"] = missing_mask
    working["_is_invalid"] = (
        invalid_mask
        | non_finite_mask
        | non_positive_mask
    )
    working["_is_valid"] = valid_mask

    results = []

    for group_key, group in working.groupby(
        grouping_columns,
        sort=True,
        dropna=False,
    ):
        isolate_id, fungicide = group_key

        valid_values = group.loc[
            group["_is_valid"],
            "_numeric_ec50",
        ]

        total_count = len(group)
        valid_count = len(valid_values)

        missing_count = int(
            group["_is_missing"].sum()
        )

        invalid_count = int(
            group["_is_invalid"].sum()
        )

        mean_value = (
            float(valid_values.mean())
            if valid_count > 0
            else None
        )

        sd_value = (
            float(valid_values.std(ddof=1))
            if valid_count >= 2
            else None
        )

        cv_percent = (
            float((sd_value / mean_value) * 100)
            if valid_count >= 2 and mean_value > 0
            else None
        )

        incomplete = (
            valid_count < expected_replicates
        )

        unexpected_count = max(
            total_count - expected_replicates,
            0,
        )

        unexpected_replicates = (
            unexpected_count > 0
        )

        high_variability = (
            max_cv_percent is not None
            and cv_percent is not None
            and cv_percent > max_cv_percent
        )

        requires_review = (
            incomplete
            or unexpected_replicates
            or missing_count > 0
            or invalid_count > 0
            or high_variability
        )

        results.append({
            isolate_col: isolate_id,
            fungicide_col: fungicide,
            "observation_count": total_count,
            "valid_replicate_count": valid_count,
            "missing_count": missing_count,
            "invalid_count": invalid_count,
            "expected_replicates": expected_replicates,
            "mean_ec50": mean_value,
            "sd_ec50": sd_value,
            "cv_percent": cv_percent,
            "incomplete_replicates": bool(incomplete),
            "unexpected_replicates": bool(
                unexpected_replicates
            ),
            "unexpected_replicate_count": int(
                unexpected_count
            ),
            "high_variability": bool(high_variability),
            "requires_review": bool(requires_review),
        })

    return pd.DataFrame(results)