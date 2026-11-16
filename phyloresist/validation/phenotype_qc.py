"""Quality-control utilities for fungicide EC50 phenotype data."""

import numpy as np
import pandas as pd


def validate_ec50_values(
    values: pd.Series,
    unit: str,
    expected_unit: str | None = None,
    iqr_multiplier: float = 1.5,
) -> dict:
    """
    Validate EC50 measurements and flag potential outliers.

    Missing and outlier values are flagged for review.
    Values are never automatically removed or modified.
    """

    if not isinstance(values, pd.Series):
        raise TypeError("values must be a pandas Series.")

    if not isinstance(unit, str) or not unit.strip():
        raise ValueError("unit must be a non-empty string.")

    if expected_unit is not None and (
        not isinstance(expected_unit, str) or not expected_unit.strip()
    ):
        raise ValueError("expected_unit must be None or a non-empty string.")

    if (
        isinstance(iqr_multiplier, bool)
        or not isinstance(iqr_multiplier, (int, float))
        or not np.isfinite(iqr_multiplier)
        or iqr_multiplier <= 0
    ):
        raise ValueError("iqr_multiplier must be a positive finite number.")

    numeric = pd.to_numeric(values, errors="coerce")

    missing_mask = values.isna()
    invalid_mask = numeric.isna() & ~missing_mask
    non_finite_mask = numeric.notna() & ~np.isfinite(numeric)
    non_positive_mask = (
        numeric.notna()
        & np.isfinite(numeric)
        & (numeric <= 0)
    )

    valid_mask = (
        numeric.notna()
        & np.isfinite(numeric)
        & (numeric > 0)
    )

    valid_values = numeric[valid_mask]

    outlier_indices = []
    lower_fence = None
    upper_fence = None

    if len(valid_values) >= 4:
        log_values = np.log10(valid_values)

        q1 = log_values.quantile(0.25)
        q3 = log_values.quantile(0.75)
        iqr = q3 - q1

        lower_log_fence = q1 - iqr_multiplier * iqr
        upper_log_fence = q3 + iqr_multiplier * iqr

        lower_fence = float(10 ** lower_log_fence)
        upper_fence = float(10 ** upper_log_fence)

        outlier_mask = (
            (log_values < lower_log_fence)
            | (log_values > upper_log_fence)
        )

        outlier_indices = log_values.index[outlier_mask].tolist()

    unit_matches = (
        None
        if expected_unit is None
        else unit.strip() == expected_unit.strip()
    )

    total = len(values)

    return {
        "measurement": "EC50",
        "reported_unit": unit,
        "expected_unit": expected_unit,
        "unit_matches": unit_matches,
        "total_count": total,
        "missing_count": int(missing_mask.sum()),
        "invalid_count": int(invalid_mask.sum()),
        "non_finite_count": int(non_finite_mask.sum()),
        "non_positive_count": int(non_positive_mask.sum()),
        "valid_count": int(valid_mask.sum()),
        "outlier_count": len(outlier_indices),
        "outlier_indices": outlier_indices,
        "outlier_method": "log10_IQR",
        "iqr_multiplier": iqr_multiplier,
        "outlier_lower_fence": lower_fence,
        "outlier_upper_fence": upper_fence,
        "statistics": {
            "minimum": (
                float(valid_values.min())
                if not valid_values.empty else None
            ),
            "median": (
                float(valid_values.median())
                if not valid_values.empty else None
            ),
            "mean": (
                float(valid_values.mean())
                if not valid_values.empty else None
            ),
            "maximum": (
                float(valid_values.max())
                if not valid_values.empty else None
            ),
        },
        "is_structurally_valid": bool(
            invalid_mask.sum() == 0
            and non_finite_mask.sum() == 0
            and non_positive_mask.sum() == 0
            and (
                unit_matches is None
                or unit_matches
            )
        ),
        "requires_review": bool(
            len(outlier_indices) > 0
            or missing_mask.sum() > 0
        ),
    }