"""Quality-control utilities for fungicide dose-response experimental design."""

import numpy as np
import pandas as pd


def audit_dose_response_design(
    data: pd.DataFrame,
    isolate_col: str = "isolate_id",
    fungicide_col: str = "fungicide",
    concentration_col: str = "concentration",
    replicate_col: str = "replicate",
    response_col: str = "growth",
    control_concentration: float = 0.0,
    expected_replicates: int = 3,
    min_treatment_levels: int = 4,
    max_control_cv_percent: float | None = None,
) -> pd.DataFrame:

    if not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a pandas DataFrame.")

    if data.empty:
        raise ValueError("data must contain at least one observation.")

    columns = [
        isolate_col,
        fungicide_col,
        concentration_col,
        replicate_col,
        response_col,
    ]

    if len(set(columns)) != len(columns):
        raise ValueError("Column names must be distinct.")

    missing_columns = [
        column for column in columns
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

    if (
        isinstance(min_treatment_levels, bool)
        or not isinstance(min_treatment_levels, int)
        or min_treatment_levels < 1
    ):
        raise ValueError(
            "min_treatment_levels must be a positive integer."
        )

    if (
        isinstance(control_concentration, bool)
        or not isinstance(control_concentration, (int, float))
        or not np.isfinite(control_concentration)
        or control_concentration < 0
    ):
        raise ValueError(
            "control_concentration must be a non-negative finite number."
        )

    if max_control_cv_percent is not None:
        if (
            isinstance(max_control_cv_percent, bool)
            or not isinstance(max_control_cv_percent, (int, float))
            or not np.isfinite(max_control_cv_percent)
            or max_control_cv_percent < 0
        ):
            raise ValueError(
                "max_control_cv_percent must be None or a "
                "non-negative finite number."
            )

    working = data[columns].copy()

    if working[[isolate_col, fungicide_col]].isna().any().any():
        raise ValueError(
            "Isolate and fungicide identifiers cannot be missing."
        )

    if working[replicate_col].isna().any():
        raise ValueError(
            "Replicate identifiers cannot be missing."
        )

    concentration_raw = working[concentration_col]
    response_raw = working[response_col]

    concentration = pd.to_numeric(
        concentration_raw,
        errors="coerce",
    )

    response = pd.to_numeric(
        response_raw,
        errors="coerce",
    )

    concentration_missing = concentration_raw.isna()

    concentration_valid = (
        concentration.notna()
        & np.isfinite(concentration)
        & (concentration >= 0)
    )

    concentration_invalid = (
        ~concentration_missing & ~concentration_valid
    )

    response_missing = response_raw.isna()

    response_valid = (
        response.notna()
        & np.isfinite(response)
    )

    response_invalid = (
        ~response_missing & ~response_valid
    )

    working["_concentration"] = concentration
    working["_response"] = response

    working["_concentration_valid"] = concentration_valid
    working["_concentration_missing"] = concentration_missing
    working["_concentration_invalid"] = concentration_invalid

    working["_response_valid"] = response_valid
    working["_response_missing"] = response_missing
    working["_response_invalid"] = response_invalid

    results = []

    group_columns = [isolate_col, fungicide_col]

    for group_key, group in working.groupby(
        group_columns,
        sort=True,
    ):
        isolate_id, fungicide = group_key

        valid_concentration_mask = group["_concentration_valid"]

        control_mask = (
            valid_concentration_mask
            & (group["_concentration"] == control_concentration)
        )

        treatment_mask = (
            valid_concentration_mask
            & (group["_concentration"] != control_concentration)
        )

        control = group.loc[control_mask]
        treatments = group.loc[treatment_mask]

        treatment_levels = sorted(
            treatments["_concentration"].unique().tolist()
        )

        control_count = len(control)
        treatment_count = len(treatments)

        control_valid = control.loc[
            control["_response_valid"],
            "_response",
        ]

        treatment_valid = treatments.loc[
            treatments["_response_valid"],
            "_response",
        ]

        control_present = control_count > 0

        insufficient_treatment_levels = (
            len(treatment_levels) < min_treatment_levels
        )

        # Count valid response measurements per concentration.
        control_replicates_complete = (
            len(control_valid) >= expected_replicates
        )

        treatment_valid_counts = (
            treatments.loc[
                treatments["_response_valid"]
            ]
            .groupby("_concentration")
            .size()
        )

        treatment_replicates_complete = True

        if not treatment_valid_counts.empty:
            treatment_replicates_complete = bool(
                (
                    treatment_valid_counts >= expected_replicates
                ).all()
            )

        # A treatment level with zero valid responses is incomplete.
        if len(treatment_levels) > 0:
            treatment_replicates_complete = (
                treatment_replicates_complete
                and len(treatment_valid_counts) == len(treatment_levels)
            )

        control_mean = (
            float(control_valid.mean())
            if not control_valid.empty else None
        )

        control_sd = (
            float(control_valid.std(ddof=1))
            if len(control_valid) >= 2 else None
        )

        control_cv = (
            float(control_sd / abs(control_mean) * 100)
            if control_sd is not None
            and control_mean is not None
            and control_mean != 0
            else None
        )

        high_control_variability = (
            max_control_cv_percent is not None
            and control_cv is not None
            and control_cv > max_control_cv_percent
        )

        missing_concentration_count = int(
            group["_concentration_missing"].sum()
        )

        invalid_concentration_count = int(
            group["_concentration_invalid"].sum()
        )

        missing_response_count = int(
            group["_response_missing"].sum()
        )

        invalid_response_count = int(
            group["_response_invalid"].sum()
        )

        # Count additional duplicate records, not every member
        # of a duplicate group.
        duplicate_mask = group.duplicated(
            subset=[
                "_concentration",
                replicate_col,
            ],
            keep="first",
        )

        duplicate_count = int(duplicate_mask.sum())

        requires_review = (
            not control_present
            or not control_replicates_complete
            or insufficient_treatment_levels
            or not treatment_replicates_complete
            or missing_concentration_count > 0
            or invalid_concentration_count > 0
            or missing_response_count > 0
            or invalid_response_count > 0
            or duplicate_count > 0
            or high_control_variability
        )

        results.append({
            isolate_col: isolate_id,
            fungicide_col: fungicide,
            "control_present": bool(control_present),
            "control_observation_count": control_count,
            "control_valid_count": len(control_valid),
            "control_mean": control_mean,
            "control_sd": control_sd,
            "control_cv_percent": control_cv,
            "high_control_variability": bool(
                high_control_variability
            ),
            "treatment_level_count": len(treatment_levels),
            "treatment_concentrations": treatment_levels,
            "treatment_observation_count": treatment_count,
            "valid_treatment_response_count": len(treatment_valid),
            "expected_replicates": expected_replicates,
            "control_replicates_complete": bool(
                control_replicates_complete
            ),
            "treatment_replicates_complete": bool(
                treatment_replicates_complete
            ),
            "insufficient_treatment_levels": bool(
                insufficient_treatment_levels
            ),
            "missing_concentration_count": missing_concentration_count,
            "invalid_concentration_count": invalid_concentration_count,
            "missing_response_count": missing_response_count,
            "invalid_response_count": invalid_response_count,
            "duplicate_observation_count": duplicate_count,
            "requires_review": bool(requires_review),
        })

    return pd.DataFrame(results)