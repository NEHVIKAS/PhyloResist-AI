"""Shared validation utilities for PhyloResist AI.

These functions provide consistent validation of biological isolate
identifiers, tabular schemas, numerical measurements, and study metadata.
"""

from collections import Counter
from typing import Any

import pandas as pd


def normalize_identifier(value: Any) -> str:
    """Normalize an isolate identifier without removing leading zeros.

    Examples:
        125.0      -> "125"
        "  FG001 " -> "FG001"
        "00125"    -> "00125"
        None       -> ""
    """
    if value is None or pd.isna(value):
        return ""

    identifier = str(value).strip()

    if not identifier:
        return ""

    # Remove a decimal suffix only when it represents an integer.
    if identifier.endswith(".0"):
        numeric_part = identifier[:-2]
        if numeric_part.isdigit():
            identifier = numeric_part

    return identifier


def validate_identifiers(
    values: list[Any],
) -> dict[str, Any]:
    """Check identifiers for missing and duplicate values."""

    normalized = [normalize_identifier(value) for value in values]

    missing = [
        index for index, value in enumerate(normalized)
        if not value
    ]

    counts = Counter(
        value for value in normalized if value
    )

    duplicates = {
        identifier: count
        for identifier, count in counts.items()
        if count > 1
    }

    return {
        "total": len(normalized),
        "valid": len(normalized) - len(missing),
        "missing_count": len(missing),
        "missing_indices": missing,
        "duplicate_identifiers": duplicates,
        "unique_count": len(counts),
        "is_valid": not missing and not duplicates,
    }


def validate_required_columns(
    dataframe: pd.DataFrame,
    required_columns: list[str],
) -> dict[str, Any]:
    """Check whether a dataframe contains required columns."""

    missing = [
        column for column in required_columns
        if column not in dataframe.columns
    ]

    return {
        "required_columns": required_columns,
        "missing_columns": missing,
        "is_valid": not missing,
    }


def validate_numeric_values(
    values: pd.Series,
    allow_missing: bool = False,
) -> dict[str, Any]:
    """Check whether measurements are numeric and finite."""

    numeric = pd.to_numeric(values, errors="coerce")

    invalid = numeric.isna() & values.notna()
    missing = values.isna()
    non_finite = numeric.notna() & ~numeric.map(
        lambda value: abs(value) != float("inf")
    )

    invalid_count = int(invalid.sum())
    missing_count = int(missing.sum())
    non_finite_count = int(non_finite.sum())

    valid = (
        invalid_count == 0
        and non_finite_count == 0
        and (allow_missing or missing_count == 0)
    )

    return {
        "total": len(values),
        "invalid_count": invalid_count,
        "missing_count": missing_count,
        "non_finite_count": non_finite_count,
        "is_valid": valid,
    }


def validate_study_metadata(
    metadata: dict[str, Any],
) -> dict[str, Any]:
    """Perform basic structural validation of study metadata."""

    required = ["study_id", "pathogen", "fungicides"]

    missing = [
        field for field in required
        if field not in metadata or metadata[field] is None
    ]

    errors = []

    if "fungicides" in metadata:
        fungicides = metadata["fungicides"]

        if not isinstance(fungicides, list) or not fungicides:
            errors.append(
                "'fungicides' must be a non-empty list."
            )

    if "study_id" in metadata:
        if not isinstance(metadata["study_id"], str):
            errors.append("'study_id' must be a string.")

    if "pathogen" in metadata:
        if not isinstance(metadata["pathogen"], str):
            errors.append("'pathogen' must be a string.")

    return {
        "missing_fields": missing,
        "errors": errors,
        "is_valid": not missing and not errors,
    }