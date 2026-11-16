import numpy as np
import pandas as pd
import pytest

from phyloresist.validation.phenotype_qc import (
    validate_ec50_values,
)


def test_valid_ec50_values():

    values = pd.Series([0.1, 0.5, 1.0, 2.0])

    result = validate_ec50_values(
        values,
        unit="mg/L",
        expected_unit="mg/L",
    )

    assert result["is_structurally_valid"] is True
    assert result["valid_count"] == 4
    assert result["invalid_count"] == 0


def test_missing_ec50_values():

    values = pd.Series([0.1, None, 1.0])

    result = validate_ec50_values(
        values,
        unit="mg/L",
    )

    assert result["missing_count"] == 1
    assert result["valid_count"] == 2
    assert result["requires_review"] is True


def test_non_numeric_ec50_values():

    values = pd.Series([0.1, "invalid", 1.0])

    result = validate_ec50_values(
        values,
        unit="mg/L",
    )

    assert result["invalid_count"] == 1
    assert result["is_structurally_valid"] is False


def test_zero_ec50_values():

    values = pd.Series([0.1, 0, 1.0])

    result = validate_ec50_values(
        values,
        unit="mg/L",
    )

    assert result["non_positive_count"] == 1
    assert result["is_structurally_valid"] is False


def test_negative_ec50_values():

    values = pd.Series([0.1, -0.5, 1.0])

    result = validate_ec50_values(
        values,
        unit="mg/L",
    )

    assert result["non_positive_count"] == 1
    assert result["is_structurally_valid"] is False


def test_infinite_ec50_values():

    values = pd.Series([0.1, np.inf, 1.0])

    result = validate_ec50_values(
        values,
        unit="mg/L",
    )

    assert result["non_finite_count"] == 1
    assert result["is_structurally_valid"] is False


def test_unit_mismatch():

    values = pd.Series([0.1, 0.5, 1.0])

    result = validate_ec50_values(
        values,
        unit="mg/L",
        expected_unit="mM",
    )

    assert result["unit_matches"] is False
    assert result["is_structurally_valid"] is False


def test_outlier_detection():

    values = pd.Series([
        0.1, 0.2, 0.3, 0.4, 0.5,
        0.6, 0.7, 0.8, 0.9, 100.0,
    ])

    result = validate_ec50_values(
        values,
        unit="mg/L",
    )

    assert result["outlier_count"] >= 1
    assert result["requires_review"] is True

    # Outliers are flagged, not deleted.
    assert result["total_count"] == 10
    assert result["valid_count"] == 10


def test_empty_measurements():

    values = pd.Series([], dtype=float)

    result = validate_ec50_values(
        values,
        unit="mg/L",
    )

    assert result["total_count"] == 0
    assert result["valid_count"] == 0
    assert result["statistics"]["median"] is None


def test_invalid_unit():

    values = pd.Series([0.1, 0.5])

    with pytest.raises(ValueError):
        validate_ec50_values(
            values,
            unit="",
        )


def test_invalid_iqr_multiplier():

    values = pd.Series([0.1, 0.5, 1.0])

    with pytest.raises(ValueError):
        validate_ec50_values(
            values,
            unit="mg/L",
            iqr_multiplier=-1,
        )


def test_descriptive_statistics():

    values = pd.Series([1.0, 2.0, 3.0, 4.0])

    result = validate_ec50_values(
        values,
        unit="mg/L",
    )

    assert result["statistics"]["minimum"] == 1.0
    assert result["statistics"]["median"] == 2.5
    assert result["statistics"]["maximum"] == 4.0