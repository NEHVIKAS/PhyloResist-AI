import pandas as pd
import pytest

from phyloresist.core.validation import (
    normalize_identifier,
    validate_identifiers,
    validate_required_columns,
    validate_numeric_values,
    validate_study_metadata,
)


# --------------------------------------------------
# 1. Identifier normalization
# --------------------------------------------------

@pytest.mark.parametrize(
    "value, expected",
    [
        (125.0, "125"),
        (" FG001 ", "FG001"),
        ("00125", "00125"),
        ("00125.0", "00125"),
        (None, ""),
        ("", ""),
    ],
)
def test_normalize_identifier(value, expected):
    assert normalize_identifier(value) == expected


# --------------------------------------------------
# 2. Identifier validation
# --------------------------------------------------

def test_validate_identifiers_with_duplicates():
    result = validate_identifiers(["FG001", "FG002", "FG001"])

    assert result["total"] == 3
    assert result["duplicate_identifiers"] == {"FG001": 2}
    assert result["is_valid"] is False


def test_validate_identifiers_with_missing_values():
    result = validate_identifiers(["FG001", None, "FG003"])

    assert result["missing_count"] == 1
    assert result["is_valid"] is False


def test_validate_identifiers_without_errors():
    result = validate_identifiers(["FG001", "FG002", "FG003"])

    assert result["unique_count"] == 3
    assert result["is_valid"] is True


# --------------------------------------------------
# 3. Required column validation
# --------------------------------------------------

def test_required_columns_present():
    dataframe = pd.DataFrame({
        "isolate_id": ["FG001", "FG002"],
        "ec50": [0.5, 1.2],
    })

    result = validate_required_columns(
        dataframe,
        ["isolate_id", "ec50"],
    )

    assert result["is_valid"] is True
    assert result["missing_columns"] == []


def test_required_columns_missing():
    dataframe = pd.DataFrame({
        "isolate_id": ["FG001", "FG002"],
    })

    result = validate_required_columns(
        dataframe,
        ["isolate_id", "ec50"],
    )

    assert result["is_valid"] is False
    assert result["missing_columns"] == ["ec50"]


# --------------------------------------------------
# 4. Numeric validation
# --------------------------------------------------

def test_numeric_values_valid():
    values = pd.Series([0.1, 1.0, 10.0, 100.0])

    result = validate_numeric_values(values)

    assert result["is_valid"] is True
    assert result["invalid_count"] == 0


def test_numeric_values_invalid():
    values = pd.Series([0.1, "invalid", 10.0])

    result = validate_numeric_values(values)

    assert result["is_valid"] is False
    assert result["invalid_count"] == 1


def test_numeric_values_with_missing():
    values = pd.Series([0.1, None, 10.0])

    result = validate_numeric_values(values)

    assert result["is_valid"] is False
    assert result["missing_count"] == 1


def test_numeric_values_allow_missing():
    values = pd.Series([0.1, None, 10.0])

    result = validate_numeric_values(
        values,
        allow_missing=True,
    )

    assert result["is_valid"] is True


def test_numeric_values_infinite():
    values = pd.Series([0.1, float("inf"), 10.0])

    result = validate_numeric_values(values)

    assert result["is_valid"] is False
    assert result["non_finite_count"] == 1


# --------------------------------------------------
# 5. Study metadata validation
# --------------------------------------------------

def test_valid_study_metadata():
    metadata = {
        "study_id": "FG_2025",
        "pathogen": "Fusarium graminearum",
        "fungicides": ["propiconazole", "tebuconazole"],
    }

    result = validate_study_metadata(metadata)

    assert result["is_valid"] is True


def test_invalid_study_metadata():
    metadata = {
        "study_id": "FG_2025",
        "pathogen": "Fusarium graminearum",
        "fungicides": [],
    }

    result = validate_study_metadata(metadata)

    assert result["is_valid"] is False
    assert len(result["errors"]) > 0