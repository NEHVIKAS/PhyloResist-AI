import pandas as pd
import pytest

from phyloresist.core.data_loader import (
    standardize_column_names,
    load_tabular_data,
    inspect_dataset,
)


# --------------------------------------------------
# 1. Column-name standardization
# --------------------------------------------------

def test_standardize_column_names():
    dataframe = pd.DataFrame({
        "Isolate ID": ["FG001"],
        "EC50-Value": [0.5],
        " Fungicide ": ["tebuconazole"],
    })

    result = standardize_column_names(dataframe)

    assert result.columns.tolist() == [
        "isolate_id",
        "ec50_value",
        "fungicide",
    ]


def test_standardization_detects_duplicate_columns():
    dataframe = pd.DataFrame(
        [["FG001", "FG002"]],
        columns=["Isolate ID", "isolate-id"],
    )

    with pytest.raises(ValueError):
        standardize_column_names(dataframe)


# --------------------------------------------------
# 2. CSV loading
# --------------------------------------------------

def test_load_csv(tmp_path):
    file_path = tmp_path / "test_data.csv"

    dataframe = pd.DataFrame({
        "Isolate ID": ["FG001", "FG002"],
        "EC50": [0.5, 1.2],
    })

    dataframe.to_csv(file_path, index=False)

    result = load_tabular_data(file_path)

    assert result.shape == (2, 2)
    assert "isolate_id" in result.columns
    assert "ec50" in result.columns


def test_csv_preserves_leading_zeros(tmp_path):
    file_path = tmp_path / "identifiers.csv"

    file_path.write_text(
        "Isolate ID,EC50\n00125,0.5\n00126,1.2\n",
        encoding="utf-8",
    )

    result = load_tabular_data(
        file_path,
        identifier_columns=["Isolate ID"],
    )

    assert result["isolate_id"].tolist() == [
        "00125",
        "00126",
    ]


# --------------------------------------------------
# 3. Excel loading
# --------------------------------------------------

def test_load_excel(tmp_path):
    pytest.importorskip("openpyxl")

    file_path = tmp_path / "test_data.xlsx"

    dataframe = pd.DataFrame({
        "Isolate ID": ["FG001", "FG002"],
        "EC50": [0.5, 1.2],
    })

    dataframe.to_excel(file_path, index=False)

    result = load_tabular_data(file_path)

    assert result.shape == (2, 2)
    assert "isolate_id" in result.columns


# --------------------------------------------------
# 4. File validation
# --------------------------------------------------

def test_missing_file():
    with pytest.raises(FileNotFoundError):
        load_tabular_data("nonexistent_dataset.csv")


def test_unsupported_file(tmp_path):
    file_path = tmp_path / "dataset.txt"
    file_path.write_text("test", encoding="utf-8")

    with pytest.raises(ValueError):
        load_tabular_data(file_path)


# --------------------------------------------------
# 5. Dataset inspection
# --------------------------------------------------

def test_inspect_dataset():
    dataframe = pd.DataFrame({
        "isolate_id": ["FG001", "FG002", "FG003"],
        "ec50": [0.5, None, 1.2],
        "fungicide": [
            "tebuconazole",
            "tebuconazole",
            "propiconazole",
        ],
    })

    result = inspect_dataset(dataframe)

    assert result["row_count"] == 3
    assert result["column_count"] == 3
    assert result["total_missing_values"] == 1
    assert result["missing_values"]["ec50"] == 1
    assert result["duplicate_rows"] == 0


def test_inspect_duplicate_rows():
    dataframe = pd.DataFrame({
        "isolate_id": ["FG001", "FG001"],
        "ec50": [0.5, 0.5],
    })

    result = inspect_dataset(dataframe)

    assert result["duplicate_rows"] == 1