"""Shared tabular data-loading utilities for PhyloResist AI.

Supports CSV and Excel files, standardized column names,
identifier preservation, and basic dataset inspection.

The loader does not modify the original source files or
make assumptions about pathogen-specific data structures.
"""

from pathlib import Path
from typing import Any

import pandas as pd


SUPPORTED_EXTENSIONS = {".csv", ".xlsx", ".xls"}


def standardize_column_names(
    dataframe: pd.DataFrame,
) -> pd.DataFrame:
    """Standardize column names without changing their meaning.

    Operations:
    - Remove leading and trailing whitespace.
    - Convert names to lowercase.
    - Replace spaces and hyphens with underscores.

    Raises:
        ValueError: If standardization produces duplicate names.
    """

    standardized = dataframe.copy()

    new_columns = [
        str(column)
        .strip()
        .lower()
        .replace(" ", "_")
        .replace("-", "_")
        for column in standardized.columns
    ]

    if len(new_columns) != len(set(new_columns)):
        duplicates = sorted({
            column
            for column in new_columns
            if new_columns.count(column) > 1
        })

        raise ValueError(
            "Column standardization produced duplicate names: "
            f"{duplicates}"
        )

    standardized.columns = new_columns

    return standardized


def load_tabular_data(
    file_path: str | Path,
    sheet_name: str | int = 0,
    identifier_columns: list[str] | None = None,
    encoding: str = "utf-8",
) -> pd.DataFrame:
    """Load a CSV or Excel dataset into a pandas DataFrame.

    Args:
        file_path: Path to the source dataset.
        sheet_name: Excel sheet name or index. Ignored for CSV.
        identifier_columns: Original column names that must be
            imported as strings to preserve leading zeros.
        encoding: Text encoding used for CSV files.

    Returns:
        A DataFrame with standardized column names.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file extension is unsupported.
        ValueError: If the dataset has no columns.
    """

    path = Path(file_path)

    if not path.exists():
        raise FileNotFoundError(
            f"Dataset file does not exist: {path}"
        )

    if not path.is_file():
        raise ValueError(
            f"Expected a file, received a directory: {path}"
        )

    extension = path.suffix.lower()

    if extension not in SUPPORTED_EXTENSIONS:
        raise ValueError(
            f"Unsupported file extension: {extension}. "
            f"Supported formats: {sorted(SUPPORTED_EXTENSIONS)}"
        )

    dtype_mapping = {
        column: "string"
        for column in (identifier_columns or [])
    }

    if extension == ".csv":
        dataframe = pd.read_csv(
            path,
            encoding=encoding,
            dtype=dtype_mapping or None,
        )

    else:
        dataframe = pd.read_excel(
            path,
            sheet_name=sheet_name,
            dtype=dtype_mapping or None,
        )

    if not isinstance(dataframe, pd.DataFrame):
        raise ValueError(
            "Expected one worksheet. Multiple-sheet imports "
            "are not supported by this function."
        )

    if dataframe.empty and len(dataframe.columns) == 0:
        raise ValueError("The dataset contains no columns.")

    return standardize_column_names(dataframe)


def inspect_dataset(
    dataframe: pd.DataFrame,
) -> dict[str, Any]:
    """Generate a basic structural summary of a dataset.

    This function does not perform biological interpretation
    or modify the supplied DataFrame.
    """

    missing_values = dataframe.isna().sum()

    return {
        "row_count": len(dataframe),
        "column_count": len(dataframe.columns),
        "columns": dataframe.columns.tolist(),
        "missing_values": {
            column: int(count)
            for column, count in missing_values.items()
        },
        "total_missing_values": int(missing_values.sum()),
        "duplicate_rows": int(
            dataframe.duplicated().sum()
        ),
        "memory_usage_bytes": int(
            dataframe.memory_usage(deep=True).sum()
        ),
    }
