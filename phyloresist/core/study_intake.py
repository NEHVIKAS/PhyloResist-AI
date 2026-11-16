"""Unified study intake workflow for PhyloResist AI.

Combines study configuration, data loading, and structural
validation into one reproducible intake process.

Structural validation does not establish biological validity.
"""

from pathlib import Path
from typing import Any

import pandas as pd

from phyloresist.core.data_loader import (
    load_tabular_data,
    inspect_dataset,
)

from phyloresist.core.validation import (
    validate_identifiers,
    validate_required_columns,
    validate_numeric_values,
)

from phyloresist.pathogens.config import StudyConfig


def validate_study_dataset(
    dataframe: pd.DataFrame,
    config: StudyConfig,
) -> dict[str, Any]:
    """Validate a dataset against its study configuration."""

    config_result = config.validate()

    if not config_result["is_valid"]:
        raise ValueError(
            f"Invalid study configuration: "
            f"{config_result['errors']}"
        )

    required_columns = [
        config.isolate_id_column,
        *config.phenotype_columns,
    ]

    column_result = validate_required_columns(
        dataframe,
        required_columns,
    )

    report = {
        "study_id": config.study_id,
        "pathogen": config.pathogen,
        "structural_validation": {
            "configuration": config_result,
            "required_columns": column_result,
        },
        "identifiers": None,
        "phenotypes": {},
        "dataset_summary": inspect_dataset(dataframe),
        "is_structurally_valid": False,
    }

    if not column_result["is_valid"]:
        return report

    identifier_result = validate_identifiers(
        dataframe[config.isolate_id_column].tolist()
    )

    report["identifiers"] = identifier_result

    phenotype_results = {}

    for column in config.phenotype_columns:
        phenotype_results[column] = validate_numeric_values(
            dataframe[column]
        )

    report["phenotypes"] = phenotype_results

    identifiers_valid = identifier_result["is_valid"]

    phenotypes_valid = all(
        result["is_valid"]
        for result in phenotype_results.values()
    )

    report["is_structurally_valid"] = (
        identifiers_valid and phenotypes_valid
    )

    return report


def run_study_intake(
    data_path: str | Path,
    config: StudyConfig,
    identifier_columns: list[str] | None = None,
) -> dict[str, Any]:
    """Load and structurally validate a study dataset."""

    dataframe = load_tabular_data(
        file_path=data_path,
        identifier_columns=identifier_columns,
    )

    report = validate_study_dataset(
        dataframe=dataframe,
        config=config,
    )

    report["source_file"] = str(data_path)

    return report