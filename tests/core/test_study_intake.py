import pandas as pd

from phyloresist.core.study_intake import (
    validate_study_dataset,
    run_study_intake,
)

from phyloresist.pathogens.config import StudyConfig


def create_test_config():
    return StudyConfig(
        study_id="FG_TEST_2026",
        pathogen="Fusarium graminearum",
        isolate_id_column="isolate_id",
        phenotype_columns=[
            "propiconazole_ec50",
            "tebuconazole_ec50",
        ],
        fungicides=[
            "propiconazole",
            "tebuconazole",
        ],
        phenotype_unit="mM",
    )


def test_valid_study_dataset():

    dataframe = pd.DataFrame({
        "isolate_id": ["FG001", "FG002", "FG003"],
        "propiconazole_ec50": [0.5, 1.2, 0.8],
        "tebuconazole_ec50": [0.3, 0.9, 1.1],
    })

    report = validate_study_dataset(
        dataframe,
        create_test_config(),
    )

    assert report["is_structurally_valid"] is True
    assert report["identifiers"]["unique_count"] == 3


def test_missing_required_column():

    dataframe = pd.DataFrame({
        "isolate_id": ["FG001", "FG002"],
        "propiconazole_ec50": [0.5, 1.2],
    })

    report = validate_study_dataset(
        dataframe,
        create_test_config(),
    )

    assert report["is_structurally_valid"] is False

    assert "tebuconazole_ec50" in (
        report["structural_validation"]
        ["required_columns"]["missing_columns"]
    )


def test_duplicate_identifiers():

    dataframe = pd.DataFrame({
        "isolate_id": ["FG001", "FG001"],
        "propiconazole_ec50": [0.5, 1.2],
        "tebuconazole_ec50": [0.3, 0.9],
    })

    report = validate_study_dataset(
        dataframe,
        create_test_config(),
    )

    assert report["is_structurally_valid"] is False

    assert report["identifiers"]["duplicate_identifiers"] == {
        "FG001": 2
    }


def test_invalid_phenotype_values():

    dataframe = pd.DataFrame({
        "isolate_id": ["FG001", "FG002"],
        "propiconazole_ec50": [0.5, "invalid"],
        "tebuconazole_ec50": [0.3, 0.9],
    })

    report = validate_study_dataset(
        dataframe,
        create_test_config(),
    )

    assert report["is_structurally_valid"] is False

    assert (
        report["phenotypes"]["propiconazole_ec50"]
        ["invalid_count"] == 1
    )


def test_intake_from_csv(tmp_path):

    dataframe = pd.DataFrame({
        "Isolate ID": ["FG001", "FG002"],
        "Propiconazole EC50": [0.5, 1.2],
        "Tebuconazole EC50": [0.3, 0.9],
    })

    file_path = tmp_path / "study.csv"

    dataframe.to_csv(file_path, index=False)

    config = create_test_config()

    report = run_study_intake(
        file_path,
        config,
        identifier_columns=["Isolate ID"],
    )

    assert report["is_structurally_valid"] is True
    assert report["source_file"] == str(file_path)