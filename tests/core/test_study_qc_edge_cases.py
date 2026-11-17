import numpy as np
import pandas as pd
import pytest

from phyloresist.pathogens.config import StudyConfig
from phyloresist.validation.study_qc import run_study_qc


@pytest.fixture
def study_config():
    return StudyConfig(
        study_id="FG_EDGE_TEST",
        pathogen="Fusarium graminearum",
        isolate_id_column="isolate_id",
        phenotype_columns=["ec50"],
        fungicides=["propiconazole"],
        phenotype_unit="mM",
    )


@pytest.fixture
def isolate_data():
    return pd.DataFrame({
        "isolate_id": ["FG001", "FG002", "FG003", "FG004", "FG005"],
        "ec50": [0.45, 0.78, 1.20, 0.56, 1.45],
    })


@pytest.fixture
def replicate_data():
    rows = []

    values = {
        "FG001": [0.40, 0.45, 0.50],
        "FG002": [0.70, 0.75, 0.80],
    }

    for isolate_id, measurements in values.items():
        for replicate, value in enumerate(measurements, start=1):
            rows.append({
                "isolate_id": isolate_id,
                "fungicide": "propiconazole",
                "replicate": replicate,
                "ec50": value,
            })

    return pd.DataFrame(rows)


@pytest.fixture
def dose_response_data():
    rows = []

    concentrations = [0, 0.1, 1, 10, 100]

    for concentration in concentrations:
        for replicate in range(1, 4):
            growth = (
                10.0
                if concentration == 0
                else 10.0 / (1 + concentration / 10)
            )

            rows.append({
                "isolate_id": "FG001",
                "fungicide": "propiconazole",
                "concentration": concentration,
                "replicate": replicate,
                "growth": growth + replicate * 0.1,
            })

    return pd.DataFrame(rows)


@pytest.mark.parametrize(
    "invalid_value",
    [0, -0.5, np.inf, "not-a-number"],
)
def test_invalid_ec50_values_fail(
    isolate_data,
    study_config,
    invalid_value,
):
    data = isolate_data.copy()
    data["ec50"] = data["ec50"].astype(object)
    data.loc[0, "ec50"] = invalid_value

    report = run_study_qc(
        data,
        study_config,
        data_level="isolate_level",
    )

    assert report["overall_status"] == "FAIL"

    result = report["results"]["phenotype_qc"]["results"]["ec50"]

    assert result["is_structurally_valid"] is False


def test_unit_mismatch_fails(isolate_data, study_config):
    report = run_study_qc(
        isolate_data,
        study_config,
        data_level="isolate_level",
        expected_unit="mg/L",
    )

    assert report["overall_status"] == "FAIL"

    result = report["results"]["phenotype_qc"]["results"]["ec50"]

    assert result["is_structurally_valid"] is False


def test_duplicate_replicate_records_fail(
    replicate_data,
    study_config,
):
    duplicate = replicate_data.iloc[[0]].copy()

    data = pd.concat(
        [replicate_data, duplicate],
        ignore_index=True,
    )

    report = run_study_qc(
        data,
        study_config,
        data_level="replicate_level",
    )

    assert report["overall_status"] == "FAIL"
    assert "replicate_qc" in report["results"]


def test_multiple_phenotype_columns_are_checked_independently():
    config = StudyConfig(
        study_id="FG_MULTI_TEST",
        pathogen="Fusarium graminearum",
        isolate_id_column="isolate_id",
        phenotype_columns=[
            "propiconazole_ec50",
            "tebuconazole_ec50",
        ],
        fungicides=["propiconazole", "tebuconazole"],
        phenotype_unit="mM",
    )

    data = pd.DataFrame({
        "isolate_id": ["FG001", "FG002", "FG003", "FG004"],
        "propiconazole_ec50": [0.4, 0.5, 0.6, 0.7],
        "tebuconazole_ec50": [0.3, 0, 0.5, 0.6],
    })

    report = run_study_qc(
        data,
        config,
        data_level="isolate_level",
    )

    results = report["results"]["phenotype_qc"]["results"]

    assert results["propiconazole_ec50"]["is_structurally_valid"]
    assert not results["tebuconazole_ec50"]["is_structurally_valid"]
    assert report["overall_status"] == "FAIL"


def test_missing_treatment_level_triggers_review(
    dose_response_data,
    study_config,
):
    data = dose_response_data[
        dose_response_data["concentration"] != 100
    ].copy()

    report = run_study_qc(
        data,
        study_config,
        data_level="dose_response",
        min_treatment_levels=4,
    )

    assert report["overall_status"] == "REVIEW"