import numpy as np
import pandas as pd
import pytest

from phyloresist.pathogens.config import StudyConfig
from phyloresist.validation.study_qc import run_study_qc


@pytest.fixture
def study_config():
    return StudyConfig(
        study_id="FG_DESIGN_EDGE",
        pathogen="Fusarium graminearum",
        isolate_id_column="isolate_id",
        phenotype_columns=["ec50"],
        fungicides=["propiconazole"],
        phenotype_unit="mM",
    )


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


def test_missing_replicate_triggers_review(
    dose_response_data,
    study_config,
):
    data = dose_response_data.drop(
        dose_response_data[
            (dose_response_data["concentration"] == 1)
            & (dose_response_data["replicate"] == 3)
        ].index
    )

    report = run_study_qc(
        data,
        study_config,
        data_level="dose_response",
        expected_replicates=3,
    )

    result = report["results"]["experimental_design_qc"]["results"][0]

    assert report["overall_status"] == "REVIEW"
    assert result["treatment_replicates_complete"] is False


def test_negative_concentration_fails(
    dose_response_data,
    study_config,
):
    data = dose_response_data.copy()
    data.loc[3, "concentration"] = -1

    report = run_study_qc(
        data,
        study_config,
        data_level="dose_response",
    )

    result = report["results"]["experimental_design_qc"]["results"][0]

    assert result["invalid_concentration_count"] == 1
    assert report["overall_status"] == "FAIL"


def test_infinite_growth_measurement_fails(
    dose_response_data,
    study_config,
):
    data = dose_response_data.copy()
    data.loc[3, "growth"] = np.inf

    report = run_study_qc(
        data,
        study_config,
        data_level="dose_response",
    )

    result = report["results"]["experimental_design_qc"]["results"][0]

    assert result["invalid_response_count"] == 1
    assert report["overall_status"] == "FAIL"


def test_high_control_variability_triggers_review(
    dose_response_data,
    study_config,
):
    data = dose_response_data.copy()

    control_indices = data[
        data["concentration"] == 0
    ].index

    data.loc[control_indices, "growth"] = [5.0, 10.0, 15.0]

    report = run_study_qc(
        data,
        study_config,
        data_level="dose_response",
        max_control_cv_percent=20,
    )

    result = report["results"]["experimental_design_qc"]["results"][0]

    assert result["high_control_variability"] is True
    assert report["overall_status"] == "REVIEW"