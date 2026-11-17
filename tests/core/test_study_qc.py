import json

import pandas as pd
import pytest

from phyloresist.pathogens.config import StudyConfig
from phyloresist.validation.study_qc import run_study_qc


@pytest.fixture
def study_config():
    return StudyConfig(
        study_id="FG_QC_TEST",
        pathogen="Fusarium graminearum",
        isolate_id_column="isolate_id",
        phenotype_columns=["ec50"],
        fungicides=["propiconazole"],
        phenotype_unit="mM",
    )


@pytest.fixture
def isolate_data():
    return pd.DataFrame({
        "isolate_id": [
            "FG001",
            "FG002",
            "FG003",
            "FG004",
            "FG005",
        ],
        "ec50": [
            0.45,
            0.78,
            1.20,
            0.56,
            1.45,
        ],
    })


@pytest.fixture
def replicate_data():
    rows = []

    values = {
        "FG001": [0.40, 0.45, 0.50],
        "FG002": [0.70, 0.75, 0.80],
    }

    for isolate_id, measurements in values.items():
        for index, value in enumerate(measurements, start=1):
            rows.append({
                "isolate_id": isolate_id,
                "fungicide": "propiconazole",
                "replicate": index,
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


def test_isolate_level_qc_passes(
    isolate_data,
    study_config,
):
    report = run_study_qc(
        isolate_data,
        study_config,
        data_level="isolate_level",
    )

    assert report["overall_status"] == "PASS"
    assert "phenotype_qc" in report["modules_executed"]
    assert "replicate_qc" in report["modules_skipped"]


def test_replicate_level_qc_passes(
    replicate_data,
    study_config,
):
    report = run_study_qc(
        replicate_data,
        study_config,
        data_level="replicate_level",
    )

    assert report["overall_status"] == "PASS"
    assert "phenotype_qc" in report["modules_executed"]
    assert "replicate_qc" in report["modules_executed"]


def test_dose_response_qc_passes(
    dose_response_data,
    study_config,
):
    report = run_study_qc(
        dose_response_data,
        study_config,
        data_level="dose_response",
    )

    assert report["overall_status"] == "PASS"
    assert (
        "experimental_design_qc"
        in report["modules_executed"]
    )


def test_unsupported_data_level(
    isolate_data,
    study_config,
):
    with pytest.raises(ValueError):
        run_study_qc(
            isolate_data,
            study_config,
            data_level="unknown",
        )


def test_empty_dataset_rejected(study_config):
    with pytest.raises(ValueError):
        run_study_qc(
            pd.DataFrame(),
            study_config,
            data_level="isolate_level",
        )


def test_non_dataframe_rejected(study_config):
    with pytest.raises(TypeError):
        run_study_qc(
            [],
            study_config,
            data_level="isolate_level",
        )


def test_invalid_configuration_rejected(isolate_data):
    config = StudyConfig(
        study_id="",
        pathogen="Fusarium graminearum",
        isolate_id_column="isolate_id",
        phenotype_columns=["ec50"],
        fungicides=["propiconazole"],
        phenotype_unit="mM",
    )

    with pytest.raises(ValueError):
        run_study_qc(
            isolate_data,
            config,
            data_level="isolate_level",
        )


def test_missing_phenotype_column(
    isolate_data,
    study_config,
):
    data = isolate_data.drop(columns=["ec50"])

    report = run_study_qc(
        data,
        study_config,
        data_level="isolate_level",
    )

    assert report["overall_status"] == "FAIL"
    assert report["errors"]


def test_replicate_missing_values_trigger_review(
    replicate_data,
    study_config,
):
    data = replicate_data.copy()
    data.loc[0, "ec50"] = None

    report = run_study_qc(
        data,
        study_config,
        data_level="replicate_level",
    )

    assert report["overall_status"] == "REVIEW"


def test_dose_response_missing_control_triggers_review(
    dose_response_data,
    study_config,
):
    data = dose_response_data[
        dose_response_data["concentration"] != 0
    ].copy()

    report = run_study_qc(
        data,
        study_config,
        data_level="dose_response",
    )

    assert report["overall_status"] == "REVIEW"


def test_missing_replicate_columns_fail(
    isolate_data,
    study_config,
):
    report = run_study_qc(
        isolate_data,
        study_config,
        data_level="replicate_level",
    )

    assert report["overall_status"] == "FAIL"
    assert report["errors"]


def test_report_is_json_serializable(
    isolate_data,
    study_config,
):
    report = run_study_qc(
        isolate_data,
        study_config,
        data_level="isolate_level",
    )

    serialized = json.dumps(report)

    assert isinstance(serialized, str)
    assert json.loads(serialized)["study_id"] == "FG_QC_TEST"


def test_original_dataframe_not_modified(
    isolate_data,
    study_config,
):
    original = isolate_data.copy(deep=True)

    run_study_qc(
        isolate_data,
        study_config,
        data_level="isolate_level",
    )

    pd.testing.assert_frame_equal(
        isolate_data,
        original,
    )


def test_report_contains_study_metadata(
    isolate_data,
    study_config,
):
    report = run_study_qc(
        isolate_data,
        study_config,
        data_level="isolate_level",
    )

    assert report["study_id"] == "FG_QC_TEST"
    assert report["pathogen"] == "Fusarium graminearum"
    assert report["data_level"] == "isolate_level"


def test_dose_response_modules_are_marked_skipped(
    dose_response_data,
    study_config,
):
    report = run_study_qc(
        dose_response_data,
        study_config,
        data_level="dose_response",
    )

    assert "phenotype_qc" in report["modules_skipped"]
    assert "replicate_qc" in report["modules_skipped"]


def test_isolate_level_modules_are_marked_skipped(
    isolate_data,
    study_config,
):
    report = run_study_qc(
        isolate_data,
        study_config,
        data_level="isolate_level",
    )

    assert "replicate_qc" in report["modules_skipped"]
    assert "experimental_design_qc" in report["modules_skipped"]