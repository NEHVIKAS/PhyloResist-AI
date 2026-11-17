import pandas as pd

from phyloresist.pathogens.config import StudyConfig
from phyloresist.validation.study_qc import run_study_qc


def make_config():
    return StudyConfig(
        study_id="FG_NEW_FEATURES",
        pathogen="Fusarium graminearum",
        isolate_id_column="isolate_id",
        phenotype_columns=["ec50"],
        fungicides=["propiconazole", "tebuconazole"],
        phenotype_unit="mM",
    )


def test_duplicate_isolate_ids_fail():
    data = pd.DataFrame({
        "isolate_id": ["FG001", "FG002", "FG001", "FG004"],
        "ec50": [0.4, 0.5, 0.6, 0.7],
    })

    report = run_study_qc(
        data,
        make_config(),
        data_level="isolate_level",
    )

    identifier_result = report["results"]["identifier_qc"]

    assert identifier_result["status"] == "FAIL"
    assert identifier_result["results"]["duplicate_identifiers"] == {
        "FG001": 2
    }
    assert report["overall_status"] == "FAIL"


def test_missing_isolate_ids_fail():
    data = pd.DataFrame({
        "isolate_id": ["FG001", None, "FG003", "FG004"],
        "ec50": [0.4, 0.5, 0.6, 0.7],
    })

    report = run_study_qc(
        data,
        make_config(),
        data_level="isolate_level",
    )

    identifier_result = report["results"]["identifier_qc"]

    assert identifier_result["status"] == "FAIL"
    assert identifier_result["results"]["missing_count"] == 1
    assert report["overall_status"] == "FAIL"


def test_fungicide_outliers_are_screened_independently():
    data = pd.DataFrame({
        "isolate_id": [
            "P001", "P002", "P003", "P004", "P005",
            "T001", "T002", "T003", "T004", "T005",
        ],
        "fungicide": (
            ["propiconazole"] * 5
            + ["tebuconazole"] * 5
        ),
        "replicate": [1] * 10,
        "ec50": [
            1.0, 1.1, 1.2, 1.3, 100.0,
            1000.0, 1100.0, 1200.0, 1300.0, 1400.0,
        ],
    })

    report = run_study_qc(
        data,
        make_config(),
        data_level="replicate_level",
        expected_replicates=1,
    )

    phenotype_result = report["results"]["phenotype_qc"]

    prop_result = phenotype_result["results"]["propiconazole"]
    tebu_result = phenotype_result["results"]["tebuconazole"]

    assert phenotype_result["status"] == "REVIEW"

    assert prop_result["requires_review"] is True
    assert prop_result["outlier_indices"] == [4]

    assert tebu_result["requires_review"] is False
    assert tebu_result["outlier_indices"] == []

    assert report["overall_status"] == "REVIEW"


def test_replicate_level_identifier_repetition_is_allowed():
    data = pd.DataFrame({
        "isolate_id": ["FG001"] * 3,
        "fungicide": ["propiconazole"] * 3,
        "replicate": [1, 2, 3],
        "ec50": [0.4, 0.45, 0.5],
    })

    report = run_study_qc(
        data,
        make_config(),
        data_level="replicate_level",
    )

    assert report["overall_status"] == "PASS"
    assert "identifier_qc" in report["modules_skipped"]