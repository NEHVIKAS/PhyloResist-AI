import pandas as pd
import pytest

from phyloresist.validation.experimental_design_qc import (
    audit_dose_response_design,
)


@pytest.fixture
def valid_experiment():
    rows = []

    concentrations = [0, 0.1, 1, 10, 100]

    for concentration in concentrations:
        for replicate in range(1, 4):
            growth = (
                10.0
                if concentration == 0
                else 10.0 / (1 + concentration / 10)
            )

            rows.append(
                {
                    "isolate_id": "FG001",
                    "fungicide": "propiconazole",
                    "concentration": concentration,
                    "replicate": replicate,
                    "growth": growth + replicate * 0.1,
                }
            )

    return pd.DataFrame(rows)


def test_valid_experimental_design(valid_experiment):
    result = audit_dose_response_design(valid_experiment)

    assert len(result) == 1
    assert result.iloc[0]["control_present"]
    assert not result.iloc[0]["requires_review"]


def test_control_presence(valid_experiment):
    result = audit_dose_response_design(valid_experiment)

    assert result.iloc[0]["control_present"]


def test_missing_control(valid_experiment):
    data = valid_experiment[
        valid_experiment["concentration"] != 0
    ].copy()

    result = audit_dose_response_design(data)

    assert not result.iloc[0]["control_present"]
    assert result.iloc[0]["requires_review"]


def test_insufficient_treatment_levels(valid_experiment):
    data = valid_experiment[
        valid_experiment["concentration"].isin([0, 0.1, 1, 10])
    ].copy()

    result = audit_dose_response_design(
        data,
        min_treatment_levels=4,
    )

    assert result.iloc[0]["treatment_level_count"] == 3
    assert result.iloc[0]["requires_review"]


def test_treatment_concentrations(valid_experiment):
    result = audit_dose_response_design(valid_experiment)

    assert result.iloc[0]["treatment_level_count"] == 4


def test_control_statistics(valid_experiment):
    result = audit_dose_response_design(valid_experiment)

    assert result.iloc[0]["control_mean"] > 0
    assert result.iloc[0]["control_sd"] >= 0
    assert result.iloc[0]["control_cv_percent"] >= 0


def test_incomplete_control_replicates(valid_experiment):
    data = valid_experiment.drop(index=0)

    result = audit_dose_response_design(
        data,
        expected_replicates=3,
    )

    assert not result.iloc[0]["control_replicates_complete"]
    assert result.iloc[0]["requires_review"]


def test_incomplete_treatment_replicates(valid_experiment):
    data = valid_experiment.drop(index=3)

    result = audit_dose_response_design(
        data,
        expected_replicates=3,
    )

    assert not result.iloc[0]["treatment_replicates_complete"]
    assert result.iloc[0]["requires_review"]


def test_invalid_response(valid_experiment):
    data = valid_experiment.copy()
    data["growth"] = data["growth"].astype(object)
    data.loc[0, "growth"] = "invalid"

    result = audit_dose_response_design(data)

    assert result.iloc[0]["invalid_response_count"] == 1
    assert result.iloc[0]["missing_response_count"] == 0
    assert result.iloc[0]["requires_review"]


def test_missing_response(valid_experiment):
    data = valid_experiment.copy()
    data.loc[0, "growth"] = None

    result = audit_dose_response_design(data)

    assert result.iloc[0]["missing_response_count"] == 1
    assert result.iloc[0]["invalid_response_count"] == 0
    assert result.iloc[0]["requires_review"]


def test_high_control_variability(valid_experiment):
    data = valid_experiment.copy()

    data.loc[
        data["concentration"] == 0,
        "growth",
    ] = [1, 10, 30]

    result = audit_dose_response_design(
        data,
        max_control_cv_percent=20,
    )

    assert result.iloc[0]["control_cv_percent"] > 20
    assert result.iloc[0]["requires_review"]


def test_duplicate_observations(valid_experiment):
    duplicate = valid_experiment.iloc[[0]].copy()
    data = pd.concat(
        [valid_experiment, duplicate],
        ignore_index=True,
    )

    result = audit_dose_response_design(data)

    assert result.iloc[0]["duplicate_observation_count"] == 1
    assert result.iloc[0]["requires_review"]


def test_multiple_fungicides(valid_experiment):
    second = valid_experiment.copy()
    second["fungicide"] = "tebuconazole"

    data = pd.concat(
        [valid_experiment, second],
        ignore_index=True,
    )

    result = audit_dose_response_design(data)

    assert len(result) == 2


def test_missing_required_column(valid_experiment):
    data = valid_experiment.drop(columns=["growth"])

    with pytest.raises(ValueError):
        audit_dose_response_design(data)


def test_invalid_expected_replicates(valid_experiment):
    with pytest.raises(ValueError):
        audit_dose_response_design(
            valid_experiment,
            expected_replicates=0,
        )


def test_invalid_minimum_treatment_levels(valid_experiment):
    with pytest.raises(ValueError):
        audit_dose_response_design(
            valid_experiment,
            min_treatment_levels=0,
        )


def test_original_data_not_modified(valid_experiment):
    original = valid_experiment.copy(deep=True)

    audit_dose_response_design(valid_experiment)

    pd.testing.assert_frame_equal(
        valid_experiment,
        original,
    )

def test_missing_concentration(valid_experiment):
    data = valid_experiment.copy()
    data.loc[0, "concentration"] = None

    result = audit_dose_response_design(data)

    assert result.iloc[0]["missing_concentration_count"] == 1
    assert result.iloc[0]["invalid_concentration_count"] == 0
    assert result.iloc[0]["requires_review"]


def test_invalid_concentration(valid_experiment):
    data = valid_experiment.copy()
    data["concentration"] = data["concentration"].astype(object)
    data.loc[0, "concentration"] = "invalid"

    result = audit_dose_response_design(data)

    assert result.iloc[0]["invalid_concentration_count"] == 1
    assert result.iloc[0]["missing_concentration_count"] == 0
    assert result.iloc[0]["requires_review"]


def test_negative_concentration(valid_experiment):
    data = valid_experiment.copy()
    data.loc[0, "concentration"] = -1

    result = audit_dose_response_design(data)

    assert result.iloc[0]["invalid_concentration_count"] == 1
    assert result.iloc[0]["requires_review"]