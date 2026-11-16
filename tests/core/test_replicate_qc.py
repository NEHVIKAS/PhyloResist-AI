import pandas as pd
import pytest

from phyloresist.validation.replicate_qc import summarize_replicates


@pytest.fixture
def valid_data():
    return pd.DataFrame({
        "isolate_id": [
            "FG001", "FG001", "FG001",
            "FG002", "FG002", "FG002",
        ],
        "fungicide": [
            "propiconazole", "propiconazole", "propiconazole",
            "propiconazole", "propiconazole", "propiconazole",
        ],
        "replicate": [1, 2, 3, 1, 2, 3],
        "ec50": [0.42, 0.45, 0.48, 0.75, 0.78, 0.81],
    })


def test_valid_replicates(valid_data):
    result = summarize_replicates(valid_data)

    assert len(result) == 2
    assert result["valid_replicate_count"].tolist() == [3, 3]


def test_mean_calculation(valid_data):
    result = summarize_replicates(valid_data)

    fg001 = result[result["isolate_id"] == "FG001"].iloc[0]

    assert fg001["mean_ec50"] == pytest.approx(0.45)


def test_sd_calculation(valid_data):
    result = summarize_replicates(valid_data)

    fg001 = result[result["isolate_id"] == "FG001"].iloc[0]

    assert fg001["sd_ec50"] == pytest.approx(0.03)


def test_cv_calculation(valid_data):
    result = summarize_replicates(valid_data)

    fg001 = result[result["isolate_id"] == "FG001"].iloc[0]

    expected_cv = (0.03 / 0.45) * 100

    assert fg001["cv_percent"] == pytest.approx(expected_cv)


def test_incomplete_replicates():
    data = pd.DataFrame({
        "isolate_id": ["FG001", "FG001"],
        "fungicide": ["propiconazole", "propiconazole"],
        "replicate": [1, 2],
        "ec50": [0.42, 0.45],
    })

    result = summarize_replicates(data)

    assert result.iloc[0]["incomplete_replicates"]
    assert result.iloc[0]["requires_review"]


def test_missing_measurement():
    data = pd.DataFrame({
        "isolate_id": ["FG001", "FG001", "FG001"],
        "fungicide": ["propiconazole"] * 3,
        "replicate": [1, 2, 3],
        "ec50": [0.42, None, 0.48],
    })

    result = summarize_replicates(data)

    assert result.iloc[0]["missing_count"] == 1
    assert result.iloc[0]["valid_replicate_count"] == 2
    assert result.iloc[0]["requires_review"]


def test_non_numeric_measurement():
    data = pd.DataFrame({
        "isolate_id": ["FG001", "FG001", "FG001"],
        "fungicide": ["propiconazole"] * 3,
        "replicate": [1, 2, 3],
        "ec50": [0.42, "invalid", 0.48],
    })

    result = summarize_replicates(data)

    assert result.iloc[0]["invalid_count"] == 1
    assert result.iloc[0]["requires_review"]


def test_zero_and_negative_measurements():
    data = pd.DataFrame({
        "isolate_id": ["FG001", "FG001", "FG001"],
        "fungicide": ["propiconazole"] * 3,
        "replicate": [1, 2, 3],
        "ec50": [0.42, 0, -0.5],
    })

    result = summarize_replicates(data)

    assert result.iloc[0]["invalid_count"] == 2
    assert result.iloc[0]["valid_replicate_count"] == 1


def test_high_variability_flag():
    data = pd.DataFrame({
        "isolate_id": ["FG001", "FG001", "FG001"],
        "fungicide": ["propiconazole"] * 3,
        "replicate": [1, 2, 3],
        "ec50": [0.2, 0.5, 1.0],
    })

    result = summarize_replicates(
        data,
        max_cv_percent=20,
    )

    assert result.iloc[0]["high_variability"]
    assert result.iloc[0]["requires_review"]


def test_no_variability_flag_without_threshold(valid_data):
    result = summarize_replicates(valid_data)

    assert not result["high_variability"].any()


def test_duplicate_replicate_identifiers():
    data = pd.DataFrame({
        "isolate_id": ["FG001", "FG001"],
        "fungicide": ["propiconazole", "propiconazole"],
        "replicate": [1, 1],
        "ec50": [0.42, 0.45],
    })

    with pytest.raises(ValueError, match="Duplicate"):
        summarize_replicates(data)


def test_missing_required_column(valid_data):
    data = valid_data.drop(columns=["ec50"])

    with pytest.raises(ValueError, match="Missing required columns"):
        summarize_replicates(data)


def test_invalid_expected_replicates(valid_data):
    with pytest.raises(ValueError):
        summarize_replicates(
            valid_data,
            expected_replicates=0,
        )


def test_invalid_cv_threshold(valid_data):
    with pytest.raises(ValueError):
        summarize_replicates(
            valid_data,
            max_cv_percent=-5,
        )


def test_original_data_is_not_modified(valid_data):
    original = valid_data.copy(deep=True)

    summarize_replicates(valid_data)

    pd.testing.assert_frame_equal(valid_data, original)


def test_multiple_fungicides():
    data = pd.DataFrame({
        "isolate_id": ["FG001"] * 6,
        "fungicide": [
            "propiconazole"] * 3 + ["tebuconazole"] * 3,
        "replicate": [1, 2, 3, 1, 2, 3],
        "ec50": [0.42, 0.45, 0.48, 0.30, 0.32, 0.34],
    })

    result = summarize_replicates(data)

    assert len(result) == 2
    assert set(result["fungicide"]) == {
        "propiconazole",
        "tebuconazole",
    }