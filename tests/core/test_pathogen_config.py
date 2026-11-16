import pytest

from phyloresist.pathogens.config import StudyConfig


@pytest.fixture
def fusarium_config():
    return StudyConfig(
        study_id="FG_2025",
        pathogen="Fusarium graminearum",
        isolate_id_column="isolate_id",
        phenotype_columns=["propiconazole_ec50", "tebuconazole_ec50"],
        fungicides=["propiconazole", "tebuconazole"],
        phenotype_unit="mM",
        genotype_format="VCF",
        reference_genome=None,
        population_column="population",
        metadata_columns=["country", "location"],
    )


def test_valid_configuration(fusarium_config):
    result = fusarium_config.validate()

    assert result["is_valid"] is True
    assert result["errors"] == []


def test_configuration_to_dictionary(fusarium_config):
    result = fusarium_config.to_dict()

    assert result["study_id"] == "FG_2025"
    assert result["pathogen"] == "Fusarium graminearum"
    assert len(result["phenotype_columns"]) == 2


def test_configuration_from_dictionary(fusarium_config):
    data = fusarium_config.to_dict()

    restored = StudyConfig.from_dict(data)

    assert restored.study_id == fusarium_config.study_id
    assert restored.pathogen == fusarium_config.pathogen


def test_configuration_json_roundtrip(
    fusarium_config,
    tmp_path,
):
    file_path = tmp_path / "study_config.json"

    fusarium_config.save_json(file_path)

    restored = StudyConfig.load_json(file_path)

    assert restored.study_id == fusarium_config.study_id
    assert restored.fungicides == fusarium_config.fungicides


def test_missing_configuration_file():
    with pytest.raises(FileNotFoundError):
        StudyConfig.load_json("missing_config.json")


def test_unknown_configuration_field():
    with pytest.raises(ValueError):
        StudyConfig.from_dict({
            "study_id": "FG_2025",
            "unknown_field": "unexpected",
        })


def test_invalid_configuration():
    config = StudyConfig(
        study_id="",
        pathogen="Fusarium graminearum",
        isolate_id_column="isolate_id",
        phenotype_columns=[],
        fungicides=["propiconazole"],
    )

    result = config.validate()

    assert result["is_valid"] is False
    assert len(result["errors"]) >= 2