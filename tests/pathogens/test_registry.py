import json
import pytest

from phyloresist.pathogens.config import StudyConfig
from phyloresist.pathogens.registry import StudyRegistry


@pytest.fixture
def fusarium_config():
    return StudyConfig(
        study_id="FG_2025",
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
        genotype_format="VCF",
    )


@pytest.fixture
def zymoseptoria_config():
    return StudyConfig(
        study_id="ZT_2026",
        pathogen="Zymoseptoria tritici",
        isolate_id_column="isolate_id",
        phenotype_columns=[
            "epoxiconazole_ec50",
            "prothioconazole_ec50",
        ],
        fungicides=[
            "epoxiconazole",
            "prothioconazole",
        ],
        phenotype_unit="mg/L",
        genotype_format=None,
    )


def test_register_study(tmp_path, fusarium_config):
    registry = StudyRegistry(tmp_path / "registry.json")

    registry.register(fusarium_config)

    assert len(registry) == 1


def test_retrieve_study(tmp_path, fusarium_config):
    registry = StudyRegistry(tmp_path / "registry.json")

    registry.register(fusarium_config)

    result = registry.get("FG_2025")

    assert result.pathogen == "Fusarium graminearum"
    assert result.study_id == "FG_2025"


def test_list_multiple_studies(
    tmp_path,
    fusarium_config,
    zymoseptoria_config,
):
    registry = StudyRegistry(tmp_path / "registry.json")

    registry.register(fusarium_config)
    registry.register(zymoseptoria_config)

    studies = registry.list_studies()

    assert len(studies) == 2

    study_ids = {
        study["study_id"]
        for study in studies
    }

    assert study_ids == {"FG_2025", "ZT_2026"}


def test_duplicate_registration(
    tmp_path,
    fusarium_config,
):
    registry = StudyRegistry(tmp_path / "registry.json")

    registry.register(fusarium_config)

    with pytest.raises(ValueError):
        registry.register(fusarium_config)


def test_missing_study(tmp_path):
    registry = StudyRegistry(tmp_path / "registry.json")

    with pytest.raises(KeyError):
        registry.get("UNKNOWN_STUDY")


def test_invalid_configuration_registration(tmp_path):
    registry = StudyRegistry(tmp_path / "registry.json")

    invalid_config = StudyConfig(
        study_id="",
        pathogen="Fusarium graminearum",
        isolate_id_column="isolate_id",
        phenotype_columns=[],
        fungicides=[],
    )

    with pytest.raises(ValueError):
        registry.register(invalid_config)

    assert len(registry) == 0


def test_registry_persistence(
    tmp_path,
    fusarium_config,
):
    registry_path = tmp_path / "registry.json"

    registry = StudyRegistry(registry_path)
    registry.register(fusarium_config)

    restored_registry = StudyRegistry(registry_path)

    assert len(restored_registry) == 1

    restored = restored_registry.get("FG_2025")

    assert restored.pathogen == "Fusarium graminearum"


def test_registry_creates_parent_directory(
    tmp_path,
    fusarium_config,
):
    registry_path = (
        tmp_path / "new_folder" / "registry.json"
    )

    registry = StudyRegistry(registry_path)
    registry.register(fusarium_config)

    assert registry_path.exists()


def test_corrupted_registry(tmp_path):
    registry_path = tmp_path / "registry.json"

    registry_path.write_text(
        '{"registry_version": "1.0", "studies": "invalid"}',
        encoding="utf-8",
    )

    with pytest.raises(ValueError):
        StudyRegistry(registry_path)


def test_registry_file_contains_correct_data(
    tmp_path,
    fusarium_config,
):
    registry_path = tmp_path / "registry.json"

    registry = StudyRegistry(registry_path)
    registry.register(fusarium_config)

    with registry_path.open(
        "r",
        encoding="utf-8",
    ) as file:
        data = json.load(file)

    assert data["registry_version"] == "1.0"
    assert data["study_count"] == 1
    assert data["studies"][0]["study_id"] == "FG_2025"