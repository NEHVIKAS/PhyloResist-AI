"""Central study registry for PhyloResist AI.

Maintains a persistent catalogue of pathogen-specific study
configurations.

The registry stores study metadata, not raw research datasets.
"""

import json
from pathlib import Path
from typing import Any

from phyloresist.pathogens.config import StudyConfig


class StudyRegistry:
    """Register, retrieve, and list research studies."""

    REGISTRY_VERSION = "1.0"

    def __init__(self, registry_path: str | Path):
        self.registry_path = Path(registry_path)
        self._studies: dict[str, dict[str, Any]] = {}

        if self.registry_path.exists():
            self._load()

    def _load(self) -> None:
        """Load the existing registry from JSON."""

        with self.registry_path.open("r", encoding="utf-8") as file:
            data = json.load(file)

        if not isinstance(data, dict):
            raise ValueError("Registry must contain a JSON object.")

        if data.get("registry_version") != self.REGISTRY_VERSION:
            raise ValueError("Unsupported registry version.")

        studies = data.get("studies")

        if not isinstance(studies, list):
            raise ValueError("'studies' must be a list.")

        loaded = {}

        for entry in studies:
            if not isinstance(entry, dict):
                raise ValueError("Each registry entry must be an object.")

            config = StudyConfig.from_dict(entry)
            validation = config.validate()

            if not validation["is_valid"]:
                raise ValueError(
                    f"Invalid registered study: {validation['errors']}"
                )

            if config.study_id in loaded:
                raise ValueError(
                    f"Duplicate study ID in registry: {config.study_id}"
                )

            loaded[config.study_id] = config.to_dict()

        self._studies = loaded

    def save(self) -> Path:
        """Persist the registry to JSON."""

        self.registry_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        data = {
            "registry_version": self.REGISTRY_VERSION,
            "study_count": len(self._studies),
            "studies": list(self._studies.values()),
        }

        temporary_path = self.registry_path.with_name(
            self.registry_path.name + ".tmp"
        )

        try:
            with temporary_path.open(
                "w",
                encoding="utf-8",
            ) as file:
                json.dump(
                    data,
                    file,
                    indent=4,
                    ensure_ascii=False,
                )

            temporary_path.replace(self.registry_path)

        finally:
            if temporary_path.exists():
                temporary_path.unlink()

        return self.registry_path

    def register(self, config: StudyConfig) -> None:
        """Register a new study configuration."""

        if not isinstance(config, StudyConfig):
            raise TypeError("Expected a StudyConfig object.")

        validation = config.validate()

        if not validation["is_valid"]:
            raise ValueError(
                f"Cannot register invalid study: {validation['errors']}"
            )

        study_id = config.study_id

        if study_id in self._studies:
            raise ValueError(
                f"Study already registered: {study_id}"
            )

        self._studies[study_id] = config.to_dict()

        try:
            self.save()
        except Exception:
            self._studies.pop(study_id, None)
            raise

    def get(self, study_id: str) -> StudyConfig:
        """Retrieve a registered study configuration."""

        if study_id not in self._studies:
            raise KeyError(f"Study not found: {study_id}")

        return StudyConfig.from_dict(
            self._studies[study_id]
        )

    def list_studies(self) -> list[dict[str, Any]]:
        """Return a summary of all registered studies."""

        summaries = []

        for study_id, data in sorted(self._studies.items()):
            summaries.append({
                "study_id": study_id,
                "pathogen": data["pathogen"],
                "phenotype_columns": data["phenotype_columns"],
                "fungicides": data["fungicides"],
            })

        return summaries

    def __len__(self) -> int:
        """Return the number of registered studies."""

        return len(self._studies)