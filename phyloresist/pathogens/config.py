"""Pathogen-specific study configuration for PhyloResist AI.

Stores study metadata and pathogen-specific experimental settings.
Configurations are independent of the analysis methods.
"""

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class StudyConfig:
    """Configuration describing one pathogen research study."""

    study_id: str
    pathogen: str

    isolate_id_column: str

    phenotype_columns: list[str]
    fungicides: list[str]

    phenotype_unit: str | None = None
    genotype_format: str | None = None
    reference_genome: str | None = None

    population_column: str | None = None
    metadata_columns: list[str] = field(default_factory=list)

    schema_version: str = "1.0"

    def validate(self) -> dict[str, Any]:
        """Validate the structure of the study configuration."""

        errors = []

        required_strings = {
            "study_id": self.study_id,
            "pathogen": self.pathogen,
            "isolate_id_column": self.isolate_id_column,
        }

        for field_name, value in required_strings.items():
            if not isinstance(value, str) or not value.strip():
                errors.append(
                    f"'{field_name}' must be a non-empty string."
                )

        if not isinstance(self.phenotype_columns, list):
            errors.append(
                "'phenotype_columns' must be a list."
            )
        elif not self.phenotype_columns:
            errors.append(
                "'phenotype_columns' cannot be empty."
            )
        elif not all(
            isinstance(column, str) and column.strip()
            for column in self.phenotype_columns
        ):
            errors.append(
                "Every phenotype column must be a non-empty string."
            )

        if not isinstance(self.fungicides, list):
            errors.append(
                "'fungicides' must be a list."
            )
        elif not all(
            isinstance(item, str) and item.strip()
            for item in self.fungicides
        ):
            errors.append(
                "Every fungicide must be a non-empty string."
            )

        if not isinstance(self.metadata_columns, list):
            errors.append(
                "'metadata_columns' must be a list."
            )

        return {
            "study_id": self.study_id,
            "pathogen": self.pathogen,
            "errors": errors,
            "is_valid": len(errors) == 0,
        }

    def to_dict(self) -> dict[str, Any]:
        """Convert configuration into a dictionary."""

        return asdict(self)

    def save_json(self, file_path: str | Path) -> Path:
        """Save the configuration as a JSON file."""

        validation = self.validate()

        if not validation["is_valid"]:
            raise ValueError(
                f"Invalid study configuration: {validation['errors']}"
            )

        path = Path(file_path)
        path.parent.mkdir(parents=True, exist_ok=True)

        with path.open("w", encoding="utf-8") as file:
            json.dump(
                self.to_dict(),
                file,
                indent=4,
                ensure_ascii=False,
            )

        return path

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "StudyConfig":
        """Create a configuration from a dictionary."""

        if not isinstance(data, dict):
            raise TypeError("Configuration data must be a dictionary.")

        allowed_fields = set(cls.__dataclass_fields__)
        unknown_fields = set(data) - allowed_fields

        if unknown_fields:
            raise ValueError(
                f"Unknown configuration fields: {sorted(unknown_fields)}"
            )

        return cls(**data)

    @classmethod
    def load_json(cls, file_path: str | Path) -> "StudyConfig":
        """Load a study configuration from JSON."""

        path = Path(file_path)

        if not path.exists():
            raise FileNotFoundError(
                f"Configuration file not found: {path}"
            )

        with path.open("r", encoding="utf-8") as file:
            data = json.load(file)

        config = cls.from_dict(data)

        validation = config.validate()

        if not validation["is_valid"]:
            raise ValueError(
                f"Invalid study configuration: {validation['errors']}"
            )

        return config