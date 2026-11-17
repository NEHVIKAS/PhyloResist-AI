import json
from types import SimpleNamespace

import pytest

from phyloresist import cli


class MockConfig:
    def __init__(self):
        self.study_id = "FG_TEST"
        self.pathogen = "Fusarium graminearum"
        self.isolate_id_column = "isolate_id"
        self.fungicides = ["propiconazole"]

    def to_dict(self):
        return {
            "study_id": self.study_id,
            "pathogen": self.pathogen,
            "isolate_id_column": self.isolate_id_column,
            "fungicides": self.fungicides,
        }


class MockRegistry:
    def __init__(self, studies=None):
        self.studies = studies or {}

    def list_studies(self):
        # Match the actual StudyRegistry interface:
        # list_studies() returns dictionaries.
        return [
            study.to_dict()
            for study in self.studies.values()
        ]

    def get(self, study_id):
        if study_id not in self.studies:
            raise KeyError(f"Study not found: {study_id}")

        return self.studies[study_id]

    def register(self, config):
        if config.study_id in self.studies:
            raise ValueError("Study already registered")

        self.studies[config.study_id] = config


@pytest.fixture
def config():
    return MockConfig()


@pytest.fixture
def registry(config):
    return MockRegistry({"FG_TEST": config})


@pytest.fixture
def sample_qc_report():
    return {
        "study_id": "FG_TEST",
        "pathogen": "Fusarium graminearum",
        "data_level": "isolate_level",
        "dataset_summary": {
            "row_count": 5,
        },
        "modules_executed": ["phenotype_qc"],
        "modules_skipped": {},
        "warnings": [],
        "errors": [],
        "overall_status": "PASS",
    }


# --------------------------------------------------
# Parser tests
# --------------------------------------------------

def test_parser_registers_all_commands():
    parser = cli.build_parser()

    expected_commands = {
        "list",
        "show",
        "register",
        "validate",
        "qc",
    }

    command_arguments = {
        "list": ["list"],
        "show": ["show", "FG_TEST"],
        "register": ["register", "config.json"],
        "validate": ["validate", "FG_TEST", "data.csv"],
        "qc": [
            "qc",
            "FG_TEST",
            "data.csv",
            "--level",
            "isolate_level",
        ],
    }

    for command in expected_commands:
        args = parser.parse_args(command_arguments[command])
        assert args.command == command


def test_qc_parser_requires_level():
    parser = cli.build_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(["qc", "FG_TEST", "data.csv"])


def test_qc_parser_accepts_options():
    parser = cli.build_parser()

    args = parser.parse_args([
        "qc",
        "FG_TEST",
        "data.csv",
        "--level",
        "replicate_level",
        "--expected-unit",
        "mM",
        "--expected-replicates",
        "3",
        "--max-replicate-cv",
        "25",
    ])

    assert args.level == "replicate_level"
    assert args.expected_unit == "mM"
    assert args.expected_replicates == 3
    assert args.max_replicate_cv == 25.0


# --------------------------------------------------
# List command tests
# --------------------------------------------------

def test_list_empty_registry(monkeypatch, capsys):
    monkeypatch.setattr(
        cli,
        "get_registry",
        lambda: MockRegistry(),
    )

    result = cli.command_list(SimpleNamespace())

    output = capsys.readouterr().out

    assert result == 0
    assert "No studies are currently registered." in output


def test_list_registered_studies(monkeypatch, registry, capsys):
    monkeypatch.setattr(
        cli,
        "get_registry",
        lambda: registry,
    )

    result = cli.command_list(SimpleNamespace())

    output = capsys.readouterr().out

    assert result == 0
    assert "FG_TEST" in output
    assert "Fusarium graminearum" in output
    assert "propiconazole" in output


# --------------------------------------------------
# Show command tests
# --------------------------------------------------

def test_show_existing_study(monkeypatch, registry, capsys):
    monkeypatch.setattr(
        cli,
        "get_registry",
        lambda: registry,
    )

    result = cli.command_show(
        SimpleNamespace(study_id="FG_TEST")
    )

    output = capsys.readouterr().out

    assert result == 0
    assert '"study_id": "FG_TEST"' in output


def test_show_unknown_study(monkeypatch, registry, capsys):
    monkeypatch.setattr(
        cli,
        "get_registry",
        lambda: registry,
    )

    result = cli.command_show(
        SimpleNamespace(study_id="UNKNOWN")
    )

    output = capsys.readouterr().out

    assert result == 1
    assert "Error:" in output


# --------------------------------------------------
# Register command tests
# --------------------------------------------------

def test_register_study(monkeypatch, config, capsys):
    test_registry = MockRegistry()

    monkeypatch.setattr(
        cli,
        "get_registry",
        lambda: test_registry,
    )

    monkeypatch.setattr(
        cli.StudyConfig,
        "load_json",
        lambda path: config,
    )

    result = cli.command_register(
        SimpleNamespace(config="config.json")
    )

    output = capsys.readouterr().out

    assert result == 0
    assert "Successfully registered: FG_TEST" in output
    assert "FG_TEST" in test_registry.studies


def test_register_invalid_configuration(monkeypatch, capsys):
    monkeypatch.setattr(
        cli.StudyConfig,
        "load_json",
        lambda path: (_ for _ in ()).throw(
            ValueError("Invalid configuration")
        ),
    )

    result = cli.command_register(
        SimpleNamespace(config="invalid.json")
    )

    output = capsys.readouterr().out

    assert result == 1
    assert "Registration failed:" in output


# --------------------------------------------------
# Validate command tests
# --------------------------------------------------

def test_validate_success(monkeypatch, registry, capsys):
    monkeypatch.setattr(
        cli,
        "get_registry",
        lambda: registry,
    )

    monkeypatch.setattr(
        cli,
        "run_study_intake",
        lambda **kwargs: {
            "study_id": "FG_TEST",
            "pathogen": "Fusarium graminearum",
            "dataset_summary": {
                "row_count": 5,
                "column_count": 3,
            },
            "is_structurally_valid": True,
        },
    )

    result = cli.command_validate(
        SimpleNamespace(
            study_id="FG_TEST",
            data="data.csv",
            output=None,
        )
    )

    output = capsys.readouterr().out

    assert result == 0
    assert "Structural validation: PASSED" in output


def test_validate_invalid_dataset(monkeypatch, registry, capsys):
    monkeypatch.setattr(
        cli,
        "get_registry",
        lambda: registry,
    )

    monkeypatch.setattr(
        cli,
        "run_study_intake",
        lambda **kwargs: {
            "study_id": "FG_TEST",
            "pathogen": "Fusarium graminearum",
            "dataset_summary": {
                "row_count": 5,
                "column_count": 2,
            },
            "is_structurally_valid": False,
        },
    )

    result = cli.command_validate(
        SimpleNamespace(
            study_id="FG_TEST",
            data="data.csv",
            output=None,
        )
    )

    output = capsys.readouterr().out

    assert result == 1
    assert "Structural validation: FAILED" in output


def test_validate_missing_file(monkeypatch, registry, capsys):
    monkeypatch.setattr(
        cli,
        "get_registry",
        lambda: registry,
    )

    monkeypatch.setattr(
        cli,
        "run_study_intake",
        lambda **kwargs: (_ for _ in ()).throw(
            FileNotFoundError("data.csv not found")
        ),
    )

    result = cli.command_validate(
        SimpleNamespace(
            study_id="FG_TEST",
            data="data.csv",
            output=None,
        )
    )

    output = capsys.readouterr().out

    assert result == 1
    assert "Validation failed:" in output


# --------------------------------------------------
# QC command tests
# --------------------------------------------------

@pytest.mark.parametrize(
    ("status", "expected_exit_code"),
    [
        ("PASS", 0),
        ("REVIEW", 2),
        ("FAIL", 1),
    ],
)
def test_qc_status_exit_codes(
    monkeypatch,
    registry,
    sample_qc_report,
    status,
    expected_exit_code,
):
    monkeypatch.setattr(
        cli,
        "get_registry",
        lambda: registry,
    )

    monkeypatch.setattr(
        cli,
        "load_tabular_data",
        lambda *args, **kwargs: object(),
    )

    report = sample_qc_report.copy()
    report["overall_status"] = status

    monkeypatch.setattr(
        cli,
        "run_study_qc",
        lambda **kwargs: report,
    )

    args = SimpleNamespace(
        study_id="FG_TEST",
        data="data.csv",
        level="isolate_level",
        expected_unit=None,
        isolate_col=None,
        fungicide_col="fungicide",
        replicate_col="replicate",
        value_col=None,
        concentration_col="concentration",
        response_col="growth",
        expected_replicates=3,
        min_treatment_levels=4,
        max_replicate_cv=None,
        max_control_cv=None,
        output=None,
    )

    result = cli.command_qc(args)

    assert result == expected_exit_code


def test_qc_saves_json_report(
    monkeypatch,
    registry,
    sample_qc_report,
    tmp_path,
):
    monkeypatch.setattr(
        cli,
        "get_registry",
        lambda: registry,
    )

    monkeypatch.setattr(
        cli,
        "load_tabular_data",
        lambda *args, **kwargs: object(),
    )

    monkeypatch.setattr(
        cli,
        "run_study_qc",
        lambda **kwargs: sample_qc_report,
    )

    output_path = tmp_path / "reports" / "qc_report.json"

    args = SimpleNamespace(
        study_id="FG_TEST",
        data="data.csv",
        level="isolate_level",
        expected_unit=None,
        isolate_col=None,
        fungicide_col="fungicide",
        replicate_col="replicate",
        value_col=None,
        concentration_col="concentration",
        response_col="growth",
        expected_replicates=3,
        min_treatment_levels=4,
        max_replicate_cv=None,
        max_control_cv=None,
        output=str(output_path),
    )

    result = cli.command_qc(args)

    assert result == 0
    assert output_path.exists()

    with output_path.open("r", encoding="utf-8") as file:
        saved_report = json.load(file)

    assert saved_report["study_id"] == "FG_TEST"
    assert saved_report["overall_status"] == "PASS"


def test_qc_missing_input_file(monkeypatch, registry, capsys):
    monkeypatch.setattr(
        cli,
        "get_registry",
        lambda: registry,
    )

    monkeypatch.setattr(
        cli,
        "load_tabular_data",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            FileNotFoundError("data.csv not found")
        ),
    )

    args = SimpleNamespace(
        study_id="FG_TEST",
        data="data.csv",
        level="isolate_level",
        expected_unit=None,
        isolate_col=None,
        fungicide_col="fungicide",
        replicate_col="replicate",
        value_col=None,
        concentration_col="concentration",
        response_col="growth",
        expected_replicates=3,
        min_treatment_levels=4,
        max_replicate_cv=None,
        max_control_cv=None,
        output=None,
    )

    result = cli.command_qc(args)

    output = capsys.readouterr().out

    assert result == 1
    assert "QC failed:" in output


# --------------------------------------------------
# JSON report tests
# --------------------------------------------------

def test_save_json_report_creates_parent_directories(tmp_path):
    output_path = tmp_path / "nested" / "report.json"

    report = {
        "study_id": "FG_TEST",
        "overall_status": "PASS",
    }

    cli.save_json_report(report, output_path)

    assert output_path.exists()

    with output_path.open("r", encoding="utf-8") as file:
        saved_report = json.load(file)

    assert saved_report == report