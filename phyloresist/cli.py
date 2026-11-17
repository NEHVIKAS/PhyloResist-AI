"""Command-line interface for PhyloResist AI."""

import argparse
import json
import sys
from pathlib import Path

from phyloresist.core.data_loader import load_tabular_data
from phyloresist.core.study_intake import run_study_intake
from phyloresist.pathogens.config import StudyConfig
from phyloresist.pathogens.registry import StudyRegistry
from phyloresist.validation.study_qc import run_study_qc


DEFAULT_REGISTRY = Path("data/registry/studies.json")


def get_registry():
    """Initialize the central study registry."""
    return StudyRegistry(DEFAULT_REGISTRY)


def save_json_report(report, output_path):
    """Save a report as formatted JSON."""

    output_path = Path(output_path)

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with output_path.open("w", encoding="utf-8") as file:
        json.dump(
            report,
            file,
            indent=4,
            allow_nan=False,
        )

    print(f"Report saved to: {output_path}")


def command_list(args):
    """List all registered studies."""

    registry = get_registry()
    studies = registry.list_studies()

    if not studies:
        print("No studies are currently registered.")
        return 0

    print(f"\nRegistered studies: {len(studies)}\n")

    for study in studies:
        print(f"Study ID: {study['study_id']}")
        print(f"Pathogen: {study['pathogen']}")
        print(f"Fungicides: {', '.join(study['fungicides'])}")
        print("-" * 45)

    return 0


def command_show(args):
    """Display a study configuration."""

    registry = get_registry()

    try:
        config = registry.get(args.study_id)
    except KeyError as error:
        print(f"Error: {error}")
        return 1

    print(json.dumps(config.to_dict(), indent=4))

    return 0


def command_register(args):
    """Register a study from a JSON configuration."""

    try:
        config = StudyConfig.load_json(args.config)

        registry = get_registry()
        registry.register(config)

    except (ValueError, TypeError, FileNotFoundError) as error:
        print(f"Registration failed: {error}")
        return 1

    print(f"Successfully registered: {config.study_id}")

    return 0


def command_validate(args):
    """Validate a dataset using a registered study."""

    registry = get_registry()

    try:
        config = registry.get(args.study_id)

        report = run_study_intake(
            data_path=args.data,
            config=config,
            identifier_columns=[config.isolate_id_column],
        )

    except (KeyError, ValueError, TypeError, FileNotFoundError) as error:
        print(f"Validation failed: {error}")
        return 1

    if args.output:
        try:
            save_json_report(report, args.output)
        except (OSError, TypeError, ValueError) as error:
            print(f"Could not save validation report: {error}")
            return 1

    print("\nStudy intake report")
    print("=" * 45)
    print(f"Study: {report['study_id']}")
    print(f"Pathogen: {report['pathogen']}")
    print(f"Rows: {report['dataset_summary']['row_count']}")
    print(f"Columns: {report['dataset_summary']['column_count']}")

    print(
        "Structural validation:",
        "PASSED" if report["is_structurally_valid"] else "FAILED",
    )

    return 0 if report["is_structurally_valid"] else 1


def command_qc(args):
    """Run unified study quality control."""

    registry = get_registry()

    try:
        config = registry.get(args.study_id)

        data = load_tabular_data(
            args.data,
            identifier_columns=[config.isolate_id_column],
        )

        report = run_study_qc(
            data=data,
            config=config,
            data_level=args.level,
            expected_unit=args.expected_unit,
            isolate_col=args.isolate_col,
            fungicide_col=args.fungicide_col,
            replicate_col=args.replicate_col,
            value_col=args.value_col,
            concentration_col=args.concentration_col,
            response_col=args.response_col,
            expected_replicates=args.expected_replicates,
            min_treatment_levels=args.min_treatment_levels,
            max_replicate_cv_percent=args.max_replicate_cv,
            max_control_cv_percent=args.max_control_cv,
        )

    except (KeyError, ValueError, TypeError, FileNotFoundError) as error:
        print(f"QC failed: {error}")
        return 1

    if args.output:
        try:
            save_json_report(report, args.output)
        except (OSError, TypeError, ValueError) as error:
            print(f"Could not save QC report: {error}")
            return 1

    print("\nPhyloResist AI - Study QC Report")
    print("=" * 50)

    print(f"Study: {report['study_id']}")
    print(f"Pathogen: {report['pathogen']}")
    print(f"Data level: {report['data_level']}")
    print(f"Observations: {report['dataset_summary']['row_count']}")

    print("\nModules executed:")

    if report["modules_executed"]:
        for module in report["modules_executed"]:
            print(f"  [EXECUTED] {module}")
    else:
        print("  None")

    print("\nModules skipped:")

    for module, reason in report["modules_skipped"].items():
        print(f"  [SKIPPED] {module}: {reason}")

    print("\nWarnings:")

    if report["warnings"]:
        for warning in report["warnings"]:
            print(f"  - {warning}")
    else:
        print("  None")

    print("\nErrors:")

    if report["errors"]:
        for error in report["errors"]:
            print(f"  - {error}")
    else:
        print("  None")

    print("\nOverall QC status")
    print("=" * 50)
    print(f"  {report['overall_status']}")

    if report["overall_status"] == "PASS":
        return 0

    if report["overall_status"] == "REVIEW":
        return 2

    return 1


def build_parser():
    """Build the command-line argument parser."""

    parser = argparse.ArgumentParser(
        prog="phyloresist",
        description="Multi-pathogen fungicide resistance research platform.",
    )

    subparsers = parser.add_subparsers(
        dest="command",
        required=True,
    )

    # List studies
    list_parser = subparsers.add_parser(
        "list",
        help="List registered studies.",
    )
    list_parser.set_defaults(func=command_list)

    # Show study
    show_parser = subparsers.add_parser(
        "show",
        help="Display a study configuration.",
    )
    show_parser.add_argument("study_id")
    show_parser.set_defaults(func=command_show)

    # Register study
    register_parser = subparsers.add_parser(
        "register",
        help="Register a study from a JSON configuration.",
    )
    register_parser.add_argument("config")
    register_parser.set_defaults(func=command_register)

    # Validate dataset
    validate_parser = subparsers.add_parser(
        "validate",
        help="Validate a dataset using a registered study.",
    )

    validate_parser.add_argument("study_id")
    validate_parser.add_argument("data")

    validate_parser.add_argument(
        "--output",
        help="Optional JSON report output path.",
    )

    validate_parser.set_defaults(func=command_validate)

    # Unified quality control
    qc_parser = subparsers.add_parser(
        "qc",
        help="Run unified study quality control.",
    )

    qc_parser.add_argument("study_id")
    qc_parser.add_argument("data")

    qc_parser.add_argument(
        "--level",
        required=True,
        choices=[
            "isolate_level",
            "replicate_level",
            "dose_response",
        ],
        help="Declared data level.",
    )

    qc_parser.add_argument(
        "--output",
        help="Optional JSON QC report output path.",
    )

    qc_parser.add_argument(
        "--expected-unit",
        help="Expected EC50 measurement unit.",
    )

    qc_parser.add_argument(
        "--isolate-col",
        help="Override the isolate identifier column.",
    )

    qc_parser.add_argument(
        "--fungicide-col",
        default="fungicide",
        help="Fungicide column name.",
    )

    qc_parser.add_argument(
        "--replicate-col",
        default="replicate",
        help="Replicate identifier column.",
    )

    qc_parser.add_argument(
        "--value-col",
        help="EC50 measurement column for replicate-level data.",
    )

    qc_parser.add_argument(
        "--concentration-col",
        default="concentration",
        help="Concentration column for dose-response data.",
    )

    qc_parser.add_argument(
        "--response-col",
        default="growth",
        help="Growth response column for dose-response data.",
    )

    qc_parser.add_argument(
        "--expected-replicates",
        type=int,
        default=3,
        help="Expected number of experimental replicates.",
    )

    qc_parser.add_argument(
        "--min-treatment-levels",
        type=int,
        default=4,
        help="Minimum number of non-control concentrations.",
    )

    qc_parser.add_argument(
        "--max-replicate-cv",
        type=float,
        help="Optional maximum replicate EC50 CV percentage.",
    )

    qc_parser.add_argument(
        "--max-control-cv",
        type=float,
        help="Optional maximum control CV percentage.",
    )

    qc_parser.set_defaults(func=command_qc)

    return parser


def main():
    """Run the command-line interface."""

    parser = build_parser()
    args = parser.parse_args()

    try:
        return args.func(args)

    except Exception as error:
        print(
            f"Unexpected error: {error}",
            file=sys.stderr,
        )
        return 1


if __name__ == "__main__":
    sys.exit(main())