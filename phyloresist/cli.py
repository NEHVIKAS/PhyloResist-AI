"""Command-line interface for PhyloResist AI."""

import argparse
import json
import sys
from pathlib import Path

from phyloresist.pathogens.config import StudyConfig
from phyloresist.pathogens.registry import StudyRegistry
from phyloresist.core.study_intake import run_study_intake


DEFAULT_REGISTRY = Path("data/registry/studies.json")


def get_registry():
    """Initialize the central study registry."""
    return StudyRegistry(DEFAULT_REGISTRY)


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
        output_path = Path(args.output)
        output_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        with output_path.open("w", encoding="utf-8") as file:
            json.dump(report, file, indent=4)

        print(f"Report saved to: {output_path}")

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

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    try:
        return args.func(args)
    except Exception as error:
        print(f"Unexpected error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())