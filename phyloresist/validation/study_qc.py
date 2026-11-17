
"""Unified quality-control orchestration for PhyloResist AI.

Coordinates phenotype, replicate, and experimental-design QC
according to the explicitly declared study data level.

This module does not fit dose-response models or establish
biological fungicide resistance.
"""

from typing import Any

import pandas as pd

from phyloresist.core.validation import validate_identifiers
from phyloresist.pathogens.config import StudyConfig

from phyloresist.validation.phenotype_qc import (
    validate_ec50_values,
)

from phyloresist.validation.replicate_qc import (
    summarize_replicates,
)

from phyloresist.validation.experimental_design_qc import (
    audit_dose_response_design,
)


SUPPORTED_DATA_LEVELS = {
    "isolate_level",
    "replicate_level",
    "dose_response",
}


def _dataframe_records(dataframe: pd.DataFrame) -> list[dict]:
    """Convert a DataFrame into JSON-compatible records."""

    return dataframe.astype(object).where(
        pd.notna(dataframe),
        None,
    ).to_dict(orient="records")


def _determine_status(
    module_results: dict[str, Any],
    errors: list[str],
) -> str:
    """Determine the overall study QC status."""

    if errors:
        return "FAIL"

    for result in module_results.values():
        if result.get("status") == "FAIL":
            return "FAIL"

    for result in module_results.values():
        if result.get("status") == "REVIEW":
            return "REVIEW"

    return "PASS"


def _validate_isolate_identifiers(
    data: pd.DataFrame,
    isolate_col: str,
) -> dict[str, Any]:
    """Validate isolate identifiers for isolate-level datasets."""

    if isolate_col not in data.columns:
        return {
            "status": "FAIL",
            "error": f"Missing isolate identifier column: {isolate_col}",
        }

    result = validate_identifiers(
        data[isolate_col].tolist()
    )

    status = (
        "PASS"
        if result["is_valid"]
        else "FAIL"
    )

    return {
        "status": status,
        "results": result,
    }


def _screen_ec50_by_fungicide(
    data: pd.DataFrame,
    value_col: str,
    fungicide_col: str,
    unit: str,
    expected_unit: str | None,
) -> dict[str, Any]:
    """Screen EC50 measurements independently for each fungicide."""

    results = {}

    for fungicide, group in data.groupby(
        fungicide_col,
        sort=True,
        dropna=False,
    ):
        fungicide_label = (
            str(fungicide)
            if pd.notna(fungicide)
            else "<MISSING>"
        )

        results[fungicide_label] = validate_ec50_values(
            values=group[value_col],
            unit=unit,
            expected_unit=expected_unit,
        )

    structural_failures = [
        fungicide
        for fungicide, result in results.items()
        if not result["is_structurally_valid"]
    ]

    review_flags = [
        fungicide
        for fungicide, result in results.items()
        if result["requires_review"]
    ]

    status = (
        "FAIL"
        if structural_failures
        else "REVIEW"
        if review_flags
        else "PASS"
    )

    return {
        "status": status,
        "fungicide_count": len(results),
        "results": results,
        "structural_failures": structural_failures,
        "fungicides_requiring_review": review_flags,
    }


def run_study_qc(
    data: pd.DataFrame,
    config: StudyConfig,
    data_level: str,
    *,
    expected_unit: str | None = None,
    isolate_col: str | None = None,
    fungicide_col: str = "fungicide",
    replicate_col: str = "replicate",
    value_col: str | None = None,
    concentration_col: str = "concentration",
    response_col: str = "growth",
    expected_replicates: int = 3,
    min_treatment_levels: int = 4,
    max_replicate_cv_percent: float | None = None,
    max_control_cv_percent: float | None = None,
) -> dict[str, Any]:
    """Run the QC modules appropriate for a study data level."""

    if not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a pandas DataFrame.")

    if data.empty:
        raise ValueError("data must contain at least one observation.")

    if not isinstance(config, StudyConfig):
        raise TypeError("config must be a StudyConfig instance.")

    if data_level not in SUPPORTED_DATA_LEVELS:
        raise ValueError(
            f"Unsupported data_level: {data_level}. "
            f"Supported levels: {sorted(SUPPORTED_DATA_LEVELS)}"
        )

    config_result = config.validate()

    if not config_result["is_valid"]:
        raise ValueError(
            f"Invalid study configuration: {config_result['errors']}"
        )

    configured_isolate_col = (
        isolate_col or config.isolate_id_column
    )

    report = {
        "study_id": config.study_id,
        "pathogen": config.pathogen,
        "data_level": data_level,
        "overall_status": None,
        "modules_executed": [],
        "modules_skipped": {},
        "results": {},
        "warnings": [],
        "errors": [],
        "dataset_summary": {
            "row_count": int(len(data)),
            "column_count": int(len(data.columns)),
            "columns": list(data.columns),
        },
    }

    module_results = {}

    # ---------------------------------------------------------
    # ISOLATE-LEVEL EC50 DATA
    # ---------------------------------------------------------

    if data_level == "isolate_level":

        report["modules_skipped"] = {
            "replicate_qc": (
                "Isolate-level data do not contain individual "
                "replicate observations."
            ),
            "experimental_design_qc": (
                "Raw concentration-response measurements "
                "are not being supplied."
            ),
        }

        # Validate isolate identifiers.
        identifier_result = _validate_isolate_identifiers(
            data=data,
            isolate_col=configured_isolate_col,
        )

        module_results["identifier_qc"] = identifier_result

        report["modules_executed"].append(
            "identifier_qc"
        )

        if identifier_result["status"] == "FAIL":
            report["errors"].append(
                "Isolate identifier validation failed: "
                f"{identifier_result.get('error', identifier_result.get('results'))}"
            )

        # Validate configured phenotype columns.
        missing_columns = [
            column
            for column in config.phenotype_columns
            if column not in data.columns
        ]

        if missing_columns:
            report["errors"].append(
                f"Missing configured phenotype columns: {missing_columns}"
            )

            module_results["phenotype_qc"] = {
                "status": "FAIL",
                "error": (
                    "Missing configured phenotype columns: "
                    f"{missing_columns}"
                ),
            }

        elif config.phenotype_unit is None:
            report["errors"].append(
                "phenotype_unit must be specified for EC50 QC."
            )

            module_results["phenotype_qc"] = {
                "status": "FAIL",
                "error": "phenotype_unit is missing.",
            }

        else:
            phenotype_results = {}

            for column in config.phenotype_columns:

                result = validate_ec50_values(
                    values=data[column],
                    unit=config.phenotype_unit,
                    expected_unit=expected_unit,
                )

                phenotype_results[column] = result

            structural_failures = [
                column
                for column, result in phenotype_results.items()
                if not result["is_structurally_valid"]
            ]

            review_flags = [
                column
                for column, result in phenotype_results.items()
                if result["requires_review"]
            ]

            module_status = (
                "FAIL"
                if structural_failures
                else "REVIEW"
                if review_flags
                else "PASS"
            )

            module_results["phenotype_qc"] = {
                "status": module_status,
                "results": phenotype_results,
            }

            report["modules_executed"].append(
                "phenotype_qc"
            )

            if structural_failures:
                report["errors"].append(
                    "Phenotype structural validation failed for: "
                    f"{structural_failures}"
                )

            if review_flags:
                report["warnings"].append(
                    "Phenotype measurements require review for: "
                    f"{review_flags}"
                )

    # ---------------------------------------------------------
    # REPLICATE-LEVEL EC50 DATA
    # ---------------------------------------------------------

    elif data_level == "replicate_level":

        report["modules_skipped"] = {
            "experimental_design_qc": (
                "Raw growth measurements across fungicide "
                "concentrations were not supplied."
            ),
            "identifier_qc": (
                "Repeated isolate identifiers are expected "
                "in replicate-level data."
            ),
        }

        selected_value_col = value_col

        if selected_value_col is None:

            if len(config.phenotype_columns) == 1:
                selected_value_col = config.phenotype_columns[0]

            else:
                report["errors"].append(
                    "value_col must be supplied when the study "
                    "has multiple phenotype columns."
                )

        required_columns = [
            configured_isolate_col,
            fungicide_col,
            replicate_col,
            selected_value_col,
        ]

        missing_columns = [
            column
            for column in required_columns
            if column is not None and column not in data.columns
        ]

        if missing_columns:
            report["errors"].append(
                f"Missing replicate-level columns: {missing_columns}"
            )

        elif config.phenotype_unit is None:
            report["errors"].append(
                "phenotype_unit must be specified for EC50 QC."
            )

        elif not report["errors"]:

            # Screen each fungicide separately.
            phenotype_result = _screen_ec50_by_fungicide(
                data=data,
                value_col=selected_value_col,
                fungicide_col=fungicide_col,
                unit=config.phenotype_unit,
                expected_unit=expected_unit,
            )

            module_results["phenotype_qc"] = phenotype_result

            report["modules_executed"].append(
                "phenotype_qc"
            )

            if phenotype_result["status"] == "FAIL":
                report["errors"].append(
                    "Replicate-level EC50 phenotype validation "
                    "failed for one or more fungicides."
                )

            elif phenotype_result["status"] == "REVIEW":
                report["warnings"].append(
                    "Replicate-level EC50 measurements require "
                    "fungicide-specific review."
                )

            try:
                replicate_summary = summarize_replicates(
                    data=data,
                    isolate_col=configured_isolate_col,
                    fungicide_col=fungicide_col,
                    replicate_col=replicate_col,
                    value_col=selected_value_col,
                    expected_replicates=expected_replicates,
                    max_cv_percent=max_replicate_cv_percent,
                )

                invalid_replicates = bool(
                    (replicate_summary["invalid_count"] > 0).any()
                )

                requires_review = bool(
                    replicate_summary["requires_review"].any()
                )

                replicate_status = (
                    "FAIL"
                    if invalid_replicates
                    else "REVIEW"
                    if requires_review
                    else "PASS"
                )

                module_results["replicate_qc"] = {
                    "status": replicate_status,
                    "results": _dataframe_records(
                        replicate_summary
                    ),
                }

                report["modules_executed"].append(
                    "replicate_qc"
                )

                if invalid_replicates:
                    report["errors"].append(
                        "Invalid EC50 measurements were detected "
                        "during replicate QC."
                    )

                elif requires_review:
                    report["warnings"].append(
                        "One or more isolate-fungicide groups "
                        "require replicate-level review."
                    )

            except (ValueError, TypeError) as error:
                report["errors"].append(
                    f"Replicate QC failed: {error}"
                )

                module_results["replicate_qc"] = {
                    "status": "FAIL",
                    "error": str(error),
                }

    # ---------------------------------------------------------
    # RAW DOSE-RESPONSE DATA
    # ---------------------------------------------------------

    elif data_level == "dose_response":

        report["modules_skipped"] = {
            "phenotype_qc": (
                "Raw growth measurements are not isolate-level "
                "EC50 estimates."
            ),
            "replicate_qc": (
                "Raw dose-response observations are assessed "
                "by experimental-design QC rather than EC50 "
                "replicate summarization."
            ),
            "identifier_qc": (
                "Repeated isolate identifiers are expected "
                "in dose-response data."
            ),
        }

        try:
            design_summary = audit_dose_response_design(
                data=data,
                isolate_col=configured_isolate_col,
                fungicide_col=fungicide_col,
                concentration_col=concentration_col,
                replicate_col=replicate_col,
                response_col=response_col,
                expected_replicates=expected_replicates,
                min_treatment_levels=min_treatment_levels,
                max_control_cv_percent=max_control_cv_percent,
            )

            structural_failure = bool(
                (
                    (design_summary["invalid_concentration_count"] > 0)
                    | (design_summary["invalid_response_count"] > 0)
                    | (design_summary["duplicate_observation_count"] > 0)
                ).any()
            )

            requires_review = bool(
                design_summary["requires_review"].any()
            )

            module_status = (
                "FAIL"
                if structural_failure
                else "REVIEW"
                if requires_review
                else "PASS"
            )

            module_results["experimental_design_qc"] = {
                "status": module_status,
                "results": _dataframe_records(design_summary),
            }

            report["modules_executed"].append(
                "experimental_design_qc"
            )

            if structural_failure:
                report["errors"].append(
                    "Invalid concentration, response, or duplicate "
                    "experimental observations were detected."
                )

            elif requires_review:
                report["warnings"].append(
                    "One or more isolate-fungicide experimental "
                    "groups require design or control review."
                )

        except (ValueError, TypeError) as error:
            report["errors"].append(
                f"Experimental-design QC failed: {error}"
            )

            module_results["experimental_design_qc"] = {
                "status": "FAIL",
                "error": str(error),
            }

    report["results"] = module_results

    report["overall_status"] = _determine_status(
        module_results=module_results,
        errors=report["errors"],
    )

    return report
