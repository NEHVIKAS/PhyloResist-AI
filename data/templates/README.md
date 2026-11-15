# Intake templates

These blank CSVs describe the target shape for future study intake. They are templates, not replacements for or edits to the original files in `data/raw/`.

- `isolate_manifest_template.csv`: one row per isolate, with stable IDs and study/pathogen provenance.
- `phenotype_measurements_template.csv`: long format, one row per isolate, fungicide, trait, and observation/replicate.
- `genotype_long_template.csv`: illustrative long-format variant calls for small examples and schema discussion. Large datasets should use an efficient matrix or a supported VCF adapter rather than this expanded format.
- `study_config_template.json`: study ID, pathogen, expected fungicides, phenotype definition, and the canonical file names used by the validation step.

The first study's mapped files are derived outputs under `data/processed/fg_graminearum_2025/`. From the repository root, run `python notebooks/13_prepare_study_intake.py` to regenerate them, then `python notebooks/14_validate_study_intake.py` for the second-pass checks. For another canonical intake folder, pass `--intake-dir` to script 14 and provide its own `study_config.json`. The original inputs remain under `data/raw/`.

Unknown fields should remain blank and be documented. Do not infer phenotype units, assay meaning, resistance thresholds, reference assemblies, or allele orientation. See [the input requirements guide](../../docs/data_input_requirements.md).
