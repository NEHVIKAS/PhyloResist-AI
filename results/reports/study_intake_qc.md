# Study intake and QC summary

Generated: 2026-11-12T14:13:46+01:00
Study: `fg_graminearum_fungicide_study_2025` (*Fusarium graminearum*)

## Intake results

- Isolate manifest rows: 152
- Genotype isolates: 152; variants: 4,857
- Genotype coding: A, C, G, T; missing calls: 0
- Phenotype rows: 304 (two fungicides × 152 isolates)
- Isolate IDs match exactly across genotype, metadata, and both phenotype files.
- Phenotype `trait1` is represented as author-confirmed isolate least-squares mean EC50 in mM across three experimental replications.
- Trait values are retained on their supplied scale; no log transform or resistance cutoff is applied.
- Source fields `lower.CL u` and `pper.CL` are carried under source-labeled columns; their semantics are not inferred from the headings.
- Genome build and reference/alternate allele orientation remain unspecified and must be resolved before coordinate-based external harmonization.
- Variant identifiers contain sequence accessions: HG970332.2, HG970333.1, HG970334.1, HG970335.2. These are consistent with the PH-1 RRes v4.0 accession family, but the study repository does not explicitly name the build.

## Derived files

- `data/processed/fg_graminearum_2025/isolate_manifest.csv`
- `data/processed/fg_graminearum_2025/phenotype_measurements.csv`
- `data/processed/fg_graminearum_2025/genotype_input_manifest.json` (describes the original matrix without expanding or duplicating it)
- `data/processed/fg_graminearum_2025/study_config.json` (study-specific expectations used by the validator)
- `data/processed/fg_graminearum_2025/phenotype_numeric_checks.csv`

Original source CSVs were read only and were not modified.
