# PhyloResist AI

PhyloResist AI is a research project for reproducible, explainable fungicide-resistance analysis across fungal plant pathogens. The current starter study uses *Fusarium graminearum* isolates and two fungicides, propiconazole and tebuconazole. This is the first dataset and pipeline exercise, not a limit on the platform's long-term scope.

The multi-pathogen direction, modular architecture, shared study data contract, and staged roadmap are documented in [docs/multi_pathogen_roadmap.md](docs/multi_pathogen_roadmap.md). The current *F. graminearum* workflow remains the first study module; it does not establish cross-pathogen model performance.

The current second-pathogen candidate is *Zymoseptoria tritici*. Its public European panel has paired phenotypes and sequencing data, but a ready-to-use full genotype matrix still needs to be confirmed. See the [dataset feasibility review](docs/second_pathogen_dataset_review.md) before adopting or downloading the cohort.

For future datasets, [data_input_requirements.md](docs/data_input_requirements.md) explains what users need to provide and when a resistance category is scientifically justified. Blank intake templates are in `data/templates/`. The current project does not yet provide a validated diagnostic label from genotype data alone.

## Current study and phenotype

The project includes 152 isolates, 4,857 genotype variants, and phenotype estimates for each fungicide. The study author confirmed that `trait1` is the isolate-level least-squares mean EC50 in millimolar (mM), averaged across three experimental replications. Values are used on their supplied scale; the analysis scripts do not apply a log transformation. Some replicate curve estimates were excluded by the study authors when other model parameters indicated a poor fit.

The annotation workbook is stored under `data/reference/`. Its provenance and scaffold-name handling are described in [data/reference/README.md](data/reference/README.md).
The study supplement and its extracted 152-isolate table are stored under `data/reference/`; provenance and matching checks are documented in [data/reference/README.md](data/reference/README.md). The blank metadata template is for future datasets; use `isolate_metadata.csv` for this study.

## Repository layout

```text
data/
  raw/                 Original genotype and phenotype CSV files
  reference/           Gene annotation workbook, study supplements, isolate metadata, and provenance notes
notebooks/             Numbered inspection and analysis scripts
results/
  reports/              Data-integrity and GWAS cross-check reports
  baseline_modeling/    Exploratory prediction metrics and predictions
  gwas/mrmlm/           mrMLM results, candidate annotations, and review
```

The scripts read the original files under `data/raw/`. Do not edit those source files in place; keep any cleaned or transformed data as separately named derived outputs.

## Reproduce the workflow

Run commands from the repository root in the Windows Anaconda Prompt or Miniforge Prompt. Activate the `phyloresist` environment first:

```bat
conda activate phyloresist
cd /d C:\Users\Vikas\PhyloResist-AI
```

Python scripts need `pandas`, `numpy`, `scipy`, and `scikit-learn`:

```bat
conda install pandas numpy scipy scikit-learn
```

Run the steps in order:

```bat
python notebooks\02_dataset_inspection.py
python notebooks\03_data_integrity.py
python notebooks\04_baseline_modeling.py
```

The GWAS step requires R and the R package `mrMLM`. `Rscript` must be available on `PATH` in the terminal. From the repository root, run:

```bat
Rscript notebooks\05_mrmlm_gwas.R
```

After the GWAS results exist, generate candidate-gene annotations, local LD summaries, and a reproducible cross-check:

```bat
python notebooks\06_post_gwas_candidate_review.py
```

The candidate review also writes `candidate_evidence_matrix.csv` and a short guide under `results/gwas/mrmlm/review/`. Prepare the source-backed isolate metadata, then run the exploratory population-structure diagnostic:

```bat
python notebooks\08_prepare_isolate_metadata.py
python notebooks\07_population_structure.py
python notebooks\09_geographic_holdout_modeling.py
```

Create an evidence-separated shortlist of the study-named candidate genes after the GWAS review:

```bat
python notebooks\10_candidate_shortlist.py
python notebooks\11_external_lead_marker_check.py
python notebooks\12_freeze_validation_marker_manifest.py
```

Map the current study to the shared intake format and generate an intake QC report (derived outputs only):

```bat
python notebooks\13_prepare_study_intake.py
python notebooks\14_validate_study_intake.py
python notebooks\15_prepare_ztritici_phenotype_intake.py
```

For a future dataset already mapped to the canonical tables and accompanied by its own `study_config.json`, validate it with `python notebooks\14_validate_study_intake.py --intake-dir data\processed\<study_folder>`.

This produces a state-colored PC1-versus-PC2 SVG plot, PCA scores with joined metadata, variance summaries, and descriptive pairwise allele-sharing outputs under `results/population_structure/`. Open `pca_pc1_pc2.svg` in a web browser; hover over a point to see its isolate ID, state, and population label. The pairwise values are not formal kinship estimates, and PCA alone does not assign biological population labels.

The grouped modeling script holds out one source state at a time and tunes Ridge only within the remaining state groups. Its results are an exploratory geographic-transfer check within this study, not independent population validation. Small state groups make some fold metrics unstable.

If Windows says `Rscript` is not recognized, run the R command in an R/RStudio terminal or install/configure R so its `bin` directory is on `PATH`. Miniforge provides the Python environment; it does not automatically install R.

## Main outputs

- `results/reports/integrity_report.txt` and `integrity_checks.csv`: input, isolate-ID, genotype, metadata, and phenotype checks.
- `results/baseline_modeling/`: repeated cross-validation metrics and isolate-level predictions for the exploratory Ridge baseline.
- `results/gwas/mrmlm/`: per-fungicide mrMLM output and an input manifest.
- `results/reports/gwas_results_crosscheck.md`: method counts, unique/shared markers, reported lead-locus checks, and phenotype scale summaries.
- `results/gwas/mrmlm/candidate_loci_annotation.csv`: annotated genes near the reported GWAS markers.
- `results/gwas/mrmlm/review/`: candidate prioritization and local LD review tables and reports.
- `results/gwas/mrmlm/review/candidate_shortlist.md` and `candidate_shortlist.csv`: candidate follow-up tiers with association, proximity, and local-LD evidence kept distinct.
- `results/gwas/mrmlm/review/external_validation_marker_manifest.md` and `.csv`: frozen candidate and benchmark markers with discovery-panel allele summaries, ready for external coordinate/allele harmonization.
- `results/reports/isolate_metadata_match_report.md`: isolate ID matching and source metadata coverage.
- `results/reports/external_validation_readiness.md`: candidate external panels, known overlap risks, and a pre-specified validation sequence.
- `data/processed/fg_graminearum_2025/`: canonical isolate manifest, long-format phenotype table, genotype schema manifest, and numeric checks generated from the current source files.
- `results/reports/study_intake_qc.md`: ID reconciliation, genotype coding, phenotype scale, and genome-build caveats for the shared intake mapping.
- `results/reports/study_intake_validation.md` and `study_intake_validation_checks.csv`: second-pass validation of canonical intake files against each other and the source genotype matrix.
- `results/reports/genome_reference_audit.md`: observed PH-1 sequence accessions and evidence for a likely RRes v4.0 reference family; exact study build and allele orientation remain unconfirmed.
- `data/templates/study_config_template.json`: starting point for study-specific validation expectations.
- `data/processed/z_tritici_european_panel_2025/`: derived phenotype and metadata intake for the public *Z. tritici* European panel. It is not genotype-ready.
- `results/reports/ztritici_european_phenotype_intake.md`: EC50 coverage and current limits for the second-pathogen phenotype intake.
- `results/reports/north_dakota_lead_marker_crosscheck.md` and `.csv`: descriptive check of the published North Dakota tebuconazole lead coordinate in the current panel; this is not independent validation.
- `data/reference/north_dakota_gwas_supplementary_*.pdf`: source-preserved North Dakota study supplements; Table S1 helps check isolate aliases and provenance.
- `results/population_structure/`: exploratory state-colored PCA plot and scores, variance summaries, and pairwise allele-sharing diagnostics.
- `results/geographic_holdout_validation/`: leave-one-state-out fold metrics, predictions, and interpretation notes.

These outputs can be regenerated by rerunning the scripts. The post-GWAS review script expects the mrMLM final-result CSVs to already be present.

## How to interpret the current results

- The baseline prediction is exploratory internal validation on a small, related isolate dataset. It is not deployment performance or evidence of generalization to a new population.
- GWAS associations and nearby genes are candidate signals. A nearby annotated gene is not thereby shown to cause resistance.
- The local LD scan is a first-pass analysis of the available genotype markers. It does not reproduce all manually selected regions in the study and cannot establish causality.
- The two fungicide datasets share the same isolates. Comparisons between fungicides are not independent replications.
- The current project has one pathogen and one study population. Do not present current models as multi-pathogen models.

## Next research milestones

1. Define the shared, versioned study and phenotype data contract, then map the current *F. graminearum* dataset to it without modifying source files.
2. Reconcile candidate loci and gene interpretations with the peer-reviewed article and its final supplementary data when available.
3. Obtain and harmonize independent tebuconazole phenotype/genotype data; resolve isolate overlap before testing the frozen candidate shortlist. The current validation-readiness report documents the available cohorts and requirements.
4. Add another pathogen when a well-documented paired genotype/phenotype cohort is ready, then validate it independently before attempting cross-pathogen prediction.

Population labels in the supplement are highly imbalanced (147 NA1, four unassigned, one NA2), so they do not support population-held-out validation. State-based holdout assesses exploratory geographic transfer within this dataset; tiny state groups and the absence of an independent external panel limit its interpretation.

## Data provenance

The genotype, phenotype, and reference annotation inputs originate from the study materials documented in [`data/reference/README.md`](data/reference/README.md). Study interpretation of `trait1` was confirmed directly by the study author. See the generated reports and the study's source paper/repository for detailed provenance and methods.
