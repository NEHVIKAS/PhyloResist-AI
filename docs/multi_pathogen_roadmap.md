# PhyloResist AI: multi-pathogen roadmap and architecture

## Project direction

PhyloResist AI is intended to support fungicide-resistance research across multiple fungal plant pathogens. Each analysis must remain specific to its pathogen, fungicide, study population, assay, and genome reference. The platform can share data handling and reporting conventions, but results from different organisms or studies must not be pooled or treated as interchangeable without explicit harmonization and validation.

The current *Fusarium graminearum* study is the first dataset module. It is a useful end-to-end pipeline exercise, not evidence that the current predictors generalize to other pathogens. Its source files and completed inspection, integrity, association, and exploratory validation work remain the reference implementation while the architecture is developed.

## Design principles

1. **Keep study inputs traceable.** Raw files remain immutable. Every derived table and report records its source files, study, pathogen, fungicide, and processing choices.
2. **Represent phenotype definitions explicitly.** Store the trait name, units, assay, experimental design, transformation (if any), and source. Never infer units or apply transformations from column names.
3. **Do not mix incompatible studies by default.** Combining cohorts requires documented checks for assay comparability, phenotype scale, population structure, genome build, variant representation, and isolate overlap.
4. **Separate shared code from biological assumptions.** File loading, ID matching, QC summaries, provenance, and report generation can be shared. Genome annotation, known resistance mechanisms, allele interpretation, and model eligibility belong in pathogen-specific configuration or modules.
5. **Keep discovery and validation separate.** Freeze candidate markers and analysis choices before examining independent validation phenotypes. Report internal, geographic, and external validation distinctly.
6. **Make limitations visible.** Small or related cohorts, missing genotypes, unresolved aliases, and sparse phenotype ranges should produce explicit warnings and eligibility decisions.

## Target architecture

```text
PhyloResist AI
├── shared core
│   ├── input registry and provenance
│   ├── schema validation and isolate matching
│   ├── genotype/phenotype QC
│   ├── analysis orchestration and reports
│   └── validation and reproducibility utilities
├── pathogen modules
│   ├── pathogen identity and taxonomy
│   ├── genome assemblies and annotation mappings
│   ├── pathogen-specific marker/mechanism references
│   └── supported genotype formats and QC rules
└── study configurations
    ├── cohort and fungicide
    ├── phenotype definition and assay metadata
    ├── genome build and coordinate conventions
    ├── population/geographic metadata
    └── pre-registered analysis and validation settings
```

An analysis run should be identified by a study configuration rather than by hard-coded assumptions in a numbered script. The existing numbered scripts remain usable during the transition; refactoring should happen in small steps, preserving their outputs and documenting any changed results.

## Shared study data contract

Each future cohort should provide, where available:

- **Study record:** study ID, source citation, pathogen name and strain, cohort description, data license/access terms, and file checksums.
- **Isolate table:** one row per isolate, stable isolate ID, alternate IDs, study/cohort ID, collection date/location where permitted, host, and population labels with their sources.
- **Phenotype table:** isolate ID, fungicide/active ingredient, measured trait, value, units, assay protocol, replicate or model-summary status, uncertainty fields, and source citation.
- **Genotype record:** isolate IDs, variant matrix or sequence accessions, variant type, reference assembly/build, coordinate system, allele representation, missing-value coding, and QC history.
- **Annotation record:** annotation source/version, assembly identifier, gene identifiers and aliases, and coordinate conventions.

Missing fields should be marked as unknown or unavailable rather than filled by inference. The contract is a target for harmonized inputs; it does not imply that every published cohort supplies every field.

## Analysis layers

### Layer 1: study-level intake and QC

Validate file structure, isolate identifiers, duplicates, missingness, genotype coding, phenotype numeric ranges, metadata completeness, and build/annotation consistency. Produce a machine-readable check table and a human-readable report before modeling.

### Layer 2: within-study analyses

Run phenotype summaries, population-structure diagnostics, association analyses, candidate review, and prediction using only that study's supported design. Record decisions and keep exploratory results clearly labeled.

### Layer 3: within-pathogen cross-study validation

First reconcile isolate overlap, assay definitions, units, phenotype estimation, genome builds, marker alleles, and cohort structure. Then test frozen candidates or models in independent cohorts. If direct marker harmonization fails, define a locus-level rule before looking at validation phenotypes.

### Layer 4: cross-pathogen comparison

Compare evidence at the level that is biologically and technically supportable. Shared fungicide classes or homologous genes do not make allele effects equivalent. Cross-pathogen prediction requires explicit orthology/feature mapping and a separate validation design; it is not a default capability of the platform.

## Roadmap

### Phase A — preserve and close out the first study

- Keep the current *F. graminearum* source data and inspection workflow intact.
- Maintain the candidate shortlist, frozen marker manifest, and validation-readiness report.
- Resolve the North Dakota alias conflict and obtain per-isolate phenotype/genotype data before claiming external replication.
- Reconcile the current study's candidate interpretation against its peer-reviewed article and final supplements when available.

### Phase B — establish the reusable data contract

- Define versioned schemas for study, isolate, phenotype, genotype, and reference-build metadata.
- Add an intake/QC entry point that accepts a study configuration and writes derived outputs into a study-specific results directory.
- Keep the current scripts as a compatibility path until the new entry point reproduces their key summaries.

### Phase C — separate shared workflow from pathogen-specific settings

- Move reusable loading, ID checks, QC, reporting, and validation utilities into shared modules.
- Represent each pathogen/genome build through explicit configuration, including annotation aliases and supported genotype conventions.
- Add a small synthetic fixture for each schema so structural checks can be exercised without distributing sensitive or large study data.

### Phase D — add another pathogen only when a suitable cohort is ready

- The current preferred candidate is *Zymoseptoria tritici*, with a large public European cohort containing per-isolate fungicide phenotypes and sequence accessions. Its full genotype-matrix readiness is not established; see [the dataset review](second_pathogen_dataset_review.md).
- First audit the public isolate table, data license, sample crosswalk, and available genotype products. If only raw reads are available, decide whether a read-to-variant pipeline fits the project's time and compute limits.
- Complete an independent study intake and QC report before fitting models.
- Validate within that pathogen and population before attempting any cross-pathogen prediction.

### Phase E — comparative surveillance research

- Add cross-study and cross-pathogen analyses only when phenotype, genotype, and reference mappings support defensible comparisons.
- Report pathogen-specific performance and uncertainty alongside any comparative result.

## Immediate next work

The first Phase B deliverables are the [user input requirements](data_input_requirements.md) and blank templates in `data/templates/`. The current *F. graminearum* files have been mapped to derived intake tables, and `notebooks/14_validate_study_intake.py` reads study-specific expectations from `study_config.json` and can validate another canonical intake folder with `--intake-dir`. Next, add synthetic fixtures for schema edge cases without changing raw CSVs. The next scientific task remains the North Dakota external replication, which cannot begin until the isolate crosswalk, numeric per-isolate tebuconazole phenotypes, genotype matrix or reads, and coordinate/reference information are available.

## Current project boundary

The repository currently contains one pathogen, one primary study panel, and two fungicides. It supports exploratory within-study analyses and geographic holdout diagnostics. It does not yet provide a validated multi-pathogen predictor or external validation of the shortlisted genes.
