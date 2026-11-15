# Data needed to assess fungicide sensitivity

## What the platform can and cannot say

PhyloResist AI is currently a research workflow, not a validated diagnostic service. A user cannot upload any genome and receive a trustworthy susceptible/resistant label from the current project. The current *Fusarium graminearum* analysis supports exploratory comparisons and candidate follow-up; its models have not demonstrated dependable prediction in an independent cohort.

There are two different ways to assess a sample:

1. **Measured sensitivity:** provide a laboratory phenotype such as EC50, with its units, assay protocol, fungicide, isolate ID, and replicate/model-summary details. This supports a continuous sensitivity estimate or comparison. A resistant/susceptible category is justified only if a validated, pathogen- and fungicide-specific breakpoint or an explicitly justified study threshold exists.
2. **Genotype-based estimate:** provide variant/genotype data for identified isolates, matched to phenotype data used to train and validate a pathogen- and fungicide-specific model. This can estimate a phenotype only within the model's validated scope. A new pathogen, fungicide, population, assay, or genome build requires separate validation or harmonization.

Genotype data alone do not establish resistance. A marker association is not proof of causality, and the same allele can behave differently across genetic backgrounds or assay conditions.

## Minimum inputs for a research assessment

| Input | Required information | Why it is needed |
|---|---|---|
| Isolate manifest | Stable isolate ID; pathogen/species; study/cohort ID | Defines what each row/sample represents and prevents mixing organisms or cohorts. |
| Phenotype table | Isolate ID; fungicide/active ingredient; trait name; numeric value; units; assay/source | Defines the observed sensitivity and ensures it can be interpreted. For resistance categories, include the validated breakpoint and its source. |
| Genotype data for genetic analyses | Isolate IDs matching the manifest; variant calls or a clearly coded genotype matrix; reference assembly/build; chromosome/position and allele coding; missing-value convention | Enables variant QC, association or model input, and coordinate/allele harmonization. |
| Provenance | Citation or data source, file version/date, and permission/license where relevant | Makes results reproducible and determines whether the data may be redistributed. |

If the task is only to summarize measured sensitivity, genotype data are not necessary. If the task is a genetic association or genotype-based prediction, matched genotype and phenotype data are both needed.

## Strongly recommended metadata

- Collection year and location, host, lineage/population, and experimental batch, with source and missingness stated.
- Replicate count, replicate-level observations when available, uncertainty fields, and whether the submitted value is a raw measurement, a model estimate, or an average such as a least-squares mean. Preserve unclear source column labels until their definitions are verified.
- Assay medium, temperature, incubation duration, inoculum, concentration series, endpoint definition, and curve-fitting method.
- Reference genome and annotation versions, ploidy/call convention, variant-calling/QC method, and whether genotype alleles are reference-oriented.
- Known isolate aliases/accessions and a crosswalk when the same isolate appears in multiple studies.

Do not fill missing metadata by guessing from filenames, isolate IDs, or numeric ranges. Record unavailable information as unknown and explain its effect on interpretation.

## Supported input shapes for the first reusable intake

- **Isolate manifest:** one row per isolate; see `data/templates/isolate_manifest_template.csv`.
- **Phenotype measurements:** long format, one row per isolate × fungicide × trait × measurement/replicate; see `data/templates/phenotype_measurements_template.csv`. Existing study files can be mapped to this format as derived data without editing their source CSVs.
- **Genotypes:** a variant-by-isolate matrix with a metadata block describing the genome build, variant identity, allele encoding, and missing-value codes. The current study's A/C/G/T matrix is an example of one encoding. VCF can be supported by a later intake adapter after its build, ploidy, and filtering rules are specified.

Raw sequencing reads (FASTQ) are not directly suitable for the current analysis scripts. They first require a separate, validated bioinformatics workflow for quality control, alignment, variant calling, and harmonization to an agreed reference assembly.

## Resistance status and thresholds

“Resistance status” is not a universal property of a genotype or a single EC50 number. The label depends on the pathogen, fungicide, assay, and an agreed interpretation rule. Where no accepted clinical/agronomic breakpoint exists, the workflow should report the measured value and uncertainty, or use clearly labeled relative categories defined for a specific study. It must not silently invent a cutoff or convert an exploratory model output into a diagnostic label.

## What a future user workflow should look like

1. Select the pathogen, fungicide, and study question.
2. Provide a manifest and phenotype file; add genotype data for genetic analyses or prediction.
3. Declare assay definitions, units, genome build, and source information.
4. Run intake checks and resolve ID, missingness, coding, and compatibility warnings.
5. Review the QC report before association or prediction is enabled.
6. Receive measured sensitivity summaries or model estimates with uncertainty, scope, and validation status. Show a resistant/susceptible label only when a suitable threshold has been supplied and validated.

The current *F. graminearum* files are the first example for building this workflow, not a general model or clinical threshold for other pathogens.
