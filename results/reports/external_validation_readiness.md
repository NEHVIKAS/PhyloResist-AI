# External Validation Readiness

## Current status

The current 152-isolate panel supports exploratory model development and geographic holdout only. It is predominantly one source population, and the geographic holdout did not improve pooled RMSE over the training-mean baseline for either fungicide. Candidate genes and prediction results therefore still need validation in a genuinely independent panel.

## External panels located

| Panel | What is reported | What it could validate | Main limitations / next check |
|---|---|---|---|
| North Dakota, 183 isolates (1981–2013) | An independent published GWAS assayed tebuconazole and prothioconazole and generated 104,770 GBS SNPs, with 11,965 retained for GWAS. It reports one tebuconazole marker association. | Potential candidate-level genetic replication for tebuconazole, if the marker/genotype data and isolate-level phenotype values can be obtained and harmonized. It does not cover propiconazole. | Supplementary Table S1 confirms current isolate 23262 / Fg08-007 against Fg08_7 by county, year, and chemotype. Current isolate 23269 / Fg08-018 resembles Fg08_18 by alias/year/chemotype, but its current county (Steele) conflicts with Table S1 (Griggs). Resolve that discrepancy and check all aliases before calling the panel independent. The published lead is reported as FGRAPHPH1_01G16325; the corresponding local annotation ID is spelled FGRAMPH1_01G16325. This gene is not among the current shortlist candidates, although its reported marker coordinate is present in this panel (see cross-check below). Harmonize genome coordinates and assay/trait definitions before comparing loci or effects. |
| Huang-Huai-Hai, China, 165 isolates (2019–2023) | A published baseline study reports isolate-level tebuconazole EC50 values and supplementary data. It also reports laboratory-selected mutants. | Independent phenotype-distribution context for tebuconazole. | No population-wide genotype matrix for the 165 field isolates was identified in the article materials reviewed, so it cannot currently validate the shortlisted GWAS loci. Do not pool values directly with this project until units and assay protocols are reconciled. |

Sources: [North Dakota fungicide-sensitivity GWAS](https://pubmed.ncbi.nlm.nih.gov/38079375/) and [publisher article](https://apsjournals.apsnet.org/doi/10.1094/PHYTO-05-23-0180-KC); [China tebuconazole baseline study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11337888/). The isolate overlap warning is also visible in this repository's [matched isolate metadata](../reports/isolate_metadata_match_report.md) and `data/reference/isolate_metadata.csv`.

## Data acquisition check

The supplied [Supplementary Figures S1-S3](../../data/reference/north_dakota_gwas_supplementary_figures.pdf) show grouped phenotype comparisons for tebuconazole and prothioconazole, while [Supplementary Table S1](../../data/reference/north_dakota_gwas_supplementary_table_s1.pdf) supplies isolate IDs, collection metadata, and chemotypes. The figures do not provide a numeric per-isolate EC50 table, and these files do not include the genome-wide genotype matrix or raw read archive needed for a locus-level reanalysis. The author-provided genotype/phenotype files and crosswalk are still needed.

## Published North Dakota lead coordinate in the current panel

The North Dakota paper's tebuconazole marker `SFG3_297186` is reported in gene `FGRAPHPH1_01G16325`. The current marker set contains `SHG970334.1_297186`, and the local annotation workbook places it inside `FGRAMPH1_01G16325` at `HG970334:296989-297434`. The published and local gene IDs differ in spelling, so confirm the gene alias and genome-build mapping before treating these as the same locus. The marker is not among the current tebuconazole mrMLM final-result markers.

The current marker IDs point to PH-1 chromosome sequence accessions HG970332.2, HG970333.1, HG970334.1, and HG970335.2, which are consistent with the published RRes v4.0 reference accession family. The study repository does not explicitly confirm the assembly version or allele orientation, so these remain author-confirmation items; see [genome reference audit](genome_reference_audit.md).

In the current panel, the unadjusted mean trait1 values were 0.019949 mM for C (n=25) and 0.016827 mM for T (n=127), with Welch p=0.1417. Excluding both overlap candidates leaves the means effectively unchanged (C n=25; T n=125; p=0.1407). These are descriptive, structure-unadjusted summaries in the discovery panel, not an independent replication test; reference-allele harmonization is still required before direction can be compared. See [marker cross-check report](north_dakota_lead_marker_crosscheck.md) and its CSV.

The primary tebuconazole test set is now recorded in [external_validation_marker_manifest.md](../gwas/mrmlm/review/external_validation_marker_manifest.md): the three exact markers shared across the current propiconazole and tebuconazole analyses, plus the published North Dakota lead as a separate benchmark. Current-panel allele directions exclude one confirmed overlap and one unresolved alias conflict. Keep the marker set fixed unless a scientific reason to amend it is documented before external phenotypes are inspected.

## Recommended validation sequence

1. **Request the North Dakota data needed for reanalysis.** Obtain isolate identifiers/alias crosswalk, per-isolate tebuconazole EC50 values, and the GBS genotype matrix or raw reads. Confirm the reference build and variant coordinates.
2. **Deduplicate before analysis.** Match the 183-isolate list against `isolate_id`, `alternate_id`, collection year, and source-publication fields in `data/reference/isolate_metadata.csv`. Exclude exact overlaps with the current panel. Keep unresolved aliases out of the validation set until resolved.
3. **Freeze the current candidate test set.** Predefine the three Tier A genes/markers in `candidate_shortlist.csv`; record the marker IDs, allele coding, and effect direction from the training panel before opening external phenotypes. Do not select or retune markers using the external panel.
4. **Harmonize the phenotype and genome coordinates.** The current `trait1` is an isolate-level least-squares mean EC50 in mM. The North Dakota study used a microtiter assay and also reports EC50, but differences in assay protocol, concentration units, estimation, and reference assembly must be resolved before comparing effect sizes. Analyze tebuconazole only for this panel.
5. **Run the independent marker replication.** For each frozen marker, report allele counts, missingness, effect estimate and confidence interval, and whether the direction agrees with the current-panel estimate. Adjust the three planned marker tests for multiplicity. If exact markers are unavailable, define a coordinate-based locus replication rule before testing.
6. **Use the China cohort only for a separate phenotype check unless genotypes are obtained.** Compare distributions and ranges only after unit and protocol harmonization; do not treat phenotype-only agreement as validation of a genetic marker or predictive model.

## Decision gate

Do not describe a candidate as externally replicated until isolate overlap is excluded, marker coordinates and alleles are harmonized, and the frozen marker's effect is tested in independent genotypes with matched tebuconazole phenotype data. The 183-isolate North Dakota study is the strongest next genetic-validation target found, but the necessary data and overlap resolution are not yet in this repository. No independent test performance is reported here.

## Current-panel overlap to resolve

The North Dakota GWAS describes its 183-isolate panel as including isolates used by Puri and Zhong (2010). Supplementary Table S1 gives the following cross-check:

| Current genotype ID | Current alternate ID / county | Supplement Table S1 alias / county | Status |
|---|---|---|---|
| 23262 | Fg08-007 / Foster | Fg08_7 / Foster County | Strong match: alias, year, county, and 15ADON chemotype agree. |
| 23269 | Fg08-018 / Steele | Fg08_18 / Griggs County | Unresolved: alias/year/chemotype resemble each other, but county conflicts. |

Exclude both conservatively until the study's isolate crosswalk or author confirms the second alias. Check all other accession names against the full crosswalk before designating any samples independent.
