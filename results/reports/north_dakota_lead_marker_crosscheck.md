# North Dakota Tebuconazole Lead Marker Cross-Check

The 2024 North Dakota GWAS reports marker `SFG3_297186` for tebuconazole and annotates it to `FGRAPHPH1_01G16325`. The current genotype file contains `SHG970334.1_297186`; its normalized scaffold and position match the current annotation workbook's `FGRAMPH1_01G16325` interval. The spelling of the gene ID differs between the published table and the local annotation, so this is recorded as a coordinate-supported candidate mapping, not an assumed gene-ID alias.

The variant is inside the local annotated gene interval (HG970334:296989-297434). Exact REF/ALT harmonization against the external assembly remains pending.

The marker was not listed among the current panel's tebuconazole mrMLM final-result markers.

## Descriptive allele summaries in the current panel

| Set | n | Allele 1 (n; mean mM) | Allele 2 (n; mean mM) | Difference (mM) | Unadjusted Welch p |
|---|---:|---|---|---:|---:|
| all_current_isolates | 152 | C (25; 0.019949) | T (127; 0.016827) | 0.003121 | 0.1417 |
| exclude_overlap_candidates | 150 | C (25; 0.019949) | T (125; 0.016816) | 0.003133 | 0.1407 |

The sensitivity check conservatively excludes IDs 23262 and 23269. Supplementary Table S1 confirms 23262 / Fg08-007 against Fg08_7 by county, year, and chemotype. The 23269 / Fg08-018 alias resembles Fg08_18, but the current county (Steele) conflicts with Table S1 (Griggs). Obtain the external crosswalk to resolve the second alias and check for any additional overlaps.

These are unadjusted summaries in the current 152-isolate panel, not a replication test. They do not correct for population structure, relatedness, selection, or testing multiple markers. Do not compare effect directions until the published reference allele is mapped to the current C/T calls. The external paper's published association remains the evidence for its original cohort; this report only checks whether that coordinate is represented and how it behaves descriptively here.

No independent data were available for this check. The next validation step still requires the North Dakota isolate crosswalk and genotype/phenotype files.

Source: [North Dakota GWAS](https://pubmed.ncbi.nlm.nih.gov/38079375/). Full values: `north_dakota_lead_marker_crosscheck.csv`.
