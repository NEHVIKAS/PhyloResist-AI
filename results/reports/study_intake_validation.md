# Study intake validation

Generated: 2026-11-12T14:13:53+01:00
Study: `fg_graminearum_fungicide_study_2025` (Fusarium graminearum)
Status: **PASS WITH WARNINGS**

## Summary

- Isolates: 152
- Phenotype rows: 304
- Fungicides: propiconazole, tebuconazole
- Failed checks: 0
- Original genotype matrix checked: checked

## Checks

| Status | Check | Details |
|---|---|---|
| PASS | Manifest required columns | Missing=[] |
| PASS | Phenotype required columns | Missing=[] |
| PASS | Unique manifest isolate IDs | Duplicate rows=0 |
| PASS | Unique phenotype isolate/fungicide/trait records | Duplicate rows=0 |
| PASS | Manifest/phenotype isolate matching | Manifest-only=0; phenotype-only=0; shared=152 |
| PASS | Configured fungicide coverage | Expected=['propiconazole', 'tebuconazole']; observed=['propiconazole', 'tebuconazole'] |
| PASS | Configured study identifier | Expected='fg_graminearum_fungicide_study_2025'; observed=['fg_graminearum_fungicide_study_2025'] |
| PASS | Configured pathogen identity | Expected='Fusarium graminearum'; observed=['Fusarium graminearum'] |
| PASS | Configured phenotype trait | Expected='trait1'; observed=['trait1'] |
| PASS | Phenotype values finite and numeric | Missing=0; non-numeric=0; non-finite=0 |
| PASS | Configured phenotype units | Expected='mM'; observed=['mM'] |
| INFO | Resistance threshold supplied | No threshold supplied; report measured EC50 values and do not assign susceptible/resistant labels. |
| PASS | Genotype source dimensions match manifest | Source=4857 variants × 152 isolates; manifest=4857 × 152 |
| PASS | Genotype source isolate IDs match manifest | Source-only=0; manifest-only=0 |
| PASS | Genotype coding matches manifest | Source=['A', 'C', 'G', 'T']; manifest=['A', 'C', 'G', 'T'] |
| WARN | Genome build identified | Genome build is unknown; resolve before cross-study coordinate harmonization. |
| INFO | Coordinate accession evidence | Observed=['HG970332.2', 'HG970333.1', 'HG970334.1', 'HG970335.2']. The observed HG970330-HG970335 accession family is consistent with the PH-1 RRes v4.0 reference, but the fungicide study repository does not explicitly confirm the build. |

This is an intake-format/QC check, not an external validation or resistance classifier.
A missing threshold means no susceptible/resistant label should be assigned.
