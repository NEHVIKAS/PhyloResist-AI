# GWAS Results Cross-Check

Generated from the repository's final mrMLM result CSVs by `06_post_gwas_candidate_review.py`.

## Method-call counts and unique markers

### Propiconazole

- FASTmrMLM: 1 calls
- ISIS EM-BLASSO: 30 calls
- mrMLM: 2 calls
- pKWmEB: 22 calls
- Unique marker IDs after cross-method deduplication: 48
- Markers detected by more than one method: 6

### Tebuconazole

- ISIS EM-BLASSO: 22 calls
- pKWmEB: 17 calls
- Unique marker IDs after cross-method deduplication: 37
- Markers detected by more than one method: 2

## Cross-fungicide comparison

- Shared exact marker IDs: 13
- Distinct marker IDs across both fungicides: 72
- The supplied preprint's Results section reports 13 shared QTNs while its Abstract says 12; these outputs match the Results count.
- These are mrMLM output counts, not a new significance test or external validation.

## Reported lead loci

- propiconazole: SHG970334.1_6640556 on chromosome 3 has r²=10.097% under mrMLM; the supplied Results section reports ~10% variance in the Results text.
- tebuconazole: SHG970333.1_8322036 on chromosome 2 has r²=3.0555% under pKWmEB; the supplied Results section reports 3.06% variance in the Results text.

## Phenotype scale and range

The study author confirmed `trait1` is an isolate-level least-squares mean EC50 in mM. The GWAS uses the stored values without log transformation. Converted summaries below are for comparison only.
- propiconazole: n=152, mean=34.93 µM, range=8.62–100.74 µM.
- tebuconazole: n=152, mean=17.34 µM, range=4.58–49.37 µM.
