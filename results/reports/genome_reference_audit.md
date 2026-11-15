# Genome reference and coordinate audit

## Finding

The current genotype matrix does not declare an assembly accession or version in its metadata, and the public fungicide repository's `analysis_2025.R` does not identify the reference assembly. The exact build therefore remains **unconfirmed**.

However, the local variant identifiers carry the following sequence accession versions:

| Sequence accession | `chrom` value | Variant rows |
|---|---:|---:|
| HG970332.2 | 1 | 1,174 |
| HG970333.1 | 2 | 1,704 |
| HG970334.1 | 3 | 982 |
| HG970335.2 | 4 | 997 |

The prefix `S` in identifiers such as `SHG970332.2_12065` appears to be a local marker-label prefix; after removing it, the identifiers point to the four PH-1 chromosome accessions. The published PH-1 genome paper identifies accessions HG970330–HG970335 as the RRes v4.0 sequence set, so the current markers are consistent with that reference family. This is supporting evidence, not confirmation that the study's variant calls and annotation coordinates were all generated against that exact version.

## Why the build warning remains

- The study repository includes the genotype CSV, annotation workbook, and an mrMLM analysis script, but the script does not state the genome assembly or variant-calling reference.
- Sequence accession labels can establish which chromosome sequence identifiers appear in the marker names, but do not by themselves establish the complete pipeline's reference files, mapping version, or allele orientation.
- The annotation workbook's gene IDs and the local coordinate mapping should be checked against the exact annotation version used by the authors.

Before cross-study coordinate or allele-direction comparisons, confirm with the authors that the genotype positions use the PH-1 RRes v4.0 sequence set represented by HG970332.2, HG970333.1, HG970334.1, and HG970335.2, and ask how the leading `S` in marker names and the alleles are defined. Keep the genome-build check as a warning until this is confirmed.

## Sources

- Local observation: `data/raw/genotype.csv`; summary reproduced in `data/processed/fg_graminearum_2025/genotype_input_manifest.json`.
- Study data/code repository: [files_fungicide_manuscript](https://github.com/Toomajian-laboratory/files_fungicide_manuscript); its current [analysis script](https://github.com/Toomajian-laboratory/files_fungicide_manuscript/blob/main/analysis_2025.R) does not specify an assembly.
- The PH-1 genome paper identifies HG970330–HG970335 as the RRes v4.0 sequence set: [King et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4511438/).
- NCBI identifies [HG970334.1 as *F. graminearum* PH-1 chromosome 3](https://www.ncbi.nlm.nih.gov/nuccore/HG970334.1).
