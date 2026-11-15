# Second-pathogen dataset review

## Candidate reviewed

**Pathogen:** *Zymoseptoria tritici* (wheat pathogen; Septoria tritici blotch)

**Study:** Puccetti et al., “A large European diversity panel reveals complex azole fungicide resistance gains of a major wheat pathogen,” *mBio* (published November 2025).

## Data readiness

This is a strong candidate for a future second pathogen module because it combines a large, geographically and temporally broad isolate collection with both fungicide phenotypes and genome sequencing. The public isolate table has now been inspected directly:

- Table S1 contains 2,156 isolates across 63 country/location labels, collected from 1979–2019. All are included in the study PCA; 1,420 were included in SNP GWAS.
- It includes per-isolate EC50 columns for five DMIs: cyproconazole (1,326 values), epoxiconazole (1,061), prochloraz (1,236), tebuconazole (1,122), and prothioconazole (1,324). It also includes relative-growth assays, binary growth calls, country, year, coordinates, and analysis-membership flags.
- The EC50 assay used liquid culture and reports fungicide concentrations in mg/L. The authors calculated EC50 from non-linear curve fits after six days, using two technical and two biological replicates. This is a different assay and unit from the current *F. graminearum* study, whose EC50 values are in mM; values must remain separate.
- Each isolate is linked to a BioProject identifier; the table contains 151 distinct BioProjects. The largest are PRJNA1160445 (647 isolates), PRJNA870233 (497), and PRJNA1159698 (487).
- The Zenodo supplement includes Table S1, Tables S1–S14, and supplementary association outputs. It does not include a full genotype matrix or VCF.

Sources: [peer-reviewed study](https://journals.asm.org/doi/10.1128/mbio.02766-25), [Zenodo supplementary record](https://zenodo.org/records/17063534), and [PubMed record](https://pubmed.ncbi.nlm.nih.gov/41283675/).

## Important limitation

The public supplement provides per-isolate phenotypes and BioProject-level accessions, plus association results. It does not include a ready-to-use genome-wide genotype matrix or VCF for the full European panel. The study generated sequencing data for 1,134 newly sequenced strains, with 286 previously sequenced strains also included in the European panel. Reproducing a genotype matrix would require resolving isolate-to-run mappings within 151 BioProjects, harmonizing prior sequence data, selecting the IPO323 reference build, and applying the published filtering rules. The SRA reads are not directly usable as input to the current PhyloResist scripts.

The supplementary archive is 219.5 MB; raw sequencing data require separate SRA downloads. The paper reports alignment to the *Z. tritici* IPO323 reference and 472,103 retained biallelic SNPs for its European-panel GWAS. Before downloading raw reads, check whether the authors can share their final filtered VCF or genotype matrix; that would avoid a large and error-prone reconstruction task.

## Recommendation

Treat this as the **preferred second-pathogen candidate**. Its phenotype intake can be prepared now from Table S1, as a separate study module. Do not start a full GWAS or train a cross-pathogen predictor yet. First ask the authors for their final filtered VCF/genotype matrix and a sample-to-sequencing-run crosswalk. If these cannot be obtained, decide whether there is enough time and compute capacity for a 151-BioProject raw-read harmonization and variant-calling pipeline.

Keep the dataset as a separate *Z. tritici* study module. Its fungicides, assays, resistance labels, reference genome, population structure, and phenotype scales must remain pathogen/study specific. Do not pool its values or train a single *F. graminearum* + *Z. tritici* model without a separate harmonization and validation design.

## How it compares with the current study

| Feature | Current *F. graminearum* study | Candidate *Z. tritici* study |
|---|---|---|
| Isolates | 152 | 1,394 European strains |
| Phenotypes | Propiconazole and tebuconazole EC50, author-confirmed in mM | EC50 for up to five DMIs plus relative growth; assay-specific definitions must be extracted |
| Genotype availability | Ready variant-by-isolate matrix, 4,857 variants | Public sequence reads/accessions reported; ready full-panel genotype matrix not established yet |
| Immediate role | Existing pipeline reference module | Candidate for the next pathogen module after data-intake feasibility check |

## Decision gate before adoption

Proceed only after a feasibility inventory confirms:

1. The isolate/phenotype table can be downloaded and interpreted without ambiguous units or IDs. **Completed:** Table S1 supplies this, and EC50 values are reported in mg/L.
2. A final genotype matrix/VCF and isolate-to-sequencing-run crosswalk are obtained, or a reproducible raw-read workflow is within time and compute limits.
3. IPO323 reference build/version, variant representation, and reuse terms are documented.
4. The chosen analysis can be validated within *Z. tritici* without mixing it with the current *F. graminearum* cohort.
