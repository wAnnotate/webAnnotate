Description of columns in the cmc_export.tsv file (MODIFIED)

GENE_NAME - The gene name for which the data has been curated in COSMIC. In most cases this is the accepted HGNC symbol
ACCESSION_NUMBER - The transcript identifier of the gene
ONC_TSG - Role of gene in cancer
CGC_TIER - Cancer gene census Tier (for more details see https://cancer.sanger.ac.uk/census)
MUTATION_URL - URL of mutation page on the main COSMIC site
LEGACY_MUTATION_ID - Legacy mutation identifier (COSM) that represents existing COSM mutation identifiers
Mutation CDS - The change that has occurred in the nucleotide sequence. Formatting is identical to the method used for the peptide sequence
Mutation AA - The change that has occurred in the peptide sequence. Formatting is based on the recommendations made by the Human Genome Variation Society (HGVS)
AA_MUT_START - Start position of the peptide change
AA_MUT_STOP - Stop position of the peptide change (for frameshift variants, this will be the same as AA_MUT_START)
SHARED_AA - Number of mutations seen in this amino acid position
GENOMIC_WT_ALLELE_SEQ - Wild-type/reference allele in the genomic change (on the forward strand)
GENOMIC_MUT_ALLELE_SEQ - Mutant allele in the genomic change (on the forward strand)
AA_WT_ALLELE_SEQ - Wild-type/reference amino acid in the peptide change
AA_MUT_ALLELE_SEQ - Mutant amino acid in the peptide change
Mutation Description CDS - Type of mutation at the nucloetide level
Mutation Description AA - Type of mutation at the amino acid level
ONTOLOGY_MUTATION_CODE - Sequence Ontology (SO) code for the mutation
GENOMIC_MUTATION_ID - Genomic mutation identifier (COSV) to indicate the definitive position of the variant on the genome. This identifier is trackable and stable between different versions of the release
Mutation genome position GRCh37 - The genomic coordinates of the mutation on the GRCh37 assembly
Mutation genome position GRCh38 - The genomic coordinates of the mutation on the GRCh38 assembly
COSMIC_SAMPLE_TESTED - Number of samples in COSMIC tested for this mutation
COSMIC_SAMPLE_MUTATED - Number of samples in COSMIC with this mutation
DISEASE - Diseases with > 1% samples mutated (or frequency > 0.01), where disease = Primary site(tissue) / Primary histology / Sub-histology = Samples mutated / Samples tested = Frequency
WGS_DISEASE - Same as DISEASE, but for whole-genome screen data only
EXAC_AF - Allele frequency in all ExAC samples
EXAC_AFR_AF - Adjusted Alt allele frequency in African & African American ExAC samples
EXAC_AMR_AF - Adjusted Alt allele frequency in American ExAC samples
EXAC_EAS_AF - Adjusted Alt allele frequency in East Asian ExAC samples
EXAC_FIN_AF - Adjusted Alt allele frequency in Finnish ExAC samples
EXAC_NFE_AF - Adjusted Alt allele frequency in Non-Finnish European ExAC samples
EXAC_SAS_AF - Adjusted Alt allele frequency in South Asian ExAC samples
GNOMAD_EXOMES_AF - Alternative allele frequency in all gnomAD exome samples (123,136 samples)
GNOMAD_EXOMES_AFR_AF - Alternative allele frequency in the African/African American gnomAD exome samples (7,652 samples)
GNOMAD_EXOMES_AMR_AF - Alternative allele frequency in the Latino gnomAD exome samples (16,791 samples)
GNOMAD_EXOMES_ASJ_AF - Alternative allele frequency in the Ashkenazi Jewish gnomAD exome samples (4,925 samples)
GNOMAD_EXOMES_EAS_AF - Alternative allele frequency in the East Asian gnomAD exome samples (8,624 samples)
GNOMAD_EXOMES_FIN_AF - Alternative allele frequency in the Finnish gnomAD exome samples (11,150 samples)
GNOMAD_EXOMES_NFE_AF - Alternative allele frequency in the Non-Finnish European gnomAD exome samples (55,860 samples)
GNOMAD_EXOMES_SAS_AF - Alternative allele frequency in the South Asian gnomAD exome samples (15,391 samples)
GNOMAD_EXOMES_OTH_AF - Alternative allele frequency in other gnomAD exome samples (2,743 samples)
GNOMAD_GENOMES_AF - Alternative allele frequency in the whole gnomAD genome samples (15,496 samples)
GNOMAD_GENOMES_AFR_AF - Alternative allele frequency in the African/African American gnomAD genome samples (4,368 samples)
GNOMAD_GENOMES_AMR_AF - Alternative allele frequency in the Latino gnomAD genome samples (419 samples)
GNOMAD_GENOMES_ASJ_AF - Alternative allele frequency in the Ashkenazi Jewish gnomAD genome samples (151 samples)
GNOMAD_GENOMES_EAS_AF - Alternative allele frequency in the East Asian gnomAD genome samples (811 samples)
GNOMAD_GENOMES_FIN_AF - Alternative allele frequency in the Finnish gnomAD genome samples (1,747 samples)
GNOMAD_GENOMES_NFE_AF - Alternative allele frequency in the Non-Finnish European gnomAD genome samples (7,509 samples)
GNOMAD_GENOMES_OTH_AF - Alternative allele frequency in other gnomAD genome samples (491 samples)
CLINVAR_CLNSIG - clinical significance as to the clinvar data set. 0 - unknown, 1 - untested, 2 - Benign, 3 - Likely benign, 4 - Likely pathogenic, 5 - Pathogenic, 6 - drug response, 7 - histocompatibility. A negative score means the the score is for the ref allele
CLINVAR_TRAIT - the trait/disease the CLINVAR_CLNSIG referring to
CLINVAR_GOLDEN_STARS - ClinVar Review Status summary. 0 - no assertion criteria provided, 1 - criteria provided, single submitter, 2 - criteria provided, multiple submitters, no conflicts, 3 - reviewed by expert panel, 4 - practice guideline
GERP_RS - GERP++ RS score, the larger the score, the more conserved the site. Scores range from -12.3 to 6.17
MIN_SIFT_SCORE - Minimum SIFT score (SIFTori). Scores range from 0 to 1. The smaller the score the more likely the SNP has damaging effect
MIN_SIFT_PRED - Prediction corresponding to the minimum sift score. If SIFTori is smaller than 0.05 the corresponding nsSNV is predicted as "D(amaging)"; otherwise it is predicted as "T(olerated)"
DNDS_DISEASE_QVAL - dn/ds diseases with significant q-values (q-value < 0.05), analysed from TCGA whole-exome data in COSMIC. Diseases are classified into AML, HNSCC, NSCLC, bladder, breast, cervical, colon, endometrioid, gastric, glioma, kidney, liver, melanoma, ovary, pancreatic, prostate, sarcoma, thyroid
MUTATION_SIGNIFICANCE_TIER - Mutation significance. 1 - high significance, 2 - medium significance, 3 - low significance, Other - No predicted significance (other mutations)