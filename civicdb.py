import csv

assertionPath = "static/civicdb/01-Feb-2021-AssertionSummaries.tsv"
clinicalEvidencePath = "static/civicdb/01-Feb-2021-ClinicalEvidenceSummaries.tsv"
genePath = "static/civicdb/01-Feb-2021-GeneSummaries.tsv"
variantGroupPath = "static/civicdb/01-Feb-2021-VariantGroupSummaries.tsv"
variantPath = "static/civicdb/01-Feb-2021-VariantSummaries.tsv"

civicVcfPath = "static/civicdb/01-Feb-2021-civic_accepted_and_submitted.vcf"


class CivicDb:
    def __init__(self):
        self.assertions = {}
        self.clinicalEvidences = {}
        self.genes = {}
        self.variantGroups = {}
        self.variants = {}

        self._constructAssertions()
        self._constructClinicalEvidences()
        self._constructGenes()
        self._constructVariantGroups()
        self._constructVariants()

    def _constructAssertions(self):
        """gene	entrez_id	variant	disease	doid	phenotypes	drugs	assertion_type	assertion_direction
        clinical_significance	acmg_codes	amp_category	nccn_guideline	nccn_guideline_version	regulatory_approval
        fda_companion_test	assertion_summary	assertion_description
        assertion_id
        evidence_item_ids	variant_id  gene_id last_review_date    assertion_civic_url evidence_items_civic_url
        variant_civic_url	gene_civic_url"""
        with open(assertionPath, mode='r') as tsvFile:
            tsvReader = csv.DictReader(tsvFile, delimiter='\t')
            for row in tsvReader:
                Id = row["assertion_id"]
                del row["assertion_id"]
                self.assertions[Id] = row

    def _constructClinicalEvidences(self):
        """gene	entrez_id	variant	disease	doid	phenotypes	drugs	drug_interaction_type	evidence_type
        evidence_direction	evidence_level	clinical_significance	evidence_statement	citation_id	source_type
        asco_abstract_id	citation	nct_ids	rating	evidence_status
        evidence_id
        variant_id	gene_id	chromosome start	stop	reference_bases	variant_bases	representative_transcript
        chromosome2	start2	stop2   representative_transcript2	ensembl_version	reference_build	variant_summary
        variant_origin	last_review_date    evidence_civic_url	variant_civic_url	gene_civic_url"""
        with open(clinicalEvidencePath, mode='r') as tsvFile:
            tsvReader = csv.DictReader(tsvFile, delimiter='\t')
            for row in tsvReader:
                Id = row["evidence_id"]
                del row["evidence_id"]
                self.clinicalEvidences[Id] = row

    def _constructGenes(self):
        """gene_id
        gene_civic_url  name    entrez_id   description last_review_date"""
        with open(genePath, mode='r') as tsvFile:
            tsvReader = csv.DictReader(tsvFile, delimiter='\t')
            for row in tsvReader:
                Id = row["gene_id"]
                del row["gene_id"]
                self.genes[Id] = row

    def _constructVariantGroups(self):
        """variant_group_id
        variant_group_civic_url	variant_group	description	last_review_date"""
        with open(variantGroupPath, mode='r') as tsvFile:
            tsvReader = csv.DictReader(tsvFile, delimiter='\t')
            for row in tsvReader:
                Id = row["variant_group_id"]
                del row["variant_group_id"]
                self.variantGroups[Id] = row

    def _constructVariants(self):
        """variant_id
        variant_civic_url   gene    entrez_id   variant summary variant_groups  chromosome  start   stop
        reference_bases variant_bases   representative_transcript   ensembl_version reference_build chromosome2 start2
        stop2   representative_transcript2  variant_types   hgvs_expressions    last_review_date
        civic_variant_evidence_score    allele_registry_id  clinvar_ids variant_aliases	assertion_ids
        assertion_civic_urls"""
        with open(variantPath, mode='r') as tsvFile:
            tsvReader = csv.DictReader(tsvFile, delimiter='\t')
            for row in tsvReader:
                Id = row["variant_id"]
                del row["variant_id"]
                self.variants[Id] = row

    def searchGenes(self, geneId):
        """Test function works well"""
        return self.genes[geneId]


"""
# test code
with open(clinicalEvidencePath, mode='r') as tsv_file:
    tsv_reader = csv.DictReader(tsv_file, delimiter='\t')
    line_count = 0
    for row in tsv_reader:
        print(row["evidence_id"], row["evidence_civic_url"])
        if line_count == 0:
            # print(f'Column names are {", ".join(row)}')
            line_count += 1
        else:
            # print(f'\t{row["gene_id"]} {row["name"]} {row["entrez_id"]}')
            line_count += 1
    print(f'Processed {line_count} lines.')
"""
