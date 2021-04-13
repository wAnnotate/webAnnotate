import csv

assertionPath = "static/civicdb/01-Feb-2021-AssertionSummaries.tsv"
clinicalEvidencePath = "static/civicdb/01-Feb-2021-ClinicalEvidenceSummaries.tsv"
genePath = "static/civicdb/01-Feb-2021-GeneSummaries.tsv"
variantGroupPath = "static/civicdb/01-Feb-2021-VariantGroupSummaries.tsv"
variantPath = "static/civicdb/01-Feb-2021-VariantSummaries.tsv"

civicVcfPath = "static/civicdb/01-Feb-2021-civic_accepted_and_submitted.vcf"


class CivicDb:  # GRCh37 (Ensembl v75)
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
                self.clinicalEvidences[Id] = row

    def _constructGenes(self):
        """gene_id
        gene_civic_url  name    entrez_id   description last_review_date"""
        with open(genePath, mode='r') as tsvFile:
            tsvReader = csv.DictReader(tsvFile, delimiter='\t')
            for row in tsvReader:
                Id = row["gene_id"]
                self.genes[Id] = row

    def _constructVariantGroups(self):
        """variant_group_id
        variant_group_civic_url	variant_group	description	last_review_date"""
        with open(variantGroupPath, mode='r') as tsvFile:
            tsvReader = csv.DictReader(tsvFile, delimiter='\t')
            for row in tsvReader:
                Id = row["variant_group_id"]
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
                self.variants[Id] = row

    def findVariantsFromLocation(self, chromosome, location):  # Returns array of variant dicts
        variants = []
        for v in self.variants.items():
            if "chromosome" in v[1] and "start" in v[1] and "stop" in v[1]:
                if v[1]["chromosome"] == str(chromosome) and int(v[1]["start"]) <= location <= int(v[1]["stop"]):
                    variants.append(v[1])
            elif "chromosome2" in v[1] and "start2" in v[1] and "stop2" in v[1]:
                if v[1]["chromosome2"] == str(chromosome) and int(v[1]["start2"]) <= location <= int(v[1]["stop2"]):
                    variants.append(v[1])
        return variants

    def findGene(self, arg):  # Gets either gene name or entrez id, returns gene dict
        isGeneName = False
        isEntrezId = False
        newArg = ""
        if type(arg) is str:
            try:
                arg = int(arg)
                isEntrezId = True
            except ValueError:
                isGeneName = True
        elif type(arg) is int:
            isEntrezId = True
        else:
            return Exception("Invalid arguments.")

        gene = {}
        if isGeneName:
            for g in self.genes.items():
                if g[1]["name"] == arg:
                    gene = g[1]
                    break
        elif isEntrezId:
            for g in self.genes.items():
                if g[1]["entrez_id"] == str(arg):
                    gene = g[1]
                    break
        else:
            Exception("Unknown error.")
        return gene

    def findGeneFromLocation(self, chromosome, location):  # Returns gene dict
        variants = self.findVariantsFromLocation(chromosome, location)
        gene = self.findGene(variants[0]["gene"])
        return gene

    def findVariantGroups(self, groups):  # Gets variant_groups, returns variant groups
        variantGroups = []
        groups = groups.split(',')
        for g in groups:
            for vg in self.variantGroups.items():
                if vg[1]["variant_group"] == g:
                    variantGroups.append(vg[1])
        return variantGroups

    def findAssertions(self, variantId):  # Gets variant_id, returns assertions
        assertions = []
        for a in self.assertions.items():
            if a[1]["variant_id"] == str(variantId):
                assertions.append(a[1])
        return assertions

    def findClinicalEvidences(self, variantId):  # Gets variant_id, returns clinical evidences
        clinicalEvidences = []
        for ce in self.clinicalEvidences.items():
            if ce[1]["variant_id"] == str(variantId):
                clinicalEvidences.append(ce[1])
        return clinicalEvidences


"""
db = CivicDb()
fcee = db.findClinicalEvidences(64)
for cee in fcee:
    for value in cee.items():
        print(value[0], "->", value[1])
    break
print()"""
# variant, representative_transcript, konum
