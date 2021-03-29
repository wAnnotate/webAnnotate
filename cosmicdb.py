import csv
import sqlite3
from sqlite3 import Error

resistanceMutationsPath = "static/cosmicdb/CosmicResistanceMutations.tsv"
HGNCPath = "static/cosmicdb/CosmicHGNC.tsv"
cgcPath = "static/cosmicdb/cancer_gene_census.csv"
cmcExportDbPath = "cosmic/sqlite/cmc_export.db"


def createConnection(dbFile):
    connection = None
    try:
        connection = sqlite3.connect(dbFile)
    except Error as err:
        print(err)
    return connection


def executeQuery(connection, query, t=()):
    try:
        cursor = connection.cursor()
        cursor.execute(query, t)
        rows = cursor.fetchall()
        return rows
    except Error as err:
        print(err)
        return None


class CosmicDb:  # GRCh37 (Ensembl v75)
    def __init__(self):
        self.resistanceMutations = {}
        self.hgnc = {}
        self.cgc = {}

        self._constructResistanceMutations()
        self._constructHGNC()
        self._constructCGC()

    def _constructResistanceMutations(self):
        """Sample Name	Sample ID	Gene Name	Transcript	Census Gene	Drug Name
        MUTATION_ID	GENOMIC_MUTATION_ID	LEGACY_MUTATION_ID	AA Mutation	CDS Mutation
        Primary Tissue	Tissue Subtype 1	Tissue Subtype 2	Histology	Histology Subtype 1	Histology Subtype 2
        Pubmed Id	CGP Study	Somatic Status	Sample Type	Zygosity	Genome Coordinates (GRCh37)	Tier
        HGVSP	HGVSC	HGVSG"""
        with open(resistanceMutationsPath, mode='r') as tsvFile:
            tsvReader = csv.DictReader(tsvFile, delimiter='\t')
            for row in tsvReader:
                Id = row["Sample ID"]
                self.resistanceMutations[Id] = row

    def _constructHGNC(self):
        """COSMIC_ID	COSMIC_GENE_NAME	Entrez_id	HGNC_ID	Mutated?	Cancer_census?	Expert Curated?"""
        with open(HGNCPath, mode='r') as tsvFile:
            tsvReader = csv.DictReader(tsvFile, delimiter='\t')
            for row in tsvReader:
                Id = row["Entrez_id"]
                self.hgnc[Id] = row

    def _constructCGC(self):
        """Gene Symbol,Name,Entrez GeneId,Genome Location,Tier,Hallmark,Chr Band,Somatic,Germline,
        Tumour Types(Somatic),Tumour Types(Germline),Cancer Syndrome,Tissue Type,Molecular Genetics,Role in Cancer,
        Mutation Types,Translocation Partner,Other Germline Mut,Other Syndrome,Synonyms"""
        with open(cgcPath, mode='r') as csvFile:
            csvReader = csv.DictReader(csvFile, delimiter=',')
            for row in csvReader:
                Id = row["Entrez GeneId"]
                self.cgc[Id] = row

    @staticmethod
    def findVariantsFromLocation(assembly, chr, pos):  # "GRCh37", "20", 100000
        query = """"""
        if assembly == "GRCh37":
            query = """SELECT * FROM mutations WHERE id IN (
            (SELECT id FROM locations WHERE grch37_chr = ? AND grch37_start <= ? AND grch37_stop >= ?))"""
        elif assembly == "GRCh38":
            query = """SELECT * FROM mutations WHERE id IN (
                        (SELECT id FROM locations WHERE grch38_chr = ? AND grch38_start <= ? AND grch38_stop >= ?)
                        )"""
        else:
            return Exception("Unsupported assembly.")

        conn = createConnection(dbFile=cmcExportDbPath)
        rows = executeQuery(conn, query, (chr, pos, pos))
        return rows

    def findGeneFromLocation(self, chr, pos):  # GRCh37 position required, returns dict or None
        chr = str(chr)
        pos = int(pos)
        for g in self.cgc.items():
            chrSplit = g[1]["Genome Location"].split(':')
            if chrSplit[0] == "":  # Unknown location
                continue
            chr2 = chrSplit[0]
            posSplit = chrSplit[1].split('-')
            if posSplit[0] == "" or posSplit[1] == "":  # Unknown location
                continue
            start = int(posSplit[0])
            stop = int(posSplit[1])
            if chr == chr2 and start <= pos <= stop:
                return g[1]
        return None

    def findGene(self, entrezId):  # Gets entrez id, returns dict or None (if you cannot find with location, try this)
        for g in self.cgc.items():
            if g[1]["Entrez GeneId"] == str(entrezId):
                return g[1]
        return None

    def findHGNC(self, entrezId):  # Gets entrez id, returns dict which includes HGNC id, or None
        if str(entrezId) in self.hgnc:
            return self.hgnc[str(entrezId)]
        return None

    def findResistanceMutations(self, legacyMutationId):  # Gets LEGACY_MUTATION_ID from db, returns array of dict
        resistanceMutations = []
        for resMut in self.resistanceMutations.items():
            if resMut[1]["LEGACY_MUTATION_ID"] == str(legacyMutationId):
                resistanceMutations.append(resMut[1])
        return resistanceMutations


"""
db = CosmicDb()

gene = db.findResistanceMutations("COSM49140")
for g in gene:
    print(g)
    print()
"""
