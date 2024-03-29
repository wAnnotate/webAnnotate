import csv
# import sqlite3
import requests
import base64

resistanceMutationsPath = "static/cosmicdb/CosmicResistanceMutations.tsv"
HGNCPath = "static/cosmicdb/CosmicHGNC.tsv"
cgcPath = "static/cosmicdb/cancer_gene_census.csv"
# cmcExportDbPath = "cosmic/sqlite/cmc_export.db"


"""def dict_factory(cursor, row):
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d"""


"""def createConnection(dbFile):
    connection = None
    try:
        connection = sqlite3.connect(dbFile)
    except sqlite3.Error as err:
        print(err)
    return connection"""


"""def executeQuery(connection, query, t=()):
    try:
        cursor = connection.cursor()
        connection.row_factory = sqlite3.Row
        cursor.execute(query, t)
        rows = cursor.fetchmany()
        keys = list(map(lambda x: x[0], cursor.description))
        return rows, keys
    except sqlite3.Error as err:
        print(err)
        return None, None"""


def executeQueryOnline(query):
    apikey = "1rDCrlqK9J7OHXCeeaw9MauDQ5b"
    dbowner = "burakbozdag"
    dbname = "cmc_export.db"
    sql = base64.b64encode(query.encode("ascii")).decode("ascii")
    url = "https://api.dbhub.io/v1/query"
    data = {
        "apikey": apikey,
        "dbowner": dbowner,
        "dbname": dbname,
        "sql": sql
    }
    r = requests.post(url=url, data=data)
    if not r.ok:
        return None
    decoded = r.json()
    if not decoded:  # Empty list
        return None
    rowDictList = []
    for row in decoded:
        rowDict = {}
        for column in row:
            key = column["Name"]
            value = column["Value"]
            if column["Type"] == 4:  # Integer
                value = int(value)
            elif column["Type"] == 5:  # Float
                value = float(value)
            else:  # String
                value = str(value)
            rowDict[key] = value
        rowDictList.append(rowDict)
    print("cosmic query executed")
    return rowDictList


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
    def findVariantsFromLocation(assembly, chr, pos, ref, alt):  # "GRCh37", "X", 100000, "A", "T"
        chr = str(chr)
        pos = int(pos)
        if assembly == "GRCh37":
            query = f"SELECT * FROM mutations WHERE chr = {chr} AND grch37_start <= {pos} AND grch37_stop >= {pos}"
        elif assembly == "GRCh38":
            query = f"SELECT * FROM mutations WHERE chr = {chr} AND grch38_start <= {pos} AND grch38_stop >= {pos}"
        else:
            return Exception("Unsupported assembly.")

        # conn = createConnection(dbFile=cmcExportDbPath)
        # rows, keys = executeQuery(query, t=(chr, pos, pos))
        row_dict_list = executeQueryOnline(query)
        """if rows:
            for row in rows:
                row_dict = {}
                for index in range(len(row)):
                    row_dict[keys[index]] = row[index]
                row_dict_list.append(row_dict)"""
        if row_dict_list is None:
            return []

        for rd in row_dict_list.copy():
            if rd["genomic_wt_allele_seq"] == ref and rd["genomic_mut_allele_seq"] == alt:
                return [rd]
            if rd["genomic_wt_allele_seq"] != "" or rd["genomic_mut_allele_seq"] != "":
                row_dict_list.remove(rd)

        return row_dict_list

    def findGeneFromLocation(self, chr, pos):  # GRCh37 position required, returns dict or None
        chr = str(chr)
        pos = int(pos)
        for g in self.cgc.items():
            if "Genome Location" not in g[1]:
                continue
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
            if "Entrez GeneId" not in g[1]:
                continue
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
            if "LEGACY_MUTATION_ID" in resMut[1] and resMut[1]["LEGACY_MUTATION_ID"] == str(legacyMutationId):
                resistanceMutations.append(resMut[1])
        return resistanceMutations


"""db = CosmicDb()

genes = db.findVariantsFromLocation("GRCh37", 7, 140453136, "A", "T")
for gs in genes:
    for g in gs.items():
        print(g[0], "-> ", g[1])
        print()"""
