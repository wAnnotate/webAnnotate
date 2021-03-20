import csv
import sqlite3
from sqlite3 import Error

# HGNCPath = "static/cosmicdb/CosmicHGNC.tsv"
resistanceMutationsPath = "static/cosmicdb/CosmicResistanceMutations.tsv"
mutantExportCensusPath = "static/cosmicdb/CosmicMutantExportCensus.tsv"
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
        pass

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


"""
db = CosmicDb()

db.findVariantsFromLocation("GRCh37", "1", 69224)"""
"""with open("CosmicMutantExportCensus.tsv", mode='r') as tsvFile:
    tsvReader = csv.DictReader(tsvFile, delimiter='\t')
    index = 0
    for row in tsvReader:
        for val in row.items():
            print(val)
        index += 1
        break
print()
"""
