import sqlite3
from sqlite3 import Error
from datetime import datetime

cmcExportDbPath = "cosmic/sqlite/cmc_export.db"


def createConnection(dbFile):
    connection = None
    try:
        connection = sqlite3.connect(dbFile)
    except Error as err:
        print(err)
    return connection


def executeQuery(connection, query):
    try:
        cursor = connection.cursor()
        cursor.execute(query)
    except Error as err:
        print(err)


def indexDb(connection):
    query = """CREATE UNIQUE INDEX mutation_index on mutations (id)"""
    executeQuery(connection, query)
    query = """CREATE INDEX location_index 
    on locations (id, grch37_chr, grch37_start, grch37_stop, grch38_chr, grch38_start, grch38_stop)"""
    executeQuery(connection, query)


if __name__ == '__main__':
    conn = createConnection(cmcExportDbPath)
    print("Indexing rows...")
    before = datetime.now()
    indexDb(conn)
    after = datetime.now()
    print("Rows indexed.")
    print("Passed time:", after - before)

    print("Committing changes...")
    conn.commit()
    print("Finished.")
