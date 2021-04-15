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
    query = """CREATE INDEX IF NOT EXISTS mutation_index on mutations (chr)"""
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
    conn.close()
