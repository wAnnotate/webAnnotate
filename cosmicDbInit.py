import csv
import sqlite3
from sqlite3 import Error
from datetime import datetime
import cosmicDbIndex

cmcExportTsvPath = "cosmic/cmc_export.tsv"
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


def createTables(connection):
    q1 = """ CREATE TABLE IF NOT EXISTS mutations (
                id integer PRIMARY KEY,
                gene_name text,
                accession_number text,
                onc_tsg text,
                cgc_tier text,
                mutation_url text,
                legacy_mutation_id text,
                mutation_cds text,
                mutation_aa text,
                aa_mut_start integer,
                aa_mut_stop integer,
                shared_aa integer,
                genomic_wt_allele_seq text,
                genomic_mut_allele_seq text,
                aa_wt_allele_seq text,
                aa_mut_allele_seq text,
                mutation_description_cds text,
                mutation_description_aa text,
                ontology_mutation_code text,
                genomic_mutation_id text,
                mutation_genome_position_grch37 text NOT NULL,
                mutation_genome_position_grch38 text NOT NULL,
                cosmic_sample_tested integer,
                cosmic_sample_mutated integer,
                disease text,
                wgs_disease text,
                exac_af real,
                exac_afr_af real,
                exac_amr_af real,
                exac_eas_af real,
                exac_fin_af real,
                exac_nfe_af real,
                exac_sas_af real,
                gnomad_exomes_af real,
                gnomad_exomes_afr_af real,
                gnomad_exomes_amr_af real,
                gnomad_exomes_asj_af real,
                gnomad_exomes_eas_af real,
                gnomad_exomes_fin_af real,
                gnomad_exomes_nfe_af real,
                gnomad_exomes_sas_af real,
                gnomad_exomes_oth_af real,
                gnomad_genomes_af real,
                gnomad_genomes_afr_af real,
                gnomad_genomes_amr_af real,
                gnomad_genomes_asj_af real,
                gnomad_genomes_eas_af real,
                gnomad_genomes_fin_af real,
                gnomad_genomes_nfe_af real,
                gnomad_genomes_oth_af real,
                clinvar_clnsig integer,
                clinvar_trait text,
                clinvar_golden_stars integer,
                gerp_rs real,
                min_sift_score real,
                min_sift_pred text,
                dnds_disease_qval text,
                mutation_significance_tier text
            ); """
    q2 = """ CREATE TABLE IF NOT EXISTS locations (
                id integer PRIMARY KEY,
                grch37_chr text NOT NULL,
                grch37_start integer NOT NULL,
                grch37_stop integer NOT NULL,
                grch38_chr text NOT NULL,
                grch38_start integer NOT NULL,
                grch38_stop integer NOT NULL,
                FOREIGN KEY (id) REFERENCES mutations (id)
            ); """
    if connection:
        executeQuery(connection, q1)
        executeQuery(connection, q2)
    else:
        print("Error!")


def insertMutation(connection, mutationT):  # 58
    query = """INSERT INTO mutations VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,
    ?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?) """
    cursor = connection.cursor()
    cursor.execute(query, mutationT)


def insertLocation(connection, locationT):  # 7
    query = """INSERT INTO locations VALUES(?,?,?,?,?,?,?)"""
    cursor = connection.cursor()
    cursor.execute(query, locationT)


def constructDb(conn, r, i):
    mutationTuple = (i,)
    locationTuple = (i,)
    success = False
    for val in r.items():
        try:
            value = float(val[1])
        except ValueError:
            value = val[1]
        mutationTuple += (value,)
        try:
            if val[0] == "Mutation genome position GRCh37":
                values = val[1].split(':')
                chromosome = str(values[0])
                values = values[1].split('-')
                start = int(values[0])
                stop = int(values[1])
                locationTuple += (chromosome, start, stop)
            elif val[0] == "Mutation genome position GRCh38":
                values = val[1].split(':')
                chromosome = int(values[0])
                values = values[1].split('-')
                start = int(values[0])
                stop = int(values[1])
                locationTuple += (chromosome, start, stop)
        except ValueError as e:
            success = False
            # print(e, "Missing information in the line, this row will be discarded. Index:", i)
            break
        success = True
    if success:  # No missing info about coordinates
        insertMutation(conn, mutationTuple)
        insertLocation(conn, locationTuple)
    return success


if __name__ == '__main__':
    # Reset db if exists
    c = createConnection(cmcExportDbPath)
    if c:
        q = """DROP TABLE locations"""
        executeQuery(c, q)
        q = """DROP TABLE mutations"""
        executeQuery(c, q)
        print("Tables dropped.")
    else:
        print("ERROR!")

    print("Creating tables...")
    createTables(c)
    print("Tables created.")
    c.close()

    print("Inserting values...")
    before = datetime.now()
    with open(cmcExportTsvPath, mode='r') as tsvFile:
        tsvReader = csv.DictReader(tsvFile, delimiter='\t')
        c = createConnection(cmcExportDbPath)
        index = 0
        try:
            for row in tsvReader:
                if constructDb(c, row, index):
                    index += 1
        except KeyboardInterrupt:
            after = datetime.now()
            print()
            print(index, "lines added.")
            print("Passed time:", after - before)
            raise SystemExit
    print("Values inserted.")
    after = datetime.now()
    print()
    print(index, "lines added.")
    print("Passed time:", after - before)

    print("Indexing rows...")
    before = datetime.now()
    cosmicDbIndex.indexDb(c)
    after = datetime.now()
    print("Rows indexed.")
    print("Passed time:", after - before)

    print("Committing changes...")
    c.commit()
    print("Finished.")
    c.close()
