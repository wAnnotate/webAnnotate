from flask import Flask, flash, request, jsonify, render_template, Response, redirect, send_from_directory, session
import vcf
import base64
from io import BufferedReader, TextIOWrapper
import bs4
import requests
from pyensembl import EnsemblRelease
from biothings_client import get_client
import json
import os
import csv
from multiprocessing import Process, Manager, Pool
from flask_session import Session

from civicdb import CivicDb
import mapping  # mapping.remap()

app = Flask(__name__)
app.secret_key = b'\xdd\xd6]j\xb0\xcc\xe3mNF{\x14\xaf\xa7\xb3\x18'
dbs = (102, 75, 54)
dbName = {
    "54": "NCBI36",
    "75": "GRCh37",
    "102": "GRCh38"
}
data = EnsemblRelease(102)
SESSION_TYPE = 'filesystem'
gene_client = get_client('gene')

variant_client = get_client('variant')

manager = Manager()
app.config.from_object(__name__)
Session(app)


@app.route("/")
def index():
    return render_template("index.html")


def getGeneInfo(gene_id, table):
    geneData = gene_client.getgene(gene_id)
    if "summary" not in geneData:
        table["summary"].append("No data avaliable")
    else:
        table["summary"].append(geneData["summary"])
    if "entrezgene" not in geneData:
        table["entrezgene"].append("No data avaliable")
    else:
        table["entrezgene"].append('<a href="https://www.ncbi.nlm.nih.gov/gene/%s">%s</a>'
                                   % (geneData["entrezgene"], geneData["entrezgene"]))
    clinical_data = "no data"
    if "clingen" in geneData and "clinical_validity" in geneData["clingen"]:
        clinical_data = ""
        if type(geneData["clingen"]["clinical_validity"]) == dict:
            for key in geneData["clingen"]["clinical_validity"]:
                clinical_data += ('<p>%s</p>' % geneData["clingen"]["clinical_validity"][key])
        elif type(geneData["clingen"]["clinical_validity"]) == list:
            for key in geneData["clingen"]["clinical_validity"]:
                clinical_data += ('<p>%s, %s, %s, %s, %s</p>' % (key["classification"], key["disease_label"],
                                                                 key["mondo"], key["online_report"], key["sop"]))
    table["clingen"].append(clinical_data)
    """
    geneData = requests.get("https://www.ncbi.nlm.nih.gov/snp/%s" % str(gene_id)).text
    soup = bs4.BeautifulSoup(geneData, 'html.parser')
    aS = soup.findAll('a', href=True)
    mainurl = "https://www.ncbi.nlm.nih.gov"
    endpoint = ""
    for a in aS:
        if "gene" in a["href"]:
            endpoint = a["href"]
            break
    summary = requests.get(mainurl + endpoint).text
    summarysoup = bs4.BeautifulSoup(summary, 'html.parser')
    summarysoup = summarysoup.find(id="summaryDiv")
    if not summarysoup:
        return
    dts = summarysoup.findAll('dt')
    dds = summarysoup.findAll('dd')
    for dt, dd in zip(dts, dds):
        if dt.get_text() not in table:
            table[dt.get_text()] = []
            table[dt.get_text()].append(dd.get_text())
        else:
            table[dt.get_text()].append(dd.get_text())
    local_headers = [dt.get_text() for dt in dts]
    for th in table:
        if th not in local_headers:
            table[th].append("")
    """


@app.route("/static/images/logom.png")
def logo():
    print("logo asked for")
    return send_from_directory(os.path.join(app.root_path, 'static', 'images'),
                               'logom.PNG', mimetype='image/png')


@app.route("/constructannotation", methods=["GET"])
def constructannotationGet():
    return render_template("constructannotation.html")


@app.route("/constructannotation", methods=["POST"])
def constructannotation():
    file = request.files["efile"]
    file = BufferedReader(file)
    file = TextIOWrapper(file)
    csv_input = csv.reader(file)
    table = {
        "rowid": [],
        "expression": [],
        "gene_id": [],
        "gene_name": [],
        "biotype": [],
        "contig": [],
        "start": [],
        "end": [],
        "strand": [],
        "genome": [],
        "summary": [],
        "clingen": [],
        "entrezgene": []
    }
    count = 0
    for row in csv_input:
        if count == 0:
            count += 1
            continue
        for key, value in zip(table.keys(), row):
            if key == "expression":
                table[key].append('<a href="/annotate/%s">Expression Graph</a>' % (count - 1))
            elif key == "entrezgene":
                if "No" not in value:
                    table[key].append('<a href="https://www.ncbi.nlm.nih.gov/gene/%s">%s</a>'
                                      % (value, value))
                else:
                    table[key].append(value)
            else:
                table[key].append(value)
        count += 1

    table[" "] = []
    for i in range(len(table["rowid"])):
        table[" "].append("""
            <button onclick="toggle(this)" style="color:white;font-size:20px;" name="+" class="btn btn-success btn-lg">
            +
                </button>
            """)
    ths = [" ", "rowid", "expression", "gene_id", "gene_name", "biotype", "contig",
           "start", "end", "strand", "genome", "summary", "clingen", "entrezgene", ]
    tablehtml = """<table id = "table" class="table table-bordered"><thead><tr>"""
    for th in ths:
        tablehtml += "<th>%s</th>" % th
    tablehtml += "</tr></thead><tbody>"
    count = len(list(table.values())[0])
    for c in range(count):
        tablehtml += "<tr>"
        for th in ths:
            tablehtml += "<td>%s</td>" % table[th][c]
        tablehtml += "</tr>"
    tablehtml += "</tbody></table>"
    session["table"] = table.copy()
    return render_template("annotated.html", table=tablehtml)


@app.route("/annotate", methods=["GET"])
def prevAnnotated():
    if "table" in session:
        ths = [" ", "rowid", "expression", "gene_id", "gene_name", "biotype", "contig",
               "start", "end", "strand", "genome", "summary", "clingen", "entrezgene", ]
        table = session["table"].copy()
        tablehtml = """<table id = "table" class="table table-bordered"><thead><tr>"""
        for th in ths:
            tablehtml += "<th>%s</th>" % th
        tablehtml += "</tr></thead><tbody>"
        count = len(list(table.values())[0])
        for c in range(count):
            tablehtml += "<tr>"
            for th in ths:
                tablehtml += "<td>%s</td>" % table[th][c]
            tablehtml += "</tr>"
        tablehtml += "</tbody></table>"
        return render_template("annotated.html", table=tablehtml)
    return redirect("/")


def getGeneFromLocation(chr, pos):  # Gets location of gene, returns a Gene object (v102, v75, v54)
    global data
    gene = data.genes_at_locus(contig=chr, position=pos)
    """if not gene:
        for db in dbs:
            if db == session["dbChoice"]:
                continue
            data = EnsemblRelease(db)
            chr, pos = mapping.remap(dbName[str(session["dbChoice"])], dbName[str(db)], chr, pos)
            gene = data.genes_at_locus(contig=chr, position=pos)
            if not gene:
                continue
            break"""
    if not gene:
        return Exception("Gene not found.")
    return gene


def getGeneFromGeneId(gId):  # Gets Ensembl id, returns a Gene object
    global data
    gene = data.gene_by_id(gId)
    """if not gene:
        for db in dbs:
            if db == session["dbChoice"]:
                continue
            data = EnsemblRelease(db)
            gene = data.gene_by_id(gId)
            if not gene:
                continue
            break"""
    if not gene:
        return Exception("Gene not found.")
    return gene


def getGeneFromRsId(rsId):  # Gets rsId, returns gene object
    if session["dbChoice"] == 102:
        varData = variant_client.getvariant(rsId, assembly="hg38")
        if varData is None:
            return Exception("Gene not found from rsId.")
    elif session["dbChoice"] == 75:
        varData = variant_client.getvariant(rsId, assembly="hg19")
        if varData is None:
            return Exception("Gene not found from rsId.")
    else:
        for assembly in ("hg19", "hg38"):
            varData = variant_client.getvariant(rsId, assembly=assembly)
            if varData is None:
                continue
            break
    if varData is None:
        return Exception("Gene not found from rsId.")
    if type(varData) == dict:
        if type(varData["cadd"]["gene"]) == list:
            gene = getGeneFromGeneId(varData["cadd"]["gene"][0]["gene_id"])
        else:
            gene = getGeneFromGeneId(varData["cadd"]["gene"]["gene_id"])
    elif type(varData) == list:
        gene = getGeneFromGeneId(varData[1]["cadd"]["gene"][1]["gene_id"])
    else:
        return Exception("Unknown output format.")
    return gene


@app.route("/annotate/<rowid>", methods=["GET"])
def expression(rowid):
    rowid = int(rowid)
    # print(session["table"]["entrezgene"][rowid])
    if "ncbi" in session["table"]["entrezgene"][rowid]:
        gene = session["table"]["entrezgene"][rowid].split("href=\"")[1].split("\"")[0].split("/")[-1]
        data1 = requests.get(
            "https://www.ncbi.nlm.nih.gov/projects/Gene/download_expression.cgi?PROJECT_DESC=PRJEB4337&GENE=%s" % gene).text
        data1 = data1.split("\n\n\n")[1]
        headers = data1.split("\n")[0].split("\t")
        data1 = data1.split("\n")[1].split("\t")
        tablehtml = '<table class="table table-bordered"><thead><tr>'
        dta = []
        for header, d in zip(headers, data1):
            if header and header != "#GeneID":
                dta.append({"name": str(header), "value": float(d)})
        dta = sorted(dta, key=lambda i: i['value'])
        for header in headers:
            tablehtml += "<th>%s</th>" % header
        tablehtml += "</tr></thead><tbody><tr>"
        for d in data1:
            tablehtml += "<td>%s</td>" % d
        tablehtml += "</tr><tbody></table>"
        return render_template("visualization.html", table=tablehtml, dta=json.dumps(dta))
    return render_template("visualization.html", table="No data available", dta="[]")


def processVCFRecord(record, table, index):
    foundGene = False
    gene_dict = {}
    subdict = {
        " ": [],
        "expression": [],
        "gene_id": [],
        "gene_name": [],
        "biotype": [],
        "contig": [],
        "start": [],
        "end": [],
        "strand": [],
        "genome": [],
        "summary": [],
        "clingen": [],
        "entrezgene": []
    }
    if record.ID:  # RsId exists
        # print("rsid exists")
        try:
            gene = getGeneFromRsId(record.ID)
            gene_dict = gene.__dict__
            foundGene = True
            getGeneInfo(gene.gene_id, subdict)
        except Exception as e:
            print("getGeneFromRsId: ", e)
            foundGene = False
    if not foundGene:
        # if not record.ID: print("rsid does not exist")
        try:
            gene = getGeneFromLocation(record.CHROM, record.POS)
            gene_dict = gene[0].__dict__
            foundGene = True
            getGeneInfo(gene[0].gene_id, subdict)
        except Exception as e:
            print("getGeneFromLocation: ", e)
            foundGene = False

    if foundGene:
        print(gene)
        for key in subdict.keys():
            if key in gene_dict.keys() and key not in ["summary", "clingen", "entrezgene", "rowid", "expression",
                                                       " "]:
                subdict[key].append(str(gene_dict[key]))
            elif key not in ["summary", "clingen", "entrezgene", "rowid", "expression", " "]:
                subdict[key].append("No data available")
        subdict["expression"].append('<a href="/annotate/%s">Expression Graph</a>' % index)
        subdict[" "].append("""
                <button onclick="toggle(this)" style="color:white;font-size:20px;" name="+" class="btn btn-success btn-lg">
                +
                </button>
        """)
    else:
        for key in subdict.keys():
            if key != "rowid" and key != " ":
                subdict[key].append("No data available")
        subdict[" "].append("")
    table[index] = subdict


@app.route("/annotate", methods=["POST"])
def annotate():
    pool = Pool(os.cpu_count())
    file = request.files["efile"]
    session["dbChoice"] = int(request.form["db"])
    print("DB choice:", session["dbChoice"])
    file.name = file.filename
    file = BufferedReader(file)
    file = TextIOWrapper(file)
    # print(type(file))
    vcf_reader = vcf.Reader(file)
    ttable = manager.dict()
    table = {
        " ": [],
        "rowid": [],
        "expression": [],
        "gene_id": [],
        "gene_name": [],
        "biotype": [],
        "contig": [],
        "start": [],
        "end": [],
        "strand": [],
        "genome": [],
        "summary": [],
        "clingen": [],
        "entrezgene": []
    }
    count = 0
    processes = []
    global data
    data = EnsemblRelease(session["dbChoice"])
    print(data)
    for record in vcf_reader:
        pool.apply_async(processVCFRecord, (record, ttable, count))
        count += 1
    pool.close()
    pool.join()
    count = len(list(ttable.keys()))
    for c in range(count):
        table["rowid"].append(c)
        for item2 in ttable[c].items():
            # print(len(item2[1]))
            table[item2[0]].append(item2[1][0])
    session["table"] = table.copy()
    # print(table)
    tablehtml = """<table id = "table" class="table table-bordered"><thead><tr>"""
    for th in table:
        tablehtml += "<th>%s</th>" % th
    tablehtml += "</tr></thead><tbody>"
    count = len(list(table.values())[0])
    for c in range(count):
        tablehtml += "<tr>"
        for th in table:
            tablehtml += "<td>%s</td>" % table[th][c]
        tablehtml += "</tr>"
    tablehtml += "</tbody></table>"
    return render_template("annotated.html", table=tablehtml)


if __name__ == "__main__":
    app.run()
