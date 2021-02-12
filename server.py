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

dbChoice = 102
dbs = (102, 75, 54)
data = EnsemblRelease(dbChoice)

gene_client = get_client('gene')

variant_client = get_client('variant')

app = Flask(__name__)
app.secret_key = b'\xdd\xd6]j\xb0\xcc\xe3mNF{\x14\xaf\xa7\xb3\x18'


@app.route("/")
def index():
    return render_template("index.html")


def getGeneInfo(gene_id, table):
    geneData = gene_client.getgene(gene_id, fields='summary,clingen,entrezgene')
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


@app.route("/prevannotated", methods=["GET"])
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
    data = EnsemblRelease(dbChoice)
    gene = data.genes_at_locus(contig=chr, position=pos)
    if not gene:
        for db in dbs:
            if db == dbChoice:
                continue
            data = EnsemblRelease(db)
            gene = data.genes_at_locus(contig=chr, position=pos)
            if not gene:
                continue
            break
    if not gene:
        return Exception("Gene not found.")
    return gene


def getGeneFromGeneId(gId):  # Gets Ensembl id, returns a Gene object
    global data
    data = EnsemblRelease(dbChoice)
    gene = data.gene_by_id(gId)
    if not gene:
        for db in dbs:
            if db == dbChoice:
                continue
            data = EnsemblRelease(db)
            gene = data.gene_by_id(gId)
            if not gene:
                continue
            break
    if not gene:
        return Exception("Gene not found.")
    return gene


def getGeneFromRsId(rsId):  # Gets rsId, returns gene object
    if dbChoice == 102:
        varData = variant_client.getvariant(rsId, assembly="hg38")
        if varData is None:
            return Exception("Gene not found from rsId.")
    elif dbChoice == 75:
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
    print(session["table"]["entrezgene"][rowid])
    if "ncbi" in session["table"]["entrezgene"][rowid]:
        gene = session["table"]["entrezgene"][rowid].split("href=\"")[1].split("\"")[0].split("/")[-1]
        data = requests.get(
            "https://www.ncbi.nlm.nih.gov/projects/Gene/download_expression.cgi?PROJECT_DESC=PRJEB4337&GENE=%s" % gene).text
        data = data.split("\n\n\n")[1]
        headers = data.split("\n")[0].split("\t")
        data = data.split("\n")[1].split("\t")
        tablehtml = '<table class="table table-bordered"><thead><tr>'
        dta = []
        for header, d in zip(headers, data):
            if header and header != "#GeneID":
                dta.append({"name": str(header), "value": float(d)})
        dta = sorted(dta, key=lambda i: i['value'])
        for header in headers:
            tablehtml += "<th>%s</th>" % header
        tablehtml += "</tr></thead><tbody><tr>"
        for d in data:
            tablehtml += "<td>%s</td>" % d
        tablehtml += "</tr><tbody></table>"
        return render_template("visualization.html", table=tablehtml, dta=json.dumps(dta))
    return render_template("visualization.html", table="No data available", dta="[]")


@app.route("/annotate", methods=["POST"])
def annotate():
    file = request.files["efile"]
    global dbChoice
    dbChoice = int(request.form["db"])
    file.name = file.filename
    file = BufferedReader(file)
    file = TextIOWrapper(file)
    print(type(file))
    vcf_reader = vcf.Reader(file)
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
    for record in vcf_reader:
        foundGene = False
        gene_dict = {}
        if record.ID:  # RsId exists
            print("rsid exists")
            try:
                gene = getGeneFromRsId(record.ID)
                gene_dict = gene.__dict__
                foundGene = True
                getGeneInfo(gene.gene_id, table)
            except Exception as e:
                print(e)
                foundGene = False
        if not foundGene:
            if not record.ID:
                print("rsid does not exist")
            try:
                gene = getGeneFromLocation(record.CHROM, record.POS)
                gene_dict = gene[0].__dict__
                foundGene = True
                getGeneInfo(gene[0].gene_id, table)
            except Exception as e:
                print(e)
                foundGene = False

        if foundGene:
            for key in table.keys():
                if key in gene_dict.keys() and key not in ["summary", "clingen", "entrezgene", "rowid", "expression",
                                                           " "]:
                    table[key].append(str(gene_dict[key]))
                elif key not in ["summary", "clingen", "entrezgene", "rowid", "expression", " "]:
                    table[key].append("No data available")
            table["rowid"].append(count)
            table["expression"].append('<a href="/annotate/%s">Expression Graph</a>' % count)
            table[" "].append("""
                <button onclick="toggle(this)" style="color:white;font-size:20px;" name="+" class="btn btn-success btn-lg">
                +
                </button>
            """)
        else:
            for key in table.keys():
                if key != "rowid" and key != " ":
                    table[key].append("No data available")
            table["rowid"].append(count)
            table[" "].append("")

        print(count, ", ", len(table["entrezgene"]))
        count += 1
        """
        try:
            if record.ID:  # RsId exists
                print("rsid exists")
                gene = getGeneFromRsId(record.ID)
                getGeneInfo(gene.gene_id, table)
                gene_dict = gene.__dict__
                foundGene = True
            else:
                print("rsid does not exist")
                gene = getGeneFromLocation(record.CHROM, record.POS)
                getGeneInfo(gene[0].gene_id, table)
                gene_dict = gene[0].__dict__
            for key in table.keys():
                if key in gene_dict.keys() and key not in ["summary", "clingen", "entrezgene", "rowid", "expression"]:
                    table[key].append(str(gene_dict[key]))
                elif key not in ["summary", "clingen", "entrezgene", "rowid", "expression"]:
                    table[key].append("No data available")
            table["rowid"].append(count)
            table["expression"].append('<a href="/annotate/%s">Expression Graph</a>' % count)
        except Exception as exp:
            print("Exception: ", exp)
            for key in table.keys():
                if key != "rowid":
                    table[key].append("No data available")
            table["rowid"].append(count)
        print(count, ", ", len(table["entrezgene"]))
        count += 1
        """
    # print(len(table["summary"]))
    # print(len(table["clingen"]))
    # print(len(table["entrezgene"]))
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
