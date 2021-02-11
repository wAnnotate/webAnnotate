from flask import Flask, flash, request, jsonify, render_template, Response, redirect, send_from_directory, session
import vcf
import base64
from io import BufferedReader, TextIOWrapper
import bs4
import requests
from pyensembl import EnsemblRelease
from biothings_client import get_client

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
    data = gene_client.getgene(gene_id, fields='summary,clingen,entrezgene')
    if "summary" not in data:
        table["summary"].append("No data avaliable")
    else:
        table["summary"].append(data["summary"])
    if "entrezgene" not in data:
        table["entrezgene"].append("No data avaliable")
    else:
        table["entrezgene"].append('<a href="https://www.ncbi.nlm.nih.gov/gene/%s">%s</a>'
                                   % (data["entrezgene"], data["entrezgene"]))
    clinical_data = "no data"
    if "clingen" in data and "clinical_validity" in data["clingen"]:
        clinical_data = ""
        for key in data["clingen"]["clinical_validity"]:
            clinical_data += ('<p>%s</p>' % data["clingen"]["clinical_validity"][key])
    table["clingen"].append(clinical_data)
    """
    data = requests.get("https://www.ncbi.nlm.nih.gov/snp/%s" % str(gene_id)).text
    soup = bs4.BeautifulSoup(data, 'html.parser')
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
def visualize(rowid):
    rowid = int(rowid)
    print(session["table"]["entrezgene"][rowid])
    if "ncbi" in session["table"]["entrezgene"][rowid]:
        data = requests.get(session["table"]["entrezgene"][rowid].split("href=\"")[1].split("\"")[0]).text
        soup = bs4.BeautifulSoup(data, 'html.parser')
        graph = soup.find(id="gene-expression-app")
        return str(graph)
    return "No data available"


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
        "rowid": [],
        "visualize": [],
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
        gene_dict = dict()
        if record.ID:  # RsId exists
            print("rsid exists")
            try:
                gene = getGeneFromRsId(record.ID)
                getGeneInfo(gene.gene_id, table)
                gene_dict = gene.__dict__
                foundGene = True
            except Exception as e:
                print(e)
                foundGene = False
        if not foundGene:
            if not record.ID:
                print("rsid does not exist")
            try:
                gene = getGeneFromLocation(record.CHROM, record.POS)
                getGeneInfo(gene[0].gene_id, table)
                gene_dict = gene[0].__dict__
                foundGene = True
            except Exception as e:
                print(e)
                foundGene = False

        if foundGene:
            for key in table.keys():
                if key in gene_dict.keys() and key not in ["summary", "clingen", "entrezgene", "rowid", "visualize"]:
                    table[key].append(str(gene_dict[key]))
                elif key not in ["summary", "clingen", "entrezgene", "rowid", "visualize"]:
                    table[key].append("No data available")
            table["rowid"].append(count)
            table["visualize"].append('<a href="/annotate/%s">Visualize</a>' % count)
        else:
            for key in table.keys():
                if key != "rowid":
                    table[key].append("No data available")
            table["rowid"].append(count)

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
                if key in gene_dict.keys() and key not in ["summary", "clingen", "entrezgene", "rowid", "visualize"]:
                    table[key].append(str(gene_dict[key]))
                elif key not in ["summary", "clingen", "entrezgene", "rowid", "visualize"]:
                    table[key].append("No data available")
            table["rowid"].append(count)
            table["visualize"].append('<a href="/annotate/%s">Visualize</a>' % count)
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
    print(table)
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
