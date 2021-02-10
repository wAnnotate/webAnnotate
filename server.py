from flask import Flask, flash, request, jsonify, render_template, Response, redirect, send_from_directory
import vcf
import base64
from io import BufferedReader, TextIOWrapper
import bs4
import requests
from pyensembl import EnsemblRelease
from biothings_client import get_client

data = EnsemblRelease(77)

gene_client = get_client('gene')

variant_client = get_client('variant')

app = Flask(__name__)


@app.route("/")
def index():
    return render_template("index.html")


def getGeneInfo(gene_id, table):
    data = gene_client.getgene(gene_id, fields='summary,clingen,entrezgene')
    table["summary"].insert(0, data["summary"])
    table["entrezgene"].insert(0, '<a href="https://www.ncbi.nlm.nih.gov/gene/%s">%s</a>'
                               % (data["entrezgene"], data["entrezgene"]))
    clinical_data = "no data"
    if "clingen" in data and "clinical_validity" in data["clingen"]:
        clinical_data = ""
        for key in data["clingen"]["clinical_validity"]:
            clinical_data += ('<p>%s</p>' % data["clingen"]["clinical_validity"][key])
    table["clingen"].insert(0, clinical_data)
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


def getGeneFromLocation(chr, pos):  # Gets location of gene, returns a Gene object
    global data
    gene = data.genes_at_locus(contig=chr, position=pos)
    if not gene:
        data = EnsemblRelease(76)
        gene = data.genes_at_locus(contig=chr, position=pos)
        if not gene:
            data = EnsemblRelease(75)
            gene = data.genes_at_locus(contig=chr, position=pos)
            if not gene:
                return Exception("Gene not found.")
    return gene


@app.route("/annotate", methods=["POST"])
def annotate():
    file = request.files["efile"]
    # genome = int(request.form["db"])
    file.name = file.filename
    file = BufferedReader(file)
    file = TextIOWrapper(file)
    print(type(file))
    vcf_reader = vcf.Reader(file)
    table = {
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
        "entrezgene": [],
    }
    for record in vcf_reader:
        try:
            gene = getGeneFromLocation(record.CHROM, record.POS)
            getGeneInfo(gene[0].gene_id, table)
            table["gene_id"].append(gene[0].gene_id)
            table["gene_name"].append(gene[0].gene_name)
            table["biotype"].append(gene[0].biotype)
            table["contig"].append(gene[0].contig)
            table["start"].append(gene[0].start)
            table["end"].append(gene[0].end)
            table["strand"].append(gene[0].strand)
            table["genome"].append(gene[0].genome)
            print("Successful adding.")
        except:
            table["gene_id"].append("No data available")
            table["gene_name"].append("No data available")
            table["biotype"].append("No data available")
            table["contig"].append("No data available")
            table["start"].append("No data available")
            table["end"].append("No data available")
            table["strand"].append("No data available")
            table["genome"].append("No data available")
            table["summary"].append("No data available")
            table["clingen"].append("No data available")
            table["entrezgene"].append("No data available")
            # TODO: getGeneInfo function must be adjusted in order to annotate variants with unknown RSid
    print(len(table["summary"]))
    print(len(table["clingen"]))
    print(len(table["entrezgene"]))
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
