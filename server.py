from flask import Flask, flash, request, jsonify, render_template, Response, redirect, send_from_directory
import vcf
import base64
from io import BufferedReader, TextIOWrapper
import bs4
import requests
from pyensembl import EnsemblRelease

from biothings_client import get_client

# Testing for genes
gene_client = get_client('gene')

gene_client.getgene('1017', fields='symbol,name')

gene_client.getgenes(['1017', '1018'], species='human', fields='symbol,name')

gene_client.query('uniprot:P24941', fields='symbol,name')

gene_client.querymany(['P24941', 'O14727'], scopes='uniprot', fields='symbol,name')

gene_client.metadata()
# Testing for genes

# Testing for variants
variant_client = get_client('variant')

variant_client.query('dbnsfp.genename:BTK', fields='_id')
# Testing for variants

app = Flask(__name__)


@app.route("/")
def index():
    return render_template("index.html")


def getGeneInfo(gene_id, table):
    data = gene_client.getgene(gene_id, fields='summary,clingen,entrezgene')
    table["summary"].insert(0,data["summary"])
    table["entrezgene"].insert(0,'<a href="https://www.ncbi.nlm.nih.gov/gene/%s">For details</a>' % data["entrezgene"])
    clinical_data = "no data"
    if "clingen" in data  and "clinical_validity" in data["clingen"]:
        clinical_data = ""
        for key in data["clingen"]["clinical_validity"]:
            clinical_data += ('<p>%s</p>'% data["clingen"]["clinical_validity"][key])
    table["clingen"].insert(0,clinical_data)
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

def getGeneFromLocation(dbId, chr, pos):  # Gets reference db and location of gene, returns a Gene object
    if dbId == 77 or dbId == 76 or dbId == 75:
        data = EnsemblRelease(dbId)
    else:
        return Exception("Wrong database id.")

    gene = data.genes_at_locus(contig=chr, position=pos)

    if not gene:
        return Exception("Gene not found.")

    return gene


@app.route("/annotate", methods=["POST"])
def annotate():
    file = request.files["efile"]
    db = int(request.form["db"])
    file.name = file.filename
    file = BufferedReader(file)
    file = TextIOWrapper(file)
    print(type(file))
    vcf_reader = vcf.Reader(file)
    table = {}
    table["summary"] = []
    table["clingen"] = []
    table["entrezgene"] = []
    table["db"] = [] 
    for record in vcf_reader:
        try:
            gene = getGeneFromLocation(db, record.CHROM, record.POS)
            getGeneInfo(gene[0].gene_id, table)
            table["db"].insert(0,gene.genome)
        except:
            table["summary"].append("No data available")
            table["clingen"].append("No data available")
            table["entrezgene"].append("No data available")
            table["db"].append("No data available")
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
