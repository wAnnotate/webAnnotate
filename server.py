from flask import Flask, flash, request, jsonify, render_template, Response, redirect, send_from_directory
from perlfunc import perl5lib, perlfunc, perlreq
import vcf
import base64
from io import BufferedReader,TextIOWrapper
import bs4
import requests
app = Flask(__name__)


@perlfunc
@perlreq('modules/annotate_variation.pl')
def annotate_variation():
    pass


@perlfunc
@perlreq('modules/coding_change.pl')
def coding_change():
    pass


@perlfunc
@perlreq('modules/convert2annovar.pl')
def convert2annovar():
    pass


@perlfunc
@perlreq('modules/retrieve_seq_from_fasta.pl')
def retrieve_seq_from_fasta():
    pass


@perlfunc
@perlreq('modules/table_annovar.pl')
def table_annovar():
    pass


@perlfunc
@perlreq('modules/variants_reduction.pl')
def variants_reduction():
    pass


@app.route("/")
def index():
    return render_template("index.html")

def getGeneInfo(gene_id,table):
    data = requests.get("https://www.ncbi.nlm.nih.gov/snp/%s" % str(gene_id)).text
    soup = bs4.BeautifulSoup(data, 'html.parser')
    aS = soup.findAll('a',href=True)
    mainurl = "https://www.ncbi.nlm.nih.gov"
    endpoint = ""
    for a in aS:
        if "gene" in a["href"]:
            endpoint = a["href"]
            break
    summary = requests.get(mainurl+endpoint).text
    summarysoup = bs4.BeautifulSoup(summary,'html.parser')
    summarysoup = summarysoup.find(id = "summaryDiv")
    if not summarysoup:
        return
    dts = summarysoup.findAll('dt')
    dds = summarysoup.findAll('dd')
    for dt,dd in zip(dts,dds):
        if dt.get_text() not in table:
            table[dt.get_text()] = []
            table[dt.get_text()].append(dd.get_text())
        else:
            table[dt.get_text()].append(dd.get_text())
    local_headers = [dt.get_text() for dt in dts]
    for th in table:
        if th not in local_headers:
            table[th].append("")

@app.route("/annotate",methods=["POST"])
def annotate():
    file = request.files["efile"]
    file.name = file.filename
    file = BufferedReader(file)
    file = TextIOWrapper(file)
    print(type(file))
    vcf_reader = vcf.Reader(file)
    table = {}
    alltds = {}
    for record in vcf_reader:
        if record.ID:
            getGeneInfo(record.ID,table)
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
    return render_template("annotated.html",table=tablehtml)

if __name__ == "__main__":
    app.run()
