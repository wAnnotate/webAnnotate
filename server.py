from flask import Flask, request, render_template, Response, redirect, send_from_directory, session, \
    copy_current_request_context
import vcf
# import base64
from io import BufferedReader, TextIOWrapper
# import bs4
import requests
# from pyensembl import EnsemblRelease
from ensembl import getGenesFromLocation, getGeneFromGeneId
from biothings_client import get_client
import json
import os
import csv
from multiprocessing import Process, Manager, Pool
from flask_session import Session
from json2html import *
from civicdb import CivicDb
from cosmicdb import CosmicDb
import traceback
import time
import mapping  # mapping.remap()
import universalKeys as dictKeys  # selected dict keys and their descriptions
import cadd
from collections import OrderedDict
import threading

app = Flask(__name__)
app.secret_key = b'\xdd\xd6]j\xb0\xcc\xe3mNF{\x14\xaf\xa7\xb3\x18'
dbs = (102, 75, 54)
dbName = {
    "54": "NCBI36",
    "75": "GRCh37",
    "102": "GRCh38"
}
# data = EnsemblRelease(102)
SESSION_TYPE = 'filesystem'
gene_client = get_client('gene')
variant_client = get_client('variant')
civic = CivicDb()
cosmic = CosmicDb()
manager = Manager()
app.config.from_object(__name__)
Session(app)
biothingsAssembly = "hg19"
tempSession = {}

ensembl_keys = {
    "strand":"Strand",
    "id": "Gene Id",
    "start": "Gene Start",
    "end": "Gene End",
    "description": "Gene Description",
    "biotype": "Gene Type",
    "external_name": "Gene Name"
}


@app.route("/")
def index():
    if "stamp" in session:
        print(time.time() - session["stamp"])

    if "done" in session or ("stamp" in session and str(session["stamp"]) + "done" in tempSession):
        if "done" in session:
            del session["done"]
        if "stamp" in session:
            if session["stamp"] in tempSession:
                del tempSession[session["stamp"]]
            if str(session["stamp"]) + "done" in tempSession:
                del tempSession[str(session["stamp"]) + "done"]
            del session["stamp"]
    if "stamp" in session and time.time() - session["stamp"] > 1000:
        if session["stamp"] in tempSession:
            del tempSession[session["stamp"]]
        if str(session["stamp"]) + "done" in tempSession:
            del tempSession[str(session["stamp"]) + "done"]
        del session["stamp"]
    return render_template("index.html")


def getGeneInfo(gene_id, table):
    geneData = gene_client.getgene(gene_id)
    if type(geneData) == dict:
        print(geneData.keys())
    if "summary" not in geneData:
        table["Summary"] = ("")
    else:
        table["Summary"] = (geneData["summary"])
    if "entrezgene" not in geneData:
        table["Entrezgene"] = ("")
    else:
        table["Entrezgene"] = ('<a target=\"_blank\" href="https://www.ncbi.nlm.nih.gov/gene/%s">%s</a>'
                               % (geneData["entrezgene"], geneData["entrezgene"]))
    clinical_data = ""
    if "clingen" in geneData and "clinical_validity" in geneData["clingen"]:
        clinical_data = ""
        if type(geneData["clingen"]["clinical_validity"]) == dict:
            for key in geneData["clingen"]["clinical_validity"]:
                clinical_data += ('<p>%s</p>' % geneData["clingen"]["clinical_validity"][key])
        elif type(geneData["clingen"]["clinical_validity"]) == list:
            for key in geneData["clingen"]["clinical_validity"]:
                clinical_data += ('<p>%s, %s, %s, %s, %s</p>' % (key["classification"], key["disease_label"],
                                                                 key["mondo"], key["online_report"], key["sop"]))
    table["Clingen"] = (clinical_data)


@app.route("/static/images/logom.png")
def logo():
    print("logo asked for")
    return send_from_directory(os.path.join(app.root_path, 'static', 'images'),
                               'logom.png', mimetype='image/png')


@app.route("/constructannotation", methods=["GET"])
def constructannotationGet():
    return render_template("constructannotation.html")


@app.route("/constructannotation", methods=["POST"])
def constructannotation():
    file = request.files["efile"]
    #file = BufferedReader(file)
    #file = TextIOWrapper(file)
    nnewtable = json.load(file)
    temp = {}
    keys = ["General","Civic","Cosmic","Cadd"]
    for key in nnewtable:
        temp[key] = {}
        for key2 in keys:
            temp[key][key2] = []
            for element in nnewtable[key][key2]:
                temp[key][key2].append(element)
    nnewtable = temp.copy()        
    print(nnewtable)
    session["table"] = nnewtable.copy()
    mainKeys = '''
                <label for="Cosmic">
                    <input type="checkbox" id="Cosmic" onclick="changeSelectText(this.parentElement)" onchange="resetSubkeys(this.parentElement,this)" />
                   Cosmic
                </label>
                <label for="Civic">
                    <input type="checkbox" id="Civic" onclick="changeSelectText(this.parentElement)" onchange="resetSubkeys(this.parentElement,this)" />
                   Civic
                </label>
                <label for="Cadd">
                    <input type="checkbox" id="Cadd" onclick="changeSelectText(this.parentElement)" onchange="resetSubkeys(this.parentElement,this)" />
                   Cadd
                </label>
                <label for="General">
                    <input type="checkbox" id="Cadd" onclick="changeSelectText(this.parentElement)" onchange="resetSubkeys(this.parentElement,this)" />
                   General
                </label>
                '''
    newtablehtml = ""
    newtablehtmlheader = """<thead><tr>"""
    newtablehtmlheader += "<th></th>"
    newtablehtmlheader += "<th>Row Index</th>"
    keys = {}
    for key in list(list(nnewtable.values())[0].keys()):
        keys[key] = {}
    newtablehtmlbody = "<tbody>"
    c = 0
    keyc = 2
    rowc = 0
    addedkeys = []
    popupdata = {}
    for item1 in list(nnewtable.items()):
        lenn = len(list(item1[1].items()))
        print(lenn)
        for subitem in zip(*item1[1].values()):
            newtablehtmlbody += """<tr><td>
                <br>
                <button onclick="toggle(this)" style="color:white;font-size:20px;" name="+" class="btn btn-success btn-lg">
                +
                </button></td><td>%s</td>""" % item1[0]
            for i in range(lenn):
                mainKey = list(keys.keys())[i]
                for item in list(subitem[i].items()):
                    val = item[1]
                    if type(val) == list and (type(val[0]) is OrderedDict or type(val[0]) is dict):
                        newtablehtmlbody += "<td>%s</td>" % (
                            getInnerAndHeaderHtmls(val,
                                                   "%s-%s-%s" % (item[0], rowc, c),
                                                   popupdata))
                        newtablehtmlheader += "<th>%s</th>" % item[0] if item[0] not in addedkeys else ""
                        if item[0] not in addedkeys:
                            keyc = addHeaderKeys(addedkeys, keys[mainKey], item[0], item[0], keyc)
                    elif type(val) == dict:
                        for key in val:
                            newtablehtmlbody += "<td>%s</td>" % (val[key])
                            newtablehtmlheader += "<th>%s</th>" % key if key not in addedkeys else ""
                            if key not in addedkeys:
                                keyc = addHeaderKeys(addedkeys, keys[mainKey], item[0], key, keyc, True)
                    else:
                        newtablehtmlbody += "<td>%s</td>" % (val)
                        newtablehtmlheader += "<th>%s</th>" % item[0] if item[0] not in addedkeys else ""
                        if item[0] not in addedkeys:
                            keyc = addHeaderKeys(addedkeys, keys[mainKey], item[0], item[0], keyc)
                    c += 1
            newtablehtmlbody += "</tr>"
        rowc += 1
    c = 1
    newtablehtmlheader += "</tr></thead>"
    newtablehtmlbody += "</tbody>"
    newtablehtml = """<table id = "table" class="table table-bordered">%s%s</table>""" % (
    newtablehtmlheader, newtablehtmlbody)
    return render_template("annotated.html", table=newtablehtml, mainKeys=mainKeys, subkeys=json.dumps(keys),
                               allData=nnewtable, popupdata=popupdata)


@app.route("/annotate", methods=["GET"])
def prevAnnotated():
    if "table" in session:
        nnewtable = session["table"].copy()
        mainKeys = '''
                <label for="Cosmic">
                    <input type="checkbox" id="Cosmic" onclick="changeSelectText(this.parentElement)" onchange="resetSubkeys(this.parentElement,this)" />
                   Cosmic
                </label>
                <label for="Civic">
                    <input type="checkbox" id="Civic" onclick="changeSelectText(this.parentElement)" onchange="resetSubkeys(this.parentElement,this)" />
                   Civic
                </label>
                <label for="Cadd">
                    <input type="checkbox" id="Cadd" onclick="changeSelectText(this.parentElement)" onchange="resetSubkeys(this.parentElement,this)" />
                   Cadd
                </label>
                <label for="General">
                    <input type="checkbox" id="Cadd" onclick="changeSelectText(this.parentElement)" onchange="resetSubkeys(this.parentElement,this)" />
                   General
                </label>
                '''
        newtablehtml = ""
        newtablehtmlheader = """<thead><tr>"""
        newtablehtmlheader += "<th></th>"
        newtablehtmlheader += "<th>Row Index</th>"
        keys = {}
        for key in list(list(nnewtable.values())[0].keys()):
            keys[key] = {}
        newtablehtmlbody = "<tbody>"
        c = 0
        keyc = 2
        rowc = 0
        addedkeys = []
        popupdata = {}
        for item1 in list(nnewtable.items()):
            lenn = len(list(item1[1].items()))
            print(lenn)
            for subitem in zip(*item1[1].values()):
                newtablehtmlbody += """<tr><td>
                    <br>
                    <button onclick="toggle(this)" style="color:white;font-size:20px;" name="+" class="btn btn-success btn-lg">
                    +
                    </button></td><td>%s</td>""" % item1[0]
                for i in range(lenn):
                    mainKey = list(keys.keys())[i]
                    for item in list(subitem[i].items()):
                        val = item[1]
                        if type(val) == list and (type(val[0]) is OrderedDict or type(val[0]) is dict):
                            newtablehtmlbody += "<td>%s</td>" % (
                                getInnerAndHeaderHtmls(val,
                                                       "%s-%s-%s" % (item[0], rowc, c),
                                                       popupdata))
                            newtablehtmlheader += "<th>%s</th>" % item[0] if item[0] not in addedkeys else ""
                            if item[0] not in addedkeys:
                                keyc = addHeaderKeys(addedkeys, keys[mainKey], item[0], item[0], keyc)
                        elif type(val) == dict:
                            for key in val:
                                newtablehtmlbody += "<td>%s</td>" % (val[key])
                                newtablehtmlheader += "<th>%s</th>" % key if key not in addedkeys else ""
                                if key not in addedkeys:
                                    keyc = addHeaderKeys(addedkeys, keys[mainKey], item[0], key, keyc, True)
                        else:
                            newtablehtmlbody += "<td>%s</td>" % (val)
                            newtablehtmlheader += "<th>%s</th>" % item[0] if item[0] not in addedkeys else ""
                            if item[0] not in addedkeys:
                                keyc = addHeaderKeys(addedkeys, keys[mainKey], item[0], item[0], keyc)

                        c += 1
                newtablehtmlbody += "</tr>"
            rowc += 1

        c = 1
        newtablehtmlheader += "</tr></thead>"
        newtablehtmlbody += "</tbody>"
        newtablehtml = """<table id = "table" class="table table-bordered">%s%s</table>""" % (
        newtablehtmlheader, newtablehtmlbody)
        return render_template("annotated.html", table=newtablehtml, mainKeys=mainKeys, subkeys=json.dumps(keys),
                               allData=nnewtable, popupdata=popupdata)
    return redirect("/")


"""def getGenesFromLocation(chr, pos):  # Gets location of gene, returns a Gene object (v102, v75, v54)
    global data
    gene = data.genes_at_locus(contig=chr, position=pos)
    if not gene:
        for db in dbs:
            if db == session["dbChoice"]:
                continue
            data = EnsemblRelease(db)
            chr, pos = mapping.remap(dbName[str(session["dbChoice"])], dbName[str(db)], chr, pos)
            gene = data.genes_at_locus(contig=chr, position=pos)
            if not gene:
                continue
            break
    if not gene:
        return Exception("Gene not found.")
    return gene"""


"""def getGeneFromGeneId(gId):  # Gets Ensembl id, returns a Gene object
    global data
    gene = data.gene_by_id(gId)
    if not gene:
        for db in dbs:
            if db == session["dbChoice"]:
                continue
            data = EnsemblRelease(db)
            gene = data.gene_by_id(gId)
            if not gene:
                continue
            break
    if not gene:
        return Exception("Gene not found.")
    return gene"""


def getGeneFromRsId(rsId):  # Gets rsId, returns gene object
    if session["dbChoice"] == 102:
        varData = variant_client.getvariant(rsId, assembly="hg38")
        if varData is None:
            return Exception("Gene not found from rsId."), varData
    elif session["dbChoice"] == 75:
        varData = variant_client.getvariant(rsId, assembly="hg19")
        if varData is None:
            return Exception("Gene not found from rsId."), varData
    else:
        for assembly in ("hg19", "hg38"):
            varData = variant_client.getvariant(rsId, assembly=assembly)
            if varData is None:
                continue
            break
    if varData is None:
        return Exception("Gene not found from rsId."), varData
    if type(varData) == dict:
        if type(varData["cadd"]["gene"]) == list:
            gene = getGeneFromGeneId(varData["cadd"]["gene"][0]["gene_id"])
        else:
            gene = getGeneFromGeneId(varData["cadd"]["gene"]["gene_id"])
    elif type(varData) == list:
        gene = getGeneFromGeneId(varData[1]["cadd"]["gene"][1]["gene_id"])
    else:
        return Exception("Unknown output format."), None
    return gene, varData


@app.route("/getVariantData", methods=["GET"])
def getVariantsData():
    gene = int(request.args.get("gene"))
    variant = int(request.args.get("variant"))
    variantType = str(request.args.get("type"))
    if variantType != "civic" and variantType != "cosmic":
        return Response(response=session["table"]["listofvariants"][gene][variant])
    elif variantType == "civic":
        return Response(response=session["table"]["listofvariantscivic"][gene][variant])
    return Response(response=session["table"]["listofvariantscosmic"][gene][variant])


@app.route("/annotate/<rowid>", methods=["GET"])
def expression(rowid):
    rowid = str(rowid)

    # print(session["table"]["entrezgene"][rowid])
    if "ncbi" in session["table"][rowid]["General"][0]["entrezgene"]:
        gene = session["table"][rowid]["General"][0]["entrezgene"].split("href=\"")[1].split("\"")[0].split("/")[-1]
        if requests.get(
                "https://www.ncbi.nlm.nih.gov/projects/Gene/download_expression.cgi?PROJECT_DESC=PRJEB4337&GENE=%s" % gene).status_code == 500:
            return render_template("visualization.html", table="", dta="[]")
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
    return render_template("visualization.html", table="", dta="[]")


def processVariantData(variant, count, hgsvs, index):
    variantdata = None
    if str(variant["_id"]) in hgsvs:
        variantdata = '<option value="%s-%s-%s">%s</option>' % (index, count, "biothings", variant["_id"])
        observed = "observed" if variant["observed"] else "not observed"
        varDict = {
            "header": "Biothings: " + str(variant["_id"] + " (" + observed + ")")
        }
        for key in variant:
            if variant[key] and key not in ["_id", "_version", "chrom", "hg19", "observed"]:
                if type(variant[key]) != dict or all(isinstance(x, float) for x in list(variant[key].values())):
                    if len(list(variant[key].values())) > 2 and all(
                            isinstance(x, float) for x in list(variant[key].values())):
                        variant[key]["graph"] = """
                                                    <button onclick="showGraph('%s',this)">Show Graph</button>
                                                """ % (key)
                    varDict[key] = json2html.convert(json=variant[key])
                elif type(variant[key]) == dict:
                    varDict[key] = {}
                    for key2 in variant[key]:
                        if type(variant[key][key2]) == dict and all(
                                isinstance(x, float) for x in list(variant[key][key2].values())):
                            variant[key][key2]["graph"] = """
                                                    <button onclick="showGraph('%s---%s',this)">Show Graph</button>
                                                """ % (key2, key)
                        varDict[key][key2] = json2html.convert(json=variant[key][key2], escape=False)
                        varDict[key][key2] += "<div id='chart%s---%s'></div>" % (key2, key)
        return json.dumps(varDict), variantdata
    return None, None


def getVariantData(id, assembly, index, count, hgsvs, variantdata):
    variant = variant_client.getvariant(id, assembly=assembly)
    if type(variant) == dict and "_id" in variant:
        vardict, vardata = processVariantData(variant, count, hgsvs, index)
        if vardict:
            variantdata += vardata
            count += 1
    elif type(variant) == list:
        for var in variant:
            if "_id" not in var:
                continue
            print(str(var["_id"]))
            vardict, vardata = processVariantData(var, count, hgsvs, index)
            if vardict:
                variantdata += vardata
                count += 1
    return count, variantdata


def formatListOfLists(listofdicts):
    for dictionary in listofdicts:
        keys = list(dictionary.keys())
        for key in keys:
            if not dictionary[key]:
                del dictionary[key]
    return listofdicts


def formatDataForKeys(data, mainkeys, dictDescFunc):
    for cs in data:
        keys = list(cs.keys())
        for key in keys:
            if key not in mainkeys:
                del cs[key]
            else:
                data = cs[key]
                cs[dictDescFunc(key)] = data
                del cs[key]


def processVCFRecord(record, index, nnewtable, value):
    mappedChr, mappedPos = mapping.remap(dbName[str(session["dbChoice"])], "GRCh37", record.CHROM, record.POS)
    foundGene = False
    gene_dict = {}
    main_sub_dict = {
        index: {
            "General": [],
            "Civic": [],
            "Cosmic": [],
            "Cadd": [],
        }
    }
    try:
        civic_variants_template = {}
        for key in dictKeys.civicVariants:
            civic_variants_template[dictKeys.civicDesc(key)] = ""
        cosmic_CMC_template = {}
        for key in dictKeys.cosmicCMC:
            cosmic_CMC_template[dictKeys.cosmicDesc(key)] = ""
    except:
        traceback.print_exc()
        return
    if record.ALT:
        for i in record.ALT:
            if len(str(i)) == 1 and len(str(record.REF)) == 1:
                main_sub_dict[index]["General"].append({
                    "Chromosome":str(record.CHROM),
                    "Position":str(record.POS),
                    "Reference Bases":str(record.REF),
                    "Alternative Bases":str(i),
                    "Summary": "",
                    "Entrezgene": "",
                    "Strand": "",
                    "Gene Start": "",
                    "Gene End": "",
                    "Gene Type": "",
                    "Gene Name": "",
                    "Expression": "",
                    "Gene Description":"",
                    "Clingen": "",
                    "Gene Id":"",
            })
    if record.ID and "rs" in record.ID:  # RsId exists
        # print("rsid exists")
        try:
            gene, variant = getGeneFromRsId(record.ID)
            gene_dict = gene.__dict__
            foundGene = True
            cnt = 0
            if record.ALT:
                for i in record.ALT:
                    if len(str(i)) == 1 and len(str(record.REF)) == 1:
                        getGeneInfo(gene_dict["id"], main_sub_dict[index]["General"][cnt])
                        for key in gene_dict.keys():
                            if key in list(main_sub_dict[index]["General"][cnt].keys()) or key in ensembl_keys:
                                if key in ensembl_keys:
                                    if key == "id":
                                        main_sub_dict[index]["General"][cnt][ensembl_keys[key]] = "<a target=\"_blank\" href=\"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=%s\">%s</a>" % (str(gene_dict[key]),str(gene_dict[key]))
                                    else:
                                        main_sub_dict[index]["General"][cnt][ensembl_keys[key]] = str(gene_dict[key])
                                else:
                                    main_sub_dict[index]["General"][cnt][key.replace("_"," ")] = str(gene_dict[key])
                        main_sub_dict[index]["General"][cnt]["Expression"] = '<a target=\"_blank\" href="/annotate/%s">Expression Graph</a>' % (index)
                        cnt += 1
        except Exception as e:
            print("getGeneFromRsId: ", e)
            foundGene = False

    if not foundGene:
        # if not record.ID: print("rsid does not exist")
        try:
            gene = getGenesFromLocation(record.CHROM, record.POS, session["dbChoice"])
            gene_dict = gene[0].__dict__
            foundGene = True
            cnt = 0
            if record.ALT:
                for i in record.ALT:
                    if len(str(i)) == 1 and len(str(record.REF)) == 1:
                        getGeneInfo(gene_dict["id"], main_sub_dict[index]["General"][cnt])
                        for key in gene_dict.keys():
                            if key in list(main_sub_dict[index]["General"][cnt].keys()) or key in ensembl_keys:
                                if key in ensembl_keys:
                                    if key == "id":
                                        main_sub_dict[index]["General"][cnt][ensembl_keys[key]] = "<a target=\"_blank\" href=\"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=%s\">%s</a>" % (str(gene_dict[key]),str(gene_dict[key]))
                                    else:
                                        main_sub_dict[index]["General"][cnt][ensembl_keys[key]] = str(gene_dict[key])
                                else:
                                    main_sub_dict[index]["General"][cnt][key.replace("_"," ")] = str(gene_dict[key])
                        main_sub_dict[index]["General"][cnt]["Expression"] = '<a target=\"_blank\" href="/annotate/%s">Expression Graph</a>' % (index)
                        cnt += 1
        except Exception as e:
            print("getGenesFromLocation: ", e)
            foundGene = False
    try:
        count = 0
        mappedChr, mappedPos = mapping.remap(dbName[str(session["dbChoice"])], "GRCh37", record.CHROM, record.POS)
        variantdata = ""
        hgsvs = []
    except:
        print(traceback.format_exc())
    try:
        if record.ALT:
            for i in record.ALT:
                hgsvs.append(str(variant_client.format_hgvs(mappedChr, mappedPos, record.REF, str(i))))
    except:
        print(traceback.format_exc())
    try:
        cnt = 0
        if record.ALT:
            for i in record.ALT:
                if len(str(i)) == 1 and len(str(record.REF)) == 1:
                    cadd_data = cadd.getSNV("GRCh37", mappedChr, mappedPos, record.REF, str(i))
                    if cadd_data and type(cadd_data) == list:
                        d = cadd_data[0]
                        cadd_dict = {}
                        for key in cadd.keys:
                            if d[key] == "NA":
                                d[key] = ""
                            if key in d:
                                cadd_dict[key.replace("-", " ").replace("_", " ")] = d[key]
                            else:
                                cadd_dict[key] = ""
                        main_sub_dict[index]["Cadd"].append(cadd_dict)
                    else:
                        temp1 = {}
                        for key2 in cadd.keys:
                            temp1[key2] = ""
                        main_sub_dict[index]["Cadd"].append(temp1)
                    main_sub_dict[index]["General"][cnt]["Alternative Bases"] = str(i)
                    cnt += 1
    except:
        print(traceback.format_exc())
    try:
        if record.ID:
            count, variantdata = getVariantData(record.ID, biothingsAssembly, index, count, hgsvs,
                                                variantdata)
        elif hgsvs:
            for hgsv in hgsvs:
                count, variantdata = getVariantData(hgsv, biothingsAssembly, index, count, hgsvs,
                                                    variantdata)
    except:
        print(traceback.format_exc())
    try:

        cnt = 0
        if record.ALT:
            for i in record.ALT:
                if len(str(i)) != 1 or len(str(record.REF)) != 1:
                    continue
                civicdata = civic.findVariantsFromLocation(mappedChr, mappedPos, str(record.REF), str(i)).copy()
                if civicdata:
                    count = 0
                    for var in civicdata:
                        variant_groups = formatListOfLists(civic.findVariantGroups(var["variant_groups"])).copy()
                        assertions = formatListOfLists(civic.findAssertions(var["variant_id"])).copy()
                        clinical_significances = formatListOfLists(
                            civic.findClinicalEvidences(var["variant_id"])).copy()
                        civicdatagene = civic.findGeneFromLocation(mappedChr, mappedPos).copy()
                        civicdatagene = formatListOfLists([civicdatagene])
                        formatDataForKeys(clinical_significances, dictKeys.civicClinicalEvidences, dictKeys.civicDesc)
                        formatDataForKeys(assertions, dictKeys.civicAssertions, dictKeys.civicDesc)
                        formatDataForKeys(variant_groups, dictKeys.civicVariantGroups, dictKeys.civicDesc)
                        formatDataForKeys(civicdatagene, dictKeys.civicGenes, dictKeys.civicDesc)
                        genehtml = civicdatagene
                        keys = list(var.keys())
                        for key in keys:
                            if key in ["variant"]:
                                continue
                            if key not in dictKeys.civicVariants:
                                del var[key]
                                continue

                            data = var[key]
                            var[dictKeys.civicDesc(key)] = data
                            del var[key]
                        temp = {}
                        for key in list(civic_variants_template.keys()):
                            temp[key] = var[key]
                        var = temp.copy()
                        html = var
                        CivicDict = {
                            "Variants": html if html else civic_variants_template,
                            "Variant Groups": variant_groups if variant_groups else "",
                            "Genes": genehtml if genehtml else "",
                            "Assertions": assertions if assertions else "",
                            "Clinical Evidences": clinical_significances if clinical_significances else "",
                        }
                        main_sub_dict[index]["Civic"].append(CivicDict)
                        count += 1
                else:
                    main_sub_dict[index]["Civic"].append(
                        {
                            "Variants": civic_variants_template,
                            "Variant Groups": "",
                            "Genes": "",
                            "Assertions": "",
                            "Clinical Evidences": "",
                        }
                    )
                main_sub_dict[index]["General"][cnt]["Alternative Bases"] = str(i)
                cnt += 1
    except:
        print(traceback.format_exc())
    try:
        cnt = 0
        if record.ALT:
            for i in record.ALT:
                if len(str(i)) != 1 or len(str(record.REF)) != 1:
                    continue
                cosmicdata = cosmic.findVariantsFromLocation("GRCh37", mappedChr, mappedPos, str(record.REF), str(i))
                if not cosmicdata:
                    mappedChr, mappedPos = mapping.remap(dbName[str(session["dbChoice"])], "GRCh38", record.CHROM,
                                                         record.POS)
                    cosmicdata = cosmic.findVariantsFromLocation("GRCh38", mappedChr, mappedPos)
                if cosmicdata:
                    count = 0
                    for row in cosmicdata:
                        if "legacy_mutation_id" in row and row["legacy_mutation_id"]:
                            resistanceMutations = cosmic.findResistanceMutations(row["legacy_mutation_id"])
                            for res in resistanceMutations:
                                reskeys = list(res.keys())
                                for key in reskeys:
                                    if key not in dictKeys.cosmicResistanceMutations:
                                        del res[key]
                                    else:
                                        res[dictKeys.cosmicDesc(key)] = res[key]
                                        del res[key]
                            row["legacy_mutation_id"] = "<a target=\"_blank\" href=\"https://cancer.sanger.ac.uk/cosmic/search?q=%s\">%s</a>" % (row["legacy_mutation_id"],row["legacy_mutation_id"])
                        keys = list(row.keys())
                        for key in keys:
                            if (key not in dictKeys.cosmicCMC):
                                del row[key]
                            elif not row[key]:
                                row[key] = ""
                            else:
                                row[dictKeys.cosmicDesc(key)] = row[key]
                                del row[key]
                        temp = {}
                        for key in list(cosmic_CMC_template.keys()):
                            if key in row:
                                temp[key] = row[key]
                            else:
                                temp[key] = ""
                        row = temp
                        variantData = json2html.convert(json=row, escape=False)
                        variantData = row
                        CosmicDict = {
                            "CMC": variantData if variantData else cosmic_CMC_template,
                            "Resistance Mutations": resistanceMutations if resistanceMutations else ""
                        }
                        main_sub_dict[index]["Cosmic"].append(CosmicDict)
                        count += 1
                else:
                    main_sub_dict[index]["Cosmic"].append(
                        {
                            "CMC": cosmic_CMC_template,
                            "Resistance Mutations": ""
                        }
                    )
                main_sub_dict[index]["General"][cnt]["Alternative Bases"] = str(i)
                cnt += 1
    except Exception as exp:
        print(index, "- variant exp: ", exp)
        print(traceback.format_exc())
    try:
        greatest_len = 0
        for key in ["Cosmic", "Civic", "Cadd", "General"]:
            if len(main_sub_dict[index][key]) > greatest_len:
                greatest_len = len(main_sub_dict[index][key])
        for key in ["Cosmic", "Civic", "Cadd", "General"]:
            while greatest_len > len(main_sub_dict[index][key]):
                if key == "Cosmic":
                    main_sub_dict[index][key].append(
                        {
                            "CMC": cosmic_CMC_template,
                            "Resistance Mutations": ""
                        }
                    )
                elif key == "Cadd":
                    temp1 = {}
                    for key2 in cadd.keys:
                        temp1[key2] = ""
                    main_sub_dict[index][key].append(temp1)
                elif key == "Civic":
                    main_sub_dict[index][key].append(
                        {
                            "Variants": civic_variants_template,
                            "Variant Groups": "",
                            "Genes": "",
                            "Assertions": "",
                            "Clinical Evidences": "",
                        }
                    )
                else:
                    main_sub_dict[index][key].append({
                    "Chromosome":str(record.CHROM),
                    "Position":str(record.POS),
                    "Reference Bases":str(record.REF),
                    "Alternative Bases":str(record.ALT),
                    "Summary": "",
                    "Entrezgene": "",
                    "Strand": "",
                    "Gene Start": "",
                    "Gene End": "",
                    "Gene Type": "",
                    "Gene Name": "",
                    "Expression": "",
                    "Gene Description":"",
                    "Clingen": "",
                    "Gene Id":"",
            })
    except:
        print(traceback.format_exc())

    if greatest_len > 0:
        nnewtable[str(index)] = {
            "General": main_sub_dict[index]["General"],
            "Civic": main_sub_dict[index]["Civic"],
            "Cosmic": main_sub_dict[index]["Cosmic"],
            "Cadd": main_sub_dict[index]["Cadd"],
        }
    value.value = value.value + 1


@app.route("/showresult", methods=["GET"])
def showresult():
    global tempSession
    if "stamp" in session and session["stamp"] in tempSession:
        try:
            (newtablehtml, mainKeys, keys, nnewtable, popupdata) = tempSession[session["stamp"]]
            session["result"] = tempSession[session["stamp"]]
            session["table"] = nnewtable.copy()
            return render_template("annotated.html", table=newtablehtml, mainKeys=mainKeys, subkeys=json.dumps(keys),
                                   allData=nnewtable, popupdata=popupdata)
        except:
            if "result" in session:
                (newtablehtml, mainKeys, keys, nnewtable, popupdata) = session["result"]
                session["table"] = nnewtable.copy()
                return render_template("annotated.html", table=newtablehtml, mainKeys=mainKeys,
                                       subkeys=json.dumps(keys), allData=nnewtable, popupdata=popupdata)
            else:
                print("exception")
                return redirect("/annotate")


    elif "result" in session:
        (newtablehtml, mainKeys, keys, nnewtable, popupdata) = session["result"]
        session["table"] = nnewtable.copy()
        return render_template("annotated.html", table=newtablehtml, mainKeys=mainKeys, subkeys=json.dumps(keys),
                               allData=nnewtable, popupdata=popupdata)
    print("lastresult")
    return redirect("/annotate")


@app.route("/isresult", methods=['GET'])
def isresult():
    print(session.keys())
    progress = 0
    try:
        if "stamp" in session and session["stamp"] in tempSession:
            session["result"] = tempSession[session["stamp"]]
            session["done"] = "done"
            return Response(response=json.dumps({"status": 1}))
    except:
        if "table" in session:
            return Response(response=json.dumps({"status": 1}))
    if "stamp" in session:
        if str(session["stamp"]) + "progress" in tempSession and str(session["stamp"]) + "count" in tempSession:
            progress = tempSession[str(session["stamp"]) + "progress"].value / tempSession[
                str(session["stamp"]) + "count"] * 100
    return Response(response=json.dumps({"status": 0, "progress": progress}))


def addHeaderKeys(addedkeys, civickeys, itemkey, key, keycount, isdict=False):
    addedkeys.append(key)
    if isdict:
        if itemkey not in civickeys:
            civickeys[itemkey] = {}
        civickeys[itemkey][key] = keycount
    else:
        civickeys[itemkey] = keycount
    return keycount + 1


def getInnerAndHeaderHtmls(elements, key, popupdata):
    innerhtml = "<option value="" selected>""</option>"
    for element in elements:
        header = list(element.values())[0]
        if "CIViC: Disease" in element:
            header = element["CIViC: Disease"]
        innerhtml += """
                        <option value="%s-%s">%s</option>
                     """ % (key, header, header)
        popupdata[key + "-" + header] = json2html.convert(json=element)
    innerhtmlselect = """<select onchange="toggleModal(this)">%s</select>
                      """ % (innerhtml)
    return innerhtmlselect


@app.route("/haveaprocess", methods=["GET"])
def haveaprocess():
    if "stamp" in session:
        return Response(response=json.dumps({"status": 1}))
    return Response(response=json.dumps({"status": 0}))


@app.route("/annotate", methods=["POST"])
def annotate():
    try:
        global tempSession
        stamp = time.time()
        if "stamp" in session:
            if session["stamp"] in tempSession:
                del tempSession[session["stamp"]]
            del session["stamp"]
        if "table" in session:
            print("table in")
            del session["table"]
        if "result" in session:
            print("result deleted")
            del session["result"]
        session["stamp"] = stamp
        start = time.time()
        print(request.form["db"])
        session["dbChoice"] = int(request.form["db"])
        print("DB choice:", session["dbChoice"])
        global data
        # data = EnsemblRelease(session["dbChoice"])
        records = []
        file = request.files["efile"]
        file.name = file.filename
        file = BufferedReader(file)
        file = TextIOWrapper(file)
        tempSession[str(stamp) + "progress"] = 0
        vcf_reader = vcf.Reader(file)
        count = 0
        for record in vcf_reader:
            records.append(record)
            count += 1
        tempSession[str(stamp) + "count"] = count

        # print(type(file))
        @copy_current_request_context
        def annotatefunc(records, stamp):
            try:
                global tempSession
                mainKeys = '''
                            <label for="Cosmic">
                                <input type="checkbox" id="Cosmic" onclick="changeSelectText(this.parentElement)" onchange="resetSubkeys(this.parentElement,this)" />
                               Cosmic
                            </label>
                            <label for="Civic">
                                <input type="checkbox" id="Civic" onclick="changeSelectText(this.parentElement)" onchange="resetSubkeys(this.parentElement,this)" />
                               Civic
                            </label>
                            <label for="Cadd">
                                <input type="checkbox" id="Cadd" onclick="changeSelectText(this.parentElement)" onchange="resetSubkeys(this.parentElement,this)" />
                               Cadd
                            </label>
                            <label for="General">
                                <input type="checkbox" id="General" onclick="changeSelectText(this.parentElement)" onchange="resetSubkeys(this.parentElement,this)" />
                               General
                            </label>
                            '''
                newtable = manager.dict()
                tempSession[str(stamp) + "progress"] = manager.Value('progress', 0)
                nnewtable = {}
                pool = Pool(os.cpu_count())
                count = 0
                for record in records:
                    pool.apply_async(processVCFRecord, (record, count, newtable, tempSession[str(stamp) + "progress"]))
                    count += 1
                pool.close()
                pool.join()
                for key in newtable:
                    nnewtable[key] = {}
                    for key2 in newtable[key]:
                        nnewtable[key][key2] = []
                        for element in newtable[key][key2]:
                            nnewtable[key][key2].append(element)
                session["table"] = nnewtable.copy()
                newtablehtml = ""
                newtablehtmlheader = """<thead><tr>"""
                newtablehtmlheader += "<th></th>"
                newtablehtmlheader += "<th>Row Index</th>"
                keys = {
                    "General": {},
                    "Civic": {},
                    "Cosmic": {},
                    "Cadd": {},
                }
                newtablehtmlbody = "<tbody>"
                c = 0
                keyc = 2
                rowc = 0
                addedkeys = []
                popupdata = {}
                for item1 in list(nnewtable.items()):
                    lenn = len(list(item1[1].items()))
                    for subitem in zip(*item1[1].values()):
                        newtablehtmlbody += """<tr><td>
                            <br>
                            <button onclick="toggle(this)" style="color:white;font-size:20px;" name="+" class="btn btn-success btn-lg">
                            +
                            </button></td><td>%s</td>""" % item1[0]
                        for i in range(lenn):
                            mainKey = list(keys.keys())[i]
                            for item in list(subitem[i].items()):
                                val = item[1]
                                if type(val) == list and (type(val[0]) is OrderedDict or type(val[0]) is dict):
                                    newtablehtmlbody += "<td>%s</td>" % (
                                        getInnerAndHeaderHtmls(val,
                                                               "%s-%s-%s" % (item[0], rowc, c),
                                                               popupdata))
                                    newtablehtmlheader += "<th>%s</th>" % item[0] if item[0] not in addedkeys else ""
                                    if item[0] not in addedkeys:
                                        keyc = addHeaderKeys(addedkeys, keys[mainKey], item[0], item[0], keyc)
                                elif type(val) == dict:
                                    for key in val:
                                        newtablehtmlbody += "<td>%s</td>" % (val[key])
                                        newtablehtmlheader += "<th>%s</th>" % key if key not in addedkeys else ""
                                        if key not in addedkeys:
                                            keyc = addHeaderKeys(addedkeys, keys[mainKey], item[0], key, keyc, True)
                                else:
                                    newtablehtmlbody += "<td>%s</td>" % (val)
                                    newtablehtmlheader += "<th>%s</th>" % item[0] if item[0] not in addedkeys else ""
                                    if item[0] not in addedkeys:
                                        keyc = addHeaderKeys(addedkeys, keys[mainKey], item[0], item[0], keyc)
                                c += 1
                        newtablehtmlbody += "</tr>"
                    rowc += 1

                c = 1
                newtablehtmlheader += "</tr></thead>"
                newtablehtmlbody += "</tbody>"
                newtablehtml = """<table id = "table" class="table table-bordered">%s%s</table>""" % (
                newtablehtmlheader, newtablehtmlbody)
                tempSession[stamp] = (newtablehtml, mainKeys, keys, nnewtable, popupdata)
                tempSession[str(stamp) + "done"] = "done"
            except:
                print(traceback.format_exc())

        thread = threading.Thread(target=annotatefunc, args=[records, stamp])
        thread.start()
        print("time passed:", time.time() - start)
        return Response(response=json.dumps({"status": 1}))
    except:
        if "stamp" in session:
            if session["stamp"] in tempSession:
                del tempSession[session["stamp"]]
            del session["stamp"]
        if "table" in session:
            print("table in")
            del session["table"]
        if "result" in session:
            print("result deleted")
            del session["result"]
        print(traceback.format_exc())
        return Response(response=json.dumps({"status": 0}))


if __name__ == "__main__":
    app.run()
