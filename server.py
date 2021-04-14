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
from json2html import *
from civicdb import CivicDb
from cosmicdb import CosmicDb
import traceback
import time
import mapping  # mapping.remap()
import universalKeys as dictKeys  # selected dict keys and their descriptions
import cadd

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
civic = CivicDb()
cosmic = CosmicDb()
manager = Manager()
app.config.from_object(__name__)
Session(app)
biothingsAssembly = "hg19"




@app.route("/")
def index():
    return render_template("index.html")


def getGeneInfo(gene_id, table):
    geneData = gene_client.getgene(gene_id)
    if "Summary" not in geneData:
        table["Summary"].append("No data avaliable")
    else:
        table["Summary"].append(geneData["Summary"])
    """
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
               "start", "end", "strand", "genome", "summary", "clingen", "entrezgene",
               "variantdata"]
        table = session["table"].copy()
        tablehtml = """<table id = "table" class="table table-bordered"><thead><tr>"""
        for th in ths:
            if th != "variants":
                tablehtml += "<th>%s</th>" % th
        tablehtml += "</tr></thead><tbody>"
        count = len(list(table.values())[0])
        for c in range(count):
            tablehtml += "<tr>"
            for th in ths:
                if th != "variants":
                    tablehtml += "<td>%s</td>" % table[th][c]
            tablehtml += "</tr>"
        tablehtml += "</tbody></table>"
        return render_template("annotated.html", table=tablehtml)
    return redirect("/")


def getGeneFromLocation(chr, pos):  # Gets location of gene, returns a Gene object (v102, v75, v54)
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
    return gene


def getGeneFromGeneId(gId):  # Gets Ensembl id, returns a Gene object
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
    return gene


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
    rowid = int(rowid)
    # print(session["table"]["entrezgene"][rowid])
    if "ncbi" in session["table"]["entrezgene"][rowid]:
        gene = session["table"]["entrezgene"][rowid].split("href=\"")[1].split("\"")[0].split("/")[-1]
        if requests.get(
                "https://www.ncbi.nlm.nih.gov/projects/Gene/download_expression.cgi?PROJECT_DESC=PRJEB4337&GENE=%s" % gene).status_code == 500:
            return render_template("visualization.html", table="No data available", dta="[]")
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


def processVariantData(variant, count, hgsvs, index):
    variantdata = None
    if str(variant["_id"]) in hgsvs:
        variantdata = '<option value="%s-%s-%s">%s</option>' % (index, count, "biothings", variant["_id"])
        observed = "observed" if variant["observed"] else "not observed"
        varDict = {
            "header":"Biothings: " + str(variant["_id"] + " (" + observed + ")")
        }
        for key in variant:
            if variant[key] and key not in ["_id","_version","chrom","hg19","observed"]:
                if type(variant[key]) != dict or all(isinstance(x, float) for x in list(variant[key].values())):
                    if len(list(variant[key].values())) > 2 and all(isinstance(x, float) for x in list(variant[key].values())):
                        variant[key]["graph"] = """
                                                    <button onclick="showGraph('%s',this)">Show Graph</button>
                                                """ % (key)
                    varDict[key] = json2html.convert(json=variant[key])
                elif type(variant[key]) == dict:
                    varDict[key] = {}
                    for key2 in variant[key]:
                        if type(variant[key][key2]) == dict and all(isinstance(x, float) for x in list(variant[key][key2].values())):
                            variant[key][key2]["graph"] = """
                                                    <button onclick="showGraph('%s---%s',this)">Show Graph</button>
                                                """ % (key2,key)
                        varDict[key][key2] = json2html.convert(json=variant[key][key2],escape = False)
                        varDict[key][key2] += "<div id='chart%s---%s'></div>" % (key2,key)
        return json.dumps(varDict), variantdata
    return None, None


def getVariantData(id, assembly, index, count, subdict, hgsvs, variantdata):
    variant = variant_client.getvariant(id, assembly=assembly)
    if type(variant) == dict and "_id" in variant:
        vardict, vardata = processVariantData(variant, count, hgsvs, index)
        if vardict:
            variantdata += vardata
            subdict["listofvariants"].append(vardict)
            count += 1
    elif type(variant) == list:
        for var in variant:
            if "_id" not in var:
                continue
            print(str(var["_id"]))
            vardict, vardata = processVariantData(var, count, hgsvs, index)
            if vardict:
                variantdata += vardata
                subdict["listofvariants"].append(vardict)
                count += 1
    return count, variantdata

def formatListOfLists(listofdicts):
    for dictionary in listofdicts:
        keys = list(dictionary.keys())
        for key in keys:
            if not dictionary[key]:
                del dictionary[key]
    return listofdicts

def processVCFRecord(record, table, index, nnewtable):
    print("new record")
    foundGene = False
    gene_dict = {}
    main_sub_dict = {
        index:{
            "Civic":[],
            "Cosmic":[],
            "Summary":[],
            "Expression":[],
            "Cadd":[],
            "Clingen":[],
            "dbnsfp":[],
            "dbsnp":[],
            "clinvar":[],
            "evs":[],
            "mutdb":[],
            "gwassnps":[],
            "docm":[],
            "snpedia":[],
            "emv":[],
            "wellderly":[],
            "exac":[],
            "grasp":[]
        }
    }
    
    subdict = {
        " ": [],
        "expression": [],
        "gene_id": [],
        "gene_name": [],
        "biotype": [],
        "contig": [],
        "start": [],
        "end": [],
        "variants": [],
        "strand": [],
        "genome": [],
        "summary": [],
        "clingen": [],
        "entrezgene": [],
        "variantdata": [],
        "listofvariants": [],
        "listofvariantscivic": [],
        "listofvariantscosmic": []
    }

    civic_variants_template = {}
    for key in dictKeys.civicVariants:
        civic_variants_template[dictKeys.civicDesc(key)] = "No data available"
    cosmic_CMC_template = {}
    for key in dictKeys.cosmicCMC:
        cosmic_CMC_template[dictKeys.cosmicDesc(key)] = "No data available"

    if record.ID:  # RsId exists
        # print("rsid exists")
        try:
            gene, variant = getGeneFromRsId(record.ID)
            gene_dict = gene.__dict__
            foundGene = True
            getGeneInfo(gene.gene_id, main_sub_dict[index])
        except Exception as e:
            print("getGeneFromRsId: ", e)
            foundGene = False

    if not foundGene:
        # if not record.ID: print("rsid does not exist")
        try:
            gene = getGeneFromLocation(record.CHROM, record.POS)
            gene_dict = gene[0].__dict__
            foundGene = True
            getGeneInfo(gene[0].gene_id, main_sub_dict[index])
        except Exception as e:
            print("getGeneFromLocation: ", e)
            foundGene = False

    if foundGene:
        try:
            count = 0
            mappedChr, mappedPos = mapping.remap(dbName[str(session["dbChoice"])], "GRCh37", record.CHROM, record.POS)
            varianthtml = "No data available"
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
            if record.ALT:
                for i in record.ALT:
                    if len(str(i)) == 1 and len(str(record.REF)) == 1:
                        cadd_data = cadd.getSNV("GRCh37",mappedChr,mappedPos,record.REF,str(i))
                        if cadd_data and type(cadd_data) == list:
                            for d in cadd_data:
                                cadd_dict = {}
                                for key in cadd.keys:
                                    if key in d:
                                        cadd_dict[key.replace("-"," ").replace("_"," ")] = d[key]
                                    else:
                                        cadd_dict[key] = "No data available"
                                main_sub_dict[index]["Cadd"].append(cadd_dict)
        except:
            print(traceback.format_exc())
        try:
            print("hgsvs ", index)
            if record.ID:
                count, variantdata = getVariantData(record.ID, biothingsAssembly, index, count, subdict, hgsvs,
                                                    variantdata)
            elif hgsvs:
                for hgsv in hgsvs:
                    count, variantdata = getVariantData(hgsv, biothingsAssembly, index, count, subdict, hgsvs,
                                                        variantdata)
        except:
            print(traceback.format_exc())
        try:
            civicdata = civic.findVariantsFromLocation(mappedChr, mappedPos).copy()
            if civicdata:
                count = 0
                for var in civicdata:
                    variant_groups = formatListOfLists(civic.findVariantGroups(var["variant_groups"])).copy()
                    assertions = formatListOfLists(civic.findAssertions(var["variant_id"])).copy()
                    clinical_significances = formatListOfLists(civic.findClinicalEvidences(var["variant_id"])).copy()
                    civicdatagene = civic.findGeneFromLocation(mappedChr,mappedPos).copy()
                    civicdatagene = formatListOfLists([civicdatagene])
                    for cs in clinical_significances:
                        keys = list(cs.keys())
                        for key in keys:
                            if key not in dictKeys.civicClinicalEvidences:
                                del cs[key]
                            else:
                                data = cs[key]
                                cs[dictKeys.civicDesc(key)] = data
                                del cs[key]
                    
                    for cs in assertions:
                        keys = list(cs.keys())
                        for key in keys:
                            if key not in dictKeys.civicAssertions:
                                del cs[key]
                            else:
                                data = cs[key]
                                cs[dictKeys.civicDesc(key)] = data
                                del cs[key]
                    for cs in variant_groups:
                        keys = list(cs.keys())
                        for key in keys:
                            if key not in dictKeys.civicVariantGroups:
                                del cs[key]
                            else:
                                data = cs[key]
                                cs[dictKeys.civicDesc(key)] = data
                                del cs[key]
                    for cs in civicdatagene:
                        keys = list(cs.keys())
                        for key in keys:
                            if key not in dictKeys.civicGenes:
                                del cs[key]
                            else:
                                data = cs[key]
                                cs[dictKeys.civicDesc(key)] = data
                                del cs[key]
                    """
                    variant_groups = json2html.convert(json = variant_groups, escape = False)
                    assertions = json2html.convert(json = assertions, escape = False)
                    clinical_significances = json2html.convert(json = clinical_significances, escape = False)
                    genehtml = json2html.convert(json = civicdatagene, escape=False)
                    """
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
                    var = temp
                    """
                    html = json2html.convert(json=var)
                    """
                    html = var
                    CivicDict = {
                        "Variants":html if html else civic_variants_template,
                        "Variant Groups":variant_groups if variant_groups else "No data available",
                        "Genes":genehtml if genehtml else "No data avaliable",
                        "Assertions":assertions if assertions else "No data available",
                        "Clinical Evidences":clinical_significances if clinical_significances else "No data available",
                    }
                    main_sub_dict[index]["Civic"].append(CivicDict)
                    count += 1
        except:
            print(traceback.format_exc())
        try:
            cosmicdata = cosmic.findVariantsFromLocation("GRCh37",mappedChr,mappedPos)
            if not cosmicdata:
                mappedChr, mappedPos = mapping.remap(dbName[str(session["dbChoice"])], "GRCh38", record.CHROM, record.POS)
                cosmicdata = cosmic.findVariantsFromLocation("GRCh38",mappedChr,mappedPos)
            if cosmicdata:
                count = 0
                for row in cosmicdata:
                    resistanceMutationsHtml = ""
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
                        """
                        if resistanceMutations:
                            resistanceMutationsHtml = json2html.convert(json=resistanceMutations)
                        """
                    keys = list(row.keys())
                    for key in keys:
                        if (key not in dictKeys.cosmicCMC):
                            del row[key]
                        elif not row[key]:
                            row[key] = "No data available"
                        else:
                            row[dictKeys.cosmicDesc(key)] = row[key]
                            del row[key]
                    temp = {}
                    for key in list(cosmic_CMC_template.keys()):
                        if key in row:
                            temp[key] = row[key]
                        else:
                            temp[key] = "No data available"
                    row = temp
                        
                    variantData = json2html.convert(json = row, escape=False)
                    variantData = row
                    CosmicDict = {
                        "CMC":variantData if variantData else cosmic_CMC_template,
                        "Resistance Mutations":resistanceMutations if resistanceMutations else "No data available"
                    }
                    main_sub_dict[index]["Cosmic"].append(CosmicDict)
                    count += 1
        except Exception as exp:
            print(index, "- variant exp: ", exp)
            print(traceback.format_exc())
        if variantdata:
            varianthtml = '<select onchange="toggleModal(this)"><option value=""></option>%s</select>' % variantdata
        print("-------")
        greatest_len = 0
        for key in ["Cosmic","Civic","Cadd"]:
            if len(main_sub_dict[index][key]) > greatest_len:
                greatest_len = len(main_sub_dict[index][key])
        for key in ["Cosmic","Civic","Cadd"]:
            while greatest_len > len(main_sub_dict[index][key]):
                if key == "Cosmic":
                    main_sub_dict[index][key].append(
                        {
                            "CMC":cosmic_CMC_template,
                            "Resistance Mutations":"No data availabe"
                        }
                    )
                elif key == "Cadd":
                    temp = {}
                    for key in cadd.keys:
                        temp[key] = "No data avilable"
                    main_sub_dict[index][key].append(temp)
                else:
                    main_sub_dict[index][key].append(
                        {
                            "Variants": civic_variants_template,
                            "Variant Groups":"No data available",
                            "Genes":"No data available",
                            "Assertions":"No data available",
                            "Clinical Evidences":"No data available",
                        }
                    )
                

        subdict["variantdata"].append(varianthtml)
        # print(gene)
        for key in subdict.keys():
            if key in gene_dict.keys() and key not in ["summary", "clingen", "entrezgene", "rowid", "expression",
                                                       "variants", "listofvariants",
                                                       "variantdata", " ", "listofvariantscivic","listofvariantscosmic"]:
                subdict[key].append(str(gene_dict[key]))
            elif key not in ["summary", "clingen", "entrezgene", "rowid", "expression", "variants", "variantdata", " ",
                             "listofvariants", "listofvariantscivic","description","listofvariantscosmic"]:
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
    if len(main_sub_dict[index]["Cosmic"]) > 0:
        nnewtable[str(index)] = {
            "Cosmic":main_sub_dict[index]["Cosmic"],
            "Civic":main_sub_dict[index]["Civic"],
            "Cadd":main_sub_dict[index]["Cadd"]
        }

@app.route("/annotate", methods=["POST"])
def annotate():
    start = time.time()
    session["dbChoice"] = int(request.form["db"])
    print("DB choice:", session["dbChoice"])
    global data
    data = EnsemblRelease(session["dbChoice"])
    file = request.files["efile"]
    file.name = file.filename
    file = BufferedReader(file)
    file = TextIOWrapper(file)
    # print(type(file))
    vcf_reader = vcf.Reader(file)
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
                ''' 
    subKeys = {
        "Cosmic":["CMC-6","Resistance Mutations-7"],
        "Civic":["Variants-1","Variant Groups-2","Genes-3","Assertions-4", "Clinical Evidences-5"],
        "Cadd": ["ConsScore-8","mirSVR Score-9","dbscSNV ada score-10","dbscSNV rf score-11","RawScore-12","PHRED-13"]
    }

    allKeys = {
        "Civic-Variants": [dictKeys.civicDesc(k) for k in dictKeys.civicVariants],
        "Civic-Variant Groups": [],
        "Civic-Genes":[dictKeys.civicDesc(k) for k in dictKeys.civicGenes],
        "Civic-Assertions":[],
        "Civic-Clinical Evidences":[],
        "Cosmic-CMC": [dictKeys.cosmicDesc(k) for k in dictKeys.cosmicCMC],
        "Cosmic-Resistance Mutations": [],
        "Cadd-ConsScore":[],
        "Cadd-mirSVR Score":[],
        "Cadd-dbscSNV ada score":[],
        "Cadd-dbscSNV rf score":[],
        "Cadd-RawScore":[],
        "Cadd-PHRED":[]
    }


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
        "entrezgene": [],
        "variants": [],
        "variantdata": [],
        "listofvariants": [],
        "listofvariantscivic": [],
        "listofvariantscosmic": []
    }
    newtable = manager.dict()
    nnewtable = {}

    pool = Pool(os.cpu_count())
    count = 0
    processes = []
    for record in vcf_reader:
        pool.apply_async(processVCFRecord, (record, ttable, count, newtable))
        count += 1
    pool.close()
    pool.join()
    count = len(list(ttable.keys()))
    for key in newtable:
        nnewtable[key] = {}
        for key2 in newtable[key]:
            nnewtable[key][key2] = []
            for element in newtable[key][key2]:
                nnewtable[key][key2].append(element)
            print(key," :",key2," :",len(nnewtable[key][key2]))
    session["table"] = table.copy()
    newtablehtml = ""
    newtablehtmlheader = """<thead><tr>"""
    newtablehtmlheader += "<th>Row Index</th>"
    """
    for key in allKeys:
        if len(allKeys[key]):
            for key2 in allKeys[key]:
                newtablehtmlheader += "<th>%s</th>" % key2
        else:
            newtablehtmlheader += "<th>%s</th>" % key
    """

    keys = {
        "Civic":{},
        "Cosmic":{},
        "Cadd":{}
    }
    newtablehtmlbody = "<tbody>"
    c = 0
    keyc = 1
    rowc = 0
    addedkeys = []
    popupdata = {}
    for item1 in list(nnewtable.items()):
        for subitem in zip(item1[1]["Cosmic"],item1[1]["Civic"],item1[1]["Cadd"]):
            newtablehtmlbody += "<tr><td>%s</td>" % item1[0]
            for item in list(subitem[1].items()):
                val = item[1]
                if type(val) == list and type(val[0]) == dict:
                    innerhtml = "<option value="" selected>""</option>"
                    for element in val:
                        innerhtml += """
                                        <option value="%s-%s-%s-%s">%s</option>
                                     """ % (item[0],rowc,c,list(element.values())[0],list(element.values())[0] )
                        popupdata["%s-%s-%s-%s" % (item[0],rowc,c,list(element.values())[0])] = json2html.convert(json = element)
                    innerhtmlselect = """<select onchange="toggleModal(this)">%s</select>
                                    """ % (innerhtml)
                    newtablehtmlbody += "<td>%s</td>" % (innerhtmlselect)
                    newtablehtmlheader += "<th>%s</th>" % item[0] if item[0] not in addedkeys else ""
                    if item[0] not in addedkeys:
                        addedkeys.append(item[0])
                        keys["Civic"][item[0]] = keyc
                        keyc += 1
                elif type(val) == dict:
                    for key in val:
                        newtablehtmlbody += "<td>%s</td>" % (val[key])
                        newtablehtmlheader += "<th>%s</th>" % key if key not in addedkeys else ""
                        if key not in addedkeys:
                            addedkeys.append(key)
                            if item[0] not in keys["Civic"]:
                                keys["Civic"][item[0]] = {}
                            keys["Civic"][item[0]][key] = keyc
                            keyc += 1
                else:
                    newtablehtmlbody += "<td>%s</td>" % (val)
                    newtablehtmlheader += "<th>%s</th>" % item[0] if item[0] not in addedkeys else ""
                    if item[0] not in addedkeys:
                        addedkeys.append(item[0])
                        keys["Civic"][item[0]] = keyc
                        keyc += 1

                c += 1
            for item in list(subitem[0].items()):
                val = item[1]
                if type(val) == dict:
                    for key in val:
                        newtablehtmlbody += "<td>%s</td>" % (val[key])
                        newtablehtmlheader += "<th>%s</th>" % key if key not in addedkeys else ""
                        if key not in addedkeys:
                            addedkeys.append(key)
                            if item[0] not in keys["Cosmic"]:
                                keys["Cosmic"][item[0]] = {}
                            keys["Cosmic"][item[0]][key] = keyc
                            keyc += 1
                else:
                    newtablehtmlbody += "<td>%s</td>" % (val)
                    newtablehtmlheader += "<th>%s</th>" % item[0] if item[0] not in addedkeys else ""
                    if item[0] not in addedkeys:
                        addedkeys.append(item[0])
                        keys["Cosmic"][item[0]] = keyc
                        keyc += 1

                c += 1
            for item in list(subitem[2].items()):
                val = item[1]
                newtablehtmlbody += "<td>%s</td>" % (val)
                newtablehtmlheader += "<th>%s</th>" % item[0] if item[0] not in addedkeys else ""
                if item[0] not in addedkeys:
                    addedkeys.append(item[0])
                    keys["Cadd"][item[0]] = keyc
                    keyc += 1
                c += 1
            newtablehtmlbody += "</tr>"
        rowc += 1

    c = 1
    newtablehtmlheader += "</tr></thead>"
    newtablehtmlbody += "</tbody>"
    newtablehtml = """<table id = "table" class="table table-bordered">%s%s</table>""" % (newtablehtmlheader,newtablehtmlbody)


    print("time passed:",time.time()-start)
    return render_template("annotated.html", table=newtablehtml, mainKeys = mainKeys, subkeys = json.dumps(keys), allKeys = json.dumps(allKeys), allData = nnewtable,popupdata = popupdata)


if __name__ == "__main__":
    app.run()
