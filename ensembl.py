import requests

server = "https://rest.ensembl.org"  # GRCh38, 102
server37 = "https://grch37.rest.ensembl.org"  # GRCh37, 75


def dict2obj(d):
    if isinstance(d, list):
        d = [dict2obj(x) for x in d]
    if not isinstance(d, dict):
        return d

    class C:
        pass

    obj = C()
    for k in d:
        obj.__dict__[k] = dict2obj(d[k])
    return obj


def returner(url):  # Do not use this!
    r = requests.get(url, headers={"Content-Type": "application/json"})
    if not r.ok:
        return Exception("Not found.")
    decoded = r.json()
    if not decoded:
        return Exception("Not found.")
    return dict2obj(decoded)


def getGeneFromGeneId(gId):
    """Retrieves features (e.g. genes, transcripts, variants and more) that overlap a region defined by the given
    identifier.

    https://rest.ensembl.org/documentation/info/overlap_id

    :param gId An Ensembl stable ID; ENSG00000157764"""

    ext = f"/overlap/id/{gId}?feature=gene"
    try:
        decoded = returner(server + ext)
        if len(decoded) == 1:
            return decoded[0]
        for g in decoded:
            if g.id == gId:
                return g
        return Exception("Not found.")
    except Exception as e:
        return Exception("Not found.")


def getGenesFromLocation(chr, pos, db):
    """Retrieves features (e.g. genes, transcripts, variants and more) that overlap a given region.

    https://rest.ensembl.org/documentation/info/overlap_region

    :param chr Chromosome
    :param pos Position
    :param db Genome version; 75 or 102"""

    ext = f"/overlap/region/human/{chr}:{pos}-{pos}?feature=gene"
    if db == 75:
        return returner(server37 + ext)
    elif db == 102:
        return returner(server + ext)
    return Exception("Select either 75 or 102.")


def getVariantFromLocation(chr, pos, ref, alt, db):  # Returns dict, not object!
    """Retrieves features (e.g. genes, transcripts, variants and more) that overlap a given region.

    https://rest.ensembl.org/documentation/info/overlap_region

    :param chr Chromosome
    :param pos Position
    :param ref Reference
    :param alt Alternative
    :param db Genome version, 75 or 102"""

    ext = f"/overlap/region/human/{chr}:{pos}-{pos}?feature=variation"
    if db == 75:
        url = server37 + ext
    elif db == 102:
        url = server + ext
    else:
        return Exception("Select either 75 or 102.")

    r = requests.get(url, headers={"Content-Type": "application/json"})
    if not r.ok:
        return Exception("Not found.")
    decoded = r.json()
    if not decoded:
        return Exception("Not found.")
    if len(decoded) == 1:
        return decoded[0]
    if len(decoded) == 0:
        return None

    chr = str(chr)
    pos = int(pos)
    ref = str(ref)
    alt = str(alt)
    for v in decoded:
        if v["start"] == pos and ref in v["alleles"] and alt in v["alleles"]:
            return v
    return None


def getVEPFromId(varId):  # dbSNP, COSMIC, HGMD
    """Fetch variant consequences based on a variant identifier.

    https://rest.ensembl.org/documentation/info/vep_id_get

    :param varId Query ID. Supports dbSNP, COSMIC and HGMD identifiers; rs56116432 COSM476"""

    ext = f"/vep/human/id/{varId}?"
    return returner(server + ext)


def getVEPFromLocation(chr, pos, alt, db):
    """Fetch variant consequences

    https://rest.ensembl.org/documentation/info/vep_region_get

    :param chr Chromosome
    :param pos Position
    :param alt Variation allele; C DUP
    :param db Genome version; 75 or 102"""

    ext = f"/vep/human/region/{chr}:{pos}-{pos}/{alt}?"
    if db == 75:
        return returner(server37 + ext)
    elif db == 102:
        return returner(server + ext)
    return Exception("Select either 75 or 102.")


"""variant = getVariantFromLocation("X", 47426121, "C", "G", 75)
print(variant)"""
