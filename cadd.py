import requests

server = "https://cadd.gs.washington.edu/api/v1.0/"

keys = ["ConsScore","mirSVR-Score","dbscSNV-ada_score","dbscSNV-rf_score","RawScore","PHRED"]

def getSNV(version, chrom, pos, ref=None, alt=None):  # GRCh37, X, 10000, A, T
    chrom = str(chrom)
    pos = str(pos)
    ext = version + "-v1.6_inclAnno/" + str(chrom) + ":" + str(pos)
    if ref and alt:
        ext += "_" + ref + "_" + alt
    print(server + ext)
    r = requests.get(server + ext, headers={"Content-Type": "application/json"})
    if not r.ok:  # Bad response
        return None
    decoded = r.json()
    if not decoded:  # Empty list
        return None
    return decoded  # List of SNVs []


"""snvList = getSNV("GRCh37", "22", 44044020)
for snv in snvList:
    for val in snv.items():
        print(val[0], ": ", val[1])
    print()"""
