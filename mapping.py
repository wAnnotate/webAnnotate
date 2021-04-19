import requests

server = "http://rest.ensembl.org"


# ext = "/map/human/GRCh37/X:1000000..1000100:1/GRCh38?"
def remap(asm_from, asm_to, chromosome, location):  # Remaps from one assembly to another
    if asm_from == asm_to:
        return chromosome, location
    ext = "/map/human/%s/%s:%s..%s/%s?" % (asm_from, chromosome, location, location + 1, asm_to)
    # print("Remapping...")
    r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        return None, None

    decoded = r.json()
    # print("org:", decoded["mappings"][0]["original"])
    # print("mapped:", decoded["mappings"][0]["mapped"])
    try:
        return decoded["mappings"][0]["mapped"]["seq_region_name"], decoded["mappings"][0]["mapped"]["start"]
    except:
        print(decoded)
        return None, None


"""
remapped = remap("NCBI36", "GRCh38", 16, 50745926)
print(remapped)
"""
