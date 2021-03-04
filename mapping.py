import requests
import sys

server = "http://rest.ensembl.org"


# ext = "/map/human/GRCh37/X:1000000..1000100:1/GRCh38?"
def remap(asm_from, asm_to, chromosome, location):  # Remaps from one assembly to another
    ext = "/map/human/%s/%s:%s..%s/%s?" % (asm_from, chromosome, location, location + 1, asm_to)
    print("Remapping...")
    r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()
    return decoded["mappings"][0]["mapped"]

# remapped = remap("NCBI36", "GRCh38", 16, 50745926)
# print(remapped)