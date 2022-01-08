from biothings_client import get_client
import mapping

def getGeneFromLocation(dbId, chr, pos):  # Gets reference db and location of gene, returns a Gene object
    if dbId == 77 or dbId == 76 or dbId == 75:
        data = EnsemblRelease(dbId)
    else:
        return Exception("Wrong database id.")

    gene = data.genes_at_locus(contig=chr, position=pos)

    if not gene:
        return Exception("Gene not found.")

    return gene


gc = get_client('gene')
vc = get_client('variant')

contig = "4"
pos = 6843853
ref = "G"
alt = "A"

#contig, pos = mapping.remap("NCBI36", "GRCh37", contig, pos)

hgvsStr = vc.format_hgvs(contig, pos, ref, alt)
print(hgvsStr)
exVar = vc.getvariant(hgvsStr, assembly='hg19')

print()
print(exVar["cosmic"])
#print(exVar[1]["vcf"])

"""
exGen = gc.getgene("ENSG00000233684")
print(exGen["ensembl"])

print()

exGen = gc.getgene("ENSG00000132155", fields='summary,entrezgene,clingen')
print(exGen["ensembl"])


exVar = vc.getvariant("rs2066844", assembly="hg19")
print (exVar)
print()
if type(exVar) == list:
    print("list")
if type(exVar) == dict:
    print('dict')

exVar = vc.getvariant("rs6054257", fields="hits")
print(exVar)
exVar = vc.getvariant("rs6040355", fields="hits")
print(exVar)
exVar = vc.getvariant("microsat1", fields="hits")
print(exVar)
"""
