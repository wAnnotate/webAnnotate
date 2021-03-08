from pyensembl import EnsemblRelease

data54 = EnsemblRelease(54)
data75 = EnsemblRelease(75)
data102 = EnsemblRelease(102)


def genesAtLocus(version, contig, position):
    if int(version) == 54:
        genes = data54.genes_at_locus(contig=contig, position=position)
    elif int(version) == 75:
        genes = data75.genes_at_locus(contig=contig, position=position)
    elif int(version) == 102:
        genes = data102.genes_at_locus(contig=contig, position=position)
    else:
        return Exception("Wrong ensembl db version.")
    if not genes:
        return None
    return genes


def geneById(version, gId):
    if int(version) == 54:
        gene = data54.gene_by_id(gene_id=gId)
    elif int(version) == 75:
        gene = data75.gene_by_id(gene_id=gId)
    elif int(version) == 102:
        gene = data102.gene_by_id(gene_id=gId)
    else:
        return Exception("Wrong ensembl db version.")
    if not gene:
        return None
    return gene
