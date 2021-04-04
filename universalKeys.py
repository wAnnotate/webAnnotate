# CIViC keys
civicVariants = [
    "summary",
    "civic_variant_evidence_score",
    "variant_civic_url",
]
civicVariantGroups = [
    "variant_group_civic_url",
]
civicGenes = [
    "gene_civic_url",
]
civicAssertions = [
    "amp_category",
    "nccn_guideline",
    "assertion_civic_url",
]
civicClinicalEvidences = [
    "disease",
    "drugs",
    "clinical_significance",
    "evidence_statement",
    "evidence_civic_url",
]


def civicDesc(key):  # Convert underscore to space and capitalize each word
    description = "CIViC:"
    for word in key.split('_'):
        description += " "
        description += word.capitalize()
    return description


# COSMIC keys
cosmicCMC = [
    "cosmic_sample_mutated",
    "gerp_rs",
    "min_sift_score",
    "min_sift_pred",
    "mutation_significance_tier",
    "mutation_url",
]
cosmicResistanceMutations = [
    "Drug Name",
]
cosmicCGC = []
cosmicHGNC = []

cosmicDescriptions = {}
f = open("cosmic/README.txt", "r")
f.readline()
f.readline()
file = f.read().split("\n")
for line in file:
    keyValue = line.split(" - ", 1)
    k = keyValue[0].lower()
    v = keyValue[1]
    v = v.split(" (", 1)[0].split(".", 1)[0].split(",", 1)[0]
    cosmicDescriptions[k] = v
f.close()


def cosmicDesc(key):
    description = "COSMIC: "
    if key in cosmicCMC:
        description += cosmicDescriptions[key]
    elif key in cosmicResistanceMutations:
        description += key
    return description
