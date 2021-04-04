def desc(key):  # Convert underscore to space and capitalize each word
    description = ""
    for word in key.split('_'):
        description += word.capitalize()
        description += " "
    return description


civicVariants = {
    "variant": [],
    "summary": [],
    "variant_groups": [],
    "variant_types": [],
    "civic_variant_evidence_score": [],
    "variant_civic_url": [],
}
civicVariantGroups = {
    "variant_group": [],
    "description": [],
    "variant_group_civic_url": [],
}
civicGenes = {
    "name": [],
    "description": [],
    "gene_civic_url": [],
}
civicAssertions = {
    "assertion_type": [],
    "amp_category": [],
    "nccn_guideline": [],
    "regulatory_approval": [],
    "fda_companion_test": [],
    "assertion_summary": [],
    "assertion_description": [],
    "assertion_civic_url": [],
}
civicClinicalEvidences = {
    "disease": [],
    "phenotypes": [],
    "drugs": [],
    "drug_interaction_type": [],
    "evidence_type": [],
    "evidence_direction": [],
    "evidence_level": [],
    "clinical_significance": [],
    "evidence_statement": [],
    "citation_id": [],
    "source_type": [],
    "citation": [],
    "rating": [],
    "evidence_status": [],
    "variant_origin": [],
    "evidence_civic_url": [],
}
