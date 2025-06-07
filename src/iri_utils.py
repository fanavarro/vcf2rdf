import re
import urllib
import uuid
import hashlib
from cyvcf2 import Variant
from rdflib import URIRef

from namespaces import *
from chromosome_lookup import ChromosomeLookup

chromosome_lookup = ChromosomeLookup()

def create_iri_from_primary_key_string(ns: str, string: str) -> URIRef:
    return URIRef(ns + str(create_uuid_from_string(string)))

def create_uuid_from_string(val: str) -> str:
    hex_string = hashlib.md5(val.encode("UTF-8")).hexdigest()
    return str(uuid.UUID(hex=hex_string))

def get_variant_primary_key_string(variant: Variant) -> str:
    return f"{variant.CHROM}-{variant.start}-{variant.REF}-{urllib.parse.quote_plus("|".join(variant.ALT))}"

def get_variant_id(variant: Variant) -> str:
    variant_id = get_variant_primary_key_string(variant)
    return create_uuid_from_string(variant_id)

def get_variant_iri(variant: Variant) -> URIRef:
    variant_primary_key_string = get_variant_primary_key_string(variant)
    return create_iri_from_primary_key_string(VARIANT_NS, variant_primary_key_string)

def get_allele_primary_key_string(variant: Variant, allele_string) -> str:
    variant_id = get_variant_id(variant)
    if '<' in allele_string:
        variant_size = variant.end - variant.start
        return urllib.parse.quote_plus(
            f"{variant_id}_{allele_string.replace('<', '').replace('>', '')}{str(variant_size)}")
    else:
        return urllib.parse.quote_plus(variant_id + '_' + allele_string)

def get_allele_id(variant: Variant, allele_string) -> str:
    allele_primary_key_string = get_allele_primary_key_string(variant, allele_string)
    return create_uuid_from_string(allele_primary_key_string)

def get_allele_iri(variant: Variant, allele_string) -> URIRef:
    return URIRef(SEQUENCE_VARIANT_NS + get_allele_id(variant, allele_string))

def get_ref_sequence_primary_key_string(variant: Variant) -> str:
    return get_variant_primary_key_string(variant)

def get_ref_sequence_id(variant: Variant) -> str:
    ref_seq_primary_key_string = get_ref_sequence_primary_key_string(variant)
    return create_uuid_from_string(ref_seq_primary_key_string)

def get_ref_sequence_iri(variant: Variant) -> URIRef:
    ref_seq_primary_key_string = get_ref_sequence_primary_key_string(variant)
    return create_iri_from_primary_key_string(SEQUENCE_REFERENCE_NS, ref_seq_primary_key_string)

def get_genotype_sequence_instance(genotype_index, variant) -> URIRef:
    if genotype_index == 0:
        return get_ref_sequence_iri(variant)

    return get_allele_iri(variant, variant.ALT[genotype_index - 1])

def get_position_iri(variant: Variant) -> URIRef:
    return URIRef(POSITION_NS + f"{variant.CHROM}-{variant.start}")

def get_chromosome_landmark_iri(variant: Variant) -> URIRef:
    return get_landmark_iri_from_chromosome_name(variant.CHROM)

def get_landmark_iri_from_chromosome_name(chromosome_name: str) -> URIRef:
    return URIRef(LANDMARK_NS + chromosome_lookup.get_refseq_identifier(chromosome_name))

def get_chromosome_feature_iri(variant: Variant) -> URIRef:
    return get_chromosome_feature_iri_from_chromosome_name(variant.CHROM)

def get_chromosome_feature_iri_from_chromosome_name(chromosome_name: str) -> URIRef:
    return URIRef(FEATURE_NS + chromosome_lookup.get_refseq_identifier(chromosome_name))

def get_sample_iri(sample_name: str) -> URIRef:
    return URIRef(SAMPLE_NS + sample_name)

def get_genotype_iri(sample_name: str, variant: Variant) -> URIRef:
    variant_id = get_variant_id(variant)
    return URIRef(GENOTYPE_NS + sample_name + '-' + variant_id)

def get_transcript_iri(transcript_id: str) -> URIRef:
    return URIRef(TRANSCRIPT_NS + transcript_id)

def get_gene_iri(gene_id: str) -> URIRef:
    return URIRef(GENE_NS + gene_id)

def get_intergenic_region_iri(intergenic_region_name: str) -> URIRef:
    return URIRef(INTERGENIC_REGION_NS + replace_special_chars_for_iri(intergenic_region_name))

def get_annotation_property_iri(annotation_name: str) -> URIRef:
    return URIRef(GFVO_ANN_EXT_NS + replace_special_chars_for_iri(annotation_name))

def replace_special_chars_for_iri(iri: str) -> str:
    new_iri = (iri
               .replace('.', '_')
               .replace(' ','_')
               .replace('/', '_')
               .replace('"', '')
               .lower())
    return re.sub(r'_+', '_', new_iri)