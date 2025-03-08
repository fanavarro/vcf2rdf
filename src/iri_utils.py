import urllib
import uuid
import hashlib
from cyvcf2 import Variant
from rdflib import URIRef

from namespaces import *


def create_uuid_from_string(val: str):
    hex_string = hashlib.md5(val.encode("UTF-8")).hexdigest()
    return str(uuid.UUID(hex=hex_string))

def get_variant_id(variant: Variant):
    variant_id = f"{variant.CHROM}-{variant.start}-{variant.REF}-{urllib.parse.quote_plus("|".join(variant.ALT))}"
    return create_uuid_from_string(variant_id)

def get_variant_iri(variant: Variant):
    return URIRef(VARIANT_NS + get_variant_id(variant))

def get_allele_id(variant: Variant, allele_string):
    variant_id = get_variant_id(variant)
    allele_id = ""
    if '<' in allele_string:
        variant_size = variant.end - variant.start
        allele_id = urllib.parse.quote_plus(f"{variant_id}_{allele_string.replace('<', ''). replace('>', '')}{str(variant_size)}")
    else:
        allele_id = urllib.parse.quote_plus(variant_id + '_' + allele_string)

    return create_uuid_from_string(allele_id)


def get_genotype_sequence_instance(genotype_index, variant):
    if genotype_index == 0:
        return URIRef(SEQUENCE_REFERENCE_NS + get_variant_id(variant))

    allele_id = get_allele_id(variant, variant.ALT[genotype_index - 1])
    return URIRef(SEQUENCE_VARIANT_NS + allele_id)
