from cyvcf2 import Variant
from rdflib import URIRef, Graph, BNode, RDF, RDFS, Literal, XSD

from iri_utils import get_position_iri, get_chromosome_landmark_iri, get_variant_iri
from namespaces import GFVO_NS, DCTERMS_NS, FALDO_NS

IDENTIFIER_CLASS = URIRef(GFVO_NS + 'Identifier')
ALIAS_CLASS = URIRef(GFVO_NS + 'Alias')
LABEL_CLASS = URIRef(GFVO_NS + 'Label')
HAS_VALUE_PROP = URIRef(GFVO_NS + "hasValue")
HAS_ATTRIBUTE_PROP = URIRef(GFVO_NS + "hasAttribute")
IDENTIFIER_PROP = URIRef(DCTERMS_NS + 'identifier')
EXACT_POSITION_CLASS = URIRef(FALDO_NS + "ExactPosition")
REFERENCE_PROP = URIRef(FALDO_NS + 'reference')
POSITION_PROP = URIRef(FALDO_NS + 'position')
LOCATION_PROP = URIRef(FALDO_NS + 'location')
HAS_PART_PROP = URIRef(GFVO_NS + "hasPart")
LOCUS_CLASS = URIRef(GFVO_NS + "Locus")
IS_LOCATED_ON_PROP = URIRef(GFVO_NS + 'isLocatedOn')


def add_identifier(graph: Graph, instance: URIRef, identifier: str):
    identifier_bnode = BNode()
    graph.add((identifier_bnode, RDF.type, IDENTIFIER_CLASS))
    graph.add((identifier_bnode, HAS_VALUE_PROP, Literal(identifier)))
    graph.add((instance, HAS_ATTRIBUTE_PROP, identifier_bnode))

    # add dcterms:identifier too
    graph.add((instance, IDENTIFIER_PROP, Literal(identifier)))


def add_label(graph: Graph, instance: URIRef, label: str):
    label_bnode = BNode()
    graph.add((label_bnode, RDF.type, LABEL_CLASS))
    graph.add((label_bnode, HAS_VALUE_PROP, Literal(label)))
    graph.add((instance, HAS_ATTRIBUTE_PROP, label_bnode))

    # add rdfs:label too
    graph.add((instance, RDFS.label, Literal(label)))

def add_location(graph: Graph, variant: Variant):
    chromosome_instance = get_chromosome_landmark_iri(variant)
    variant_instance = get_variant_iri(variant)

    locus_instance = BNode()
    graph.add((locus_instance, RDF.type, LOCUS_CLASS))

    position_instance = BNode()
    graph.add((position_instance, RDF.type, EXACT_POSITION_CLASS))
    graph.add((position_instance, POSITION_PROP, Literal(variant.start, datatype=XSD.int)))

    graph.add((locus_instance, HAS_ATTRIBUTE_PROP, position_instance))
    graph.add((locus_instance, IS_LOCATED_ON_PROP, position_instance))

    graph.add((variant_instance, HAS_PART_PROP, locus_instance))
    graph.add((variant_instance, IS_LOCATED_ON_PROP, chromosome_instance))

def add_alias(graph: Graph, instance: URIRef, alias: str):
    alias_bnode = BNode()
    graph.add((alias_bnode, RDF.type, ALIAS_CLASS))
    graph.add((alias_bnode, HAS_VALUE_PROP, Literal(alias)))
    graph.add((instance, HAS_ATTRIBUTE_PROP, alias_bnode))