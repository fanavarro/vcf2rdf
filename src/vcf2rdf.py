import argparse
import concurrent.futures

from rdflib import Graph, URIRef, RDF, Literal, RDFS, XSD, BNode
from cyvcf2 import VCF, Variant
from rdflib.namespace import NamespaceManager, Namespace

from iri_utils import *
from namespaces import *

SAMPLE_CLASS = URIRef(GFVO_NS + 'Sample')
CHROMOSOME_CLASS = URIRef(GFVO_NS + 'Chromosome')
SEQUENCE_ALTERATION_CLASS = URIRef(SEQUENCE_ONTOLOGY_NS + "0001059")
POSITION_CLASS = URIRef(FALDO_NS + "Position")
REFERENCE_SEQUENCE_CLASS = URIRef(GFVO_NS + "ReferenceSequence")
SEQUENCE_VARIANT_CLASS = URIRef(GFVO_NS + "SequenceVariant")
GENOTYPE_CLASS = URIRef(GFVO_NS + "Genotype")
HOMOZYGOUS_CLASS =URIRef(GFVO_NS + "Homozygous")
HETEROZYGOUS_CLASS = URIRef(GFVO_NS + "Heterozygous")
HEMIZYGOUS_CLASS = URIRef(GFVO_NS + "Hemizygous")

IDENTIFIER_PROP = URIRef(DCTERMS_NS + 'identifier')
LENGTH_PROP = URIRef(UNIPROT_NS + 'length')
REFERENCE_PROP = URIRef(FALDO_NS + 'reference')
POSITION_PROP = URIRef(FALDO_NS + 'position')
LOCATION_PROP = URIRef(FALDO_NS + 'location')
HAS_VALUE_PROP = URIRef(GFVO_NS + "hasValue")
HAS_ATTRIBUTE_PROP = URIRef(GFVO_NS + "hasAttribute")
HAS_QUALITY_PROP = URIRef(GFVO_NS + "hasQuality")
HAS_FIRST_PART_PROP = URIRef(GFVO_NS + 'hasFirstPart')
HAS_LAST_PART_PROP = URIRef(GFVO_NS + 'hasLastPart')
HAS_EVIDENCE_PROP = URIRef(GFVO_NS + 'hasEvidence')
HAS_SOURCE_PROP = URIRef(GFVO_NS + 'hasSource')

def create_arg_parser():
    parser = argparse.ArgumentParser(prog='vcf2rdf')
    parser.add_argument('--input', '-i', required=True, nargs='+', help='Input VCF files.')
    parser.add_argument('--output', '-o', required=True, help='Output RDF file.')
    parser.add_argument('--threads', '-t', required=False, default=1, type=int, help='Output RDF file.')
    return parser.parse_args()


def include_samples(vcf_info: VCF, graph: Graph) -> list:
    sample_list = []
    for sample in vcf_info.samples:
        sample_instance = URIRef(SAMPLE_NS + sample)
        graph.add((sample_instance, RDF.type, SAMPLE_CLASS))
        graph.add((sample_instance, IDENTIFIER_PROP, Literal(sample)))
        sample_list.append({'sample_instance': sample_instance,
                            'sample_name': sample})
    return sample_list


def include_ref_seqs(vcf_info: VCF, graph: Graph):
    seqlens = vcf_info.seqlens
    seqnames = vcf_info.seqnames
    for i in range(len(seqnames)):
        seqname = seqnames[i]
        seqlen = seqlens[i]
        chromosome_instance = URIRef(CHROMOSOME_NS + seqname)
        graph.add((chromosome_instance, RDF.type, CHROMOSOME_CLASS))
        graph.add((chromosome_instance, RDFS.label, Literal(seqname)))
        graph.add((chromosome_instance, LENGTH_PROP, Literal(seqlen, datatype=XSD.int)))




def include_variant(variant: Variant, sample_list: list, graph: Graph):
    variant_id = get_variant_id(variant)
    variant_instance = get_variant_iri(variant)
    graph.add((variant_instance, RDF.type, SEQUENCE_ALTERATION_CLASS))

    chromosome_instance = URIRef(CHROMOSOME_NS + variant.CHROM)

    position_instance = URIRef(POSITION_NS + f"{variant.CHROM}-{variant.start}")
    graph.add((position_instance, RDF.type, POSITION_CLASS))
    graph.add((position_instance, REFERENCE_PROP, chromosome_instance))
    graph.add((position_instance, POSITION_PROP, Literal(variant.start, datatype=XSD.int)))
    graph.add((position_instance, RDFS.label, Literal(f"{variant.CHROM}:{variant.start}")))

    graph.add((variant_instance, LOCATION_PROP, position_instance))
    graph.add((variant_instance, IDENTIFIER_PROP, Literal(variant.ID)))

    # Use of reference sequence to describe REF field
    reference_sequence_instance = URIRef(SEQUENCE_REFERENCE_NS + variant_id)
    graph.add((reference_sequence_instance, RDF.type, REFERENCE_SEQUENCE_CLASS))
    graph.add((reference_sequence_instance, HAS_VALUE_PROP, Literal(variant.REF)))
    graph.add((variant_instance, HAS_ATTRIBUTE_PROP, reference_sequence_instance))

    # Use of sequence variant to describe ALT field
    for allele in variant.ALT:
        allele_id = get_allele_id(variant, allele)
        sequence_variant_instance = URIRef(SEQUENCE_VARIANT_NS + allele_id)
        graph.add((sequence_variant_instance, RDF.type, SEQUENCE_VARIANT_CLASS))
        graph.add((sequence_variant_instance, HAS_VALUE_PROP, Literal(allele)))
        graph.add((variant_instance, HAS_ATTRIBUTE_PROP, sequence_variant_instance))

    # Include information about samples
    include_sample_genotypes(variant, sample_list, graph)


def get_zygosity_class(gt_type: int):
    # HOM_REF=0, HET=1. For gts012=True HOM_ALT=2, UNKNOWN=3
    if gt_type == 0:
        return HOMOZYGOUS_CLASS
    if gt_type == 1:
        return HETEROZYGOUS_CLASS
    if gt_type == 2:
        return HOMOZYGOUS_CLASS
    return None


def get_genotype_string(genotype):
    if genotype is None or len(genotype) != 3:
        return None
    return f"{genotype[0]}{'|' if genotype[2] else '/'}{genotype[1]}"


def include_sample_genotypes(variant: Variant, sample_list: list, graph: Graph):
    variant_id = get_variant_id(variant)
    variant_instance = get_variant_iri(variant)

    for i in range(len(sample_list)):
        sample_instance = sample_list[i]['sample_instance']
        sample_name = sample_list[i]['sample_name']
        gt_type = variant.gt_types[i] # HOM_REF=0, HET=1. For gts012=True HOM_ALT=2, UNKNOWN=3
        gt_phase = variant.gt_phases[i] # boolean indicating whether each sample is phased
        genotype = variant.genotypes[i] # [0, 1, True] corresponds to 0|1 while [1, 2, False] corresponds to 1/2

        genotype_instance = None
        zygosity = get_zygosity_class(gt_type)
        if zygosity is not None:
            genotype_instance = URIRef(GENOTYPE_NS + sample_name + '-' + variant_id)
            graph.add((genotype_instance, RDF.type, GENOTYPE_CLASS))
            graph.add((genotype_instance, HAS_QUALITY_PROP, zygosity))
            genotype_string = get_genotype_string(genotype)
            if genotype_string is not None:
                graph.add((genotype_instance, RDFS.label, Literal(genotype_string)))

            first_sequence_genotype = get_genotype_sequence_instance(genotype[0], variant)
            graph.add((genotype_instance, HAS_FIRST_PART_PROP, first_sequence_genotype))
            second_sequence_genotype = get_genotype_sequence_instance(genotype[1], variant)
            graph.add((genotype_instance, HAS_LAST_PART_PROP, second_sequence_genotype))

        evidence = BNode()
        graph.add((evidence, HAS_SOURCE_PROP, sample_instance))
        if genotype_instance is not None:
            graph.add((evidence, HAS_ATTRIBUTE_PROP, genotype_instance))
        graph.add((variant_instance, HAS_EVIDENCE_PROP, evidence))


def include_prefixes(graph: Graph):
    nm = NamespaceManager(graph, bind_namespaces="rdflib")
    nm.bind(prefix='sample', namespace=Namespace("https://namespaces.org/sample#"), override=True,
            replace=True)
    nm.bind(prefix='chr', namespace=Namespace("https://namespaces.org/chromosome#"), override=True,
            replace=True)
    nm.bind(prefix='var', namespace=Namespace("https://namespaces.org/variant#"), override=True,
            replace=True)
    nm.bind(prefix='pos', namespace=Namespace("https://namespaces.org/position#"), override=True,
            replace=True)
    nm.bind(prefix='alt', namespace=Namespace("https://namespaces.org/alt#"), override=True,
            replace=True)
    nm.bind(prefix='ref', namespace=Namespace("https://namespaces.org/ref#"), override=True,
            replace=True)
    nm.bind(prefix='gt', namespace=Namespace("https://namespaces.org/gt#"), override=True,
            replace=True)
    nm.bind(prefix='faldo', namespace=Namespace("http://biohackathon.org/resource/faldo#"), override=True,
            replace=True)
    nm.bind(prefix='obo', namespace=Namespace("http://purl.obolibrary.org/obo/"), override=True, replace=True)
    nm.bind(prefix='gfvo', namespace=Namespace("https://github.com/BioInterchange/Ontologies/gfvo#"), override=True, replace=True)
    nm.bind(prefix='up', namespace=Namespace('http://purl.uniprot.org/core/'), override=True, replace=True)

def generate_graph_for_vcf(vcf_file, threads: int) -> Graph:
    vcf_graph: Graph = Graph()
    vcf_info = VCF(vcf_file, strict_gt=True, threads=threads)
    sample_list = include_samples(vcf_info, vcf_graph)
    include_ref_seqs(vcf_info, vcf_graph)
    for variant in vcf_info:
        include_variant(variant, sample_list, vcf_graph)
    return vcf_graph

def generate_rdf(vcf_files: list, output_rdf_file: str, threads: int):
    graph: Graph = Graph()

    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        future_dict = {}
        for vcf_file in vcf_files:
            future_dict[vcf_file] = executor.submit(generate_graph_for_vcf, vcf_file, 1)

        for vcf_file, future in future_dict.items():
            try:
                vcf_graph = future.result()
            except Exception as exc:
                print(f'Error transforming {vcf_file}: {str(exc)}')
            else:
                print(f"{vcf_file} transformed")
                graph = graph + vcf_graph

    include_prefixes(graph)
    graph.serialize(output_rdf_file, format="ttl")

if __name__ == '__main__':
    arguments = create_arg_parser()
    generate_rdf(arguments.input, arguments.output, arguments.threads)