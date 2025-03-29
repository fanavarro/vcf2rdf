import argparse
import concurrent.futures
import traceback
from rdflib import Graph, URIRef, RDF, Literal, XSD, BNode
from cyvcf2 import VCF, Variant
from rdflib.namespace import NamespaceManager, Namespace
from iri_utils import *
from namespaces import *
from src.gfvo_utils import add_identifier, add_label, add_location

SAMPLE_CLASS = URIRef(GFVO_NS + 'Sample')
CHROMOSOME_CLASS = URIRef(SEQUENCE_ONTOLOGY_NS + '0000340')
LANDMARK_CLASS = URIRef(GFVO_NS + 'Landmark')
SEQUENCE_ALTERATION_CLASS = URIRef(SEQUENCE_ONTOLOGY_NS + "0001059")
POSITION_CLASS = URIRef(FALDO_NS + "Position")
REFERENCE_SEQUENCE_CLASS = URIRef(GFVO_NS + "ReferenceSequence")
SEQUENCE_VARIANT_CLASS = URIRef(GFVO_NS + "SequenceVariant")
GENOTYPE_CLASS = URIRef(GFVO_NS + "Genotype")
HOMOZYGOUS_CLASS =URIRef(GFVO_NS + "Homozygous")
HETEROZYGOUS_CLASS = URIRef(GFVO_NS + "Heterozygous")
HEMIZYGOUS_CLASS = URIRef(GFVO_NS + "Hemizygous")
CONDITIONAL_GENOTYPE_QUALITY_CLASS = URIRef(GFVO_NS + 'ConditionalGenotypeQuality')
COVERAGE_CLASS = URIRef(GFVO_NS + 'Coverage')
ALLELE_FREQUENCY_CLASS = URIRef(GFVO_NS + 'AlleleFrequency')
STRUCTURAL_DELETION_CLASS = URIRef(SEQUENCE_ONTOLOGY_NS + '0000159')
STRUCTURAL_INSERTION_CLASS = URIRef(SEQUENCE_ONTOLOGY_NS + '0000667')
STRUCTURAL_DUPLICATION_CLASS = URIRef(SEQUENCE_ONTOLOGY_NS + '1000035')
STRUCTURAL_INVERSION_CLASS = URIRef(SEQUENCE_ONTOLOGY_NS + '1000036')
STRUCTURAL_COPY_NUMBER_VARIATION_CLASS = URIRef(SEQUENCE_ONTOLOGY_NS + '0001019')
STRUCTURAL_BREAKEND_CLASS = URIRef(SEQUENCE_ONTOLOGY_NS + '0002062')

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
        sample_instance = get_sample_iri(sample)
        graph.add((sample_instance, RDF.type, SAMPLE_CLASS))
        # graph.add((sample_instance, IDENTIFIER_PROP, Literal(sample)))
        add_identifier(graph, sample_instance, sample)
        sample_list.append({'sample_instance': sample_instance,
                            'sample_name': sample})
    return sample_list


def include_ref_seqs(vcf_info: VCF, graph: Graph):
    seqlens = vcf_info.seqlens
    seqnames = vcf_info.seqnames
    for i in range(len(seqnames)):
        seqname = seqnames[i]
        seqlen = seqlens[i]
        chromosome_instance = get_landmark_iri_from_chromosome_name(seqname)
        graph.add((chromosome_instance, RDF.type, CHROMOSOME_CLASS))
        graph.add((chromosome_instance, RDF.type, LANDMARK_CLASS))
        # graph.add((chromosome_instance, RDFS.label, Literal(seqname)))
        add_label(graph, chromosome_instance, seqname)
        graph.add((chromosome_instance, LENGTH_PROP, Literal(seqlen, datatype=XSD.int)))

def get_sv_type_class(variant) -> URIRef:
    sv_type = variant.INFO.get('SVTYPE')
    if sv_type is not None:
        if sv_type == 'DEL':
            return STRUCTURAL_DELETION_CLASS
        if sv_type == 'INS':
            return STRUCTURAL_INSERTION_CLASS
        if sv_type == 'DUP':
            return STRUCTURAL_DUPLICATION_CLASS
        if sv_type == 'INV':
            return STRUCTURAL_INVERSION_CLASS
        if sv_type == 'CNV':
            return STRUCTURAL_COPY_NUMBER_VARIATION_CLASS
        if sv_type == 'BND':
            return STRUCTURAL_BREAKEND_CLASS
        print(f"SVTYPE '{sv_type}' not recognised")
    return None



def include_variant(variant: Variant, sample_list: list, graph: Graph):
    variant_instance = get_variant_iri(variant)
    graph.add((variant_instance, RDF.type, SEQUENCE_ALTERATION_CLASS))
    add_identifier(graph, variant_instance, variant.ID)
    add_location(graph, variant)

    sv_type_class = get_sv_type_class(variant)
    if sv_type_class is not None:
        graph.add((variant_instance, RDF.type, sv_type_class))

    # Use of reference sequence to describe REF field
    reference_sequence_instance = get_ref_sequence_iri(variant)
    graph.add((reference_sequence_instance, RDF.type, REFERENCE_SEQUENCE_CLASS))
    graph.add((reference_sequence_instance, HAS_VALUE_PROP, Literal(variant.REF)))
    graph.add((variant_instance, HAS_ATTRIBUTE_PROP, reference_sequence_instance))

    # Use of sequence variant to describe ALT field
    for allele in variant.ALT:
        sequence_variant_instance = get_allele_iri(variant, allele)
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

def get_ref_depth(variant:Variant, sample_index: int) -> int:
    gt_ref_depth = variant.gt_ref_depths[sample_index]
    if gt_ref_depth is not None and gt_ref_depth != -1:
        return gt_ref_depth

    ## SNIFFLES FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">
    if 'DR' in variant.FORMAT:
        return variant.format('DR')[sample_index][0]
    else:
        return None

def get_alt_depth(variant: Variant, sample_index: int) -> int:
    gt_alt_depth = variant.gt_alt_depths[sample_index]
    if gt_alt_depth is not None and gt_alt_depth != -1:
        return gt_alt_depth

    ## SNIFFLES FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">
    if 'DV' in variant.FORMAT:
        return variant.format('DV')[sample_index][0]
    else:
        return None


def get_depth(variant: Variant, sample_index: int) -> int:
    gt_depth = variant.gt_depths[sample_index]  # DP field

    if gt_depth is not None and gt_depth != -1:
        return gt_depth

    gt_alt_depth = get_alt_depth(variant, sample_index)
    gt_ref_depth = get_ref_depth(variant, sample_index)

    if gt_alt_depth is not None and gt_ref_depth is not None:
        return gt_alt_depth + gt_ref_depth

    return None

def get_alt_freq(variant: Variant, sample_index: int) -> int:
    gt_alt_freq = variant.gt_alt_freqs[sample_index] # AF field

    if gt_alt_freq is not None and gt_alt_freq != -1:
        return gt_alt_freq

    number_of_reads = get_depth(variant, sample_index)
    number_of_reads_supporting_variation = get_alt_depth(variant, sample_index)

    if number_of_reads is not None and number_of_reads_supporting_variation is not None and number_of_reads > 0:
        return float(number_of_reads_supporting_variation)/float(number_of_reads)

    return None

def include_sample_genotypes(variant: Variant, sample_list: list, graph: Graph):
    variant_id = get_variant_id(variant)
    variant_instance = get_variant_iri(variant)

    for i in range(len(sample_list)):
        sample_instance = sample_list[i]['sample_instance']
        sample_name = sample_list[i]['sample_name']
        gt_type = variant.gt_types[i] # HOM_REF=0, HET=1. For gts012=True HOM_ALT=2, UNKNOWN=3
        gt_phase = variant.gt_phases[i] # boolean indicating whether each sample is phased
        genotype = variant.genotypes[i] # [0, 1, True] corresponds to 0|1 while [1, 2, False] corresponds to 1/2
        gt_depth = get_depth(variant, i)
        gt_alt_freq = get_alt_freq(variant, i)
        gt_genotype_quality = variant.gt_quals[i] # GQ field


        genotype_instance = None
        zygosity = get_zygosity_class(gt_type)
        if zygosity is not None:
            genotype_instance = get_genotype_iri(sample_name, variant)
            graph.add((genotype_instance, RDF.type, GENOTYPE_CLASS))
            graph.add((genotype_instance, HAS_QUALITY_PROP, zygosity))
            genotype_string = get_genotype_string(genotype)
            if genotype_string is not None:
                # graph.add((genotype_instance, RDFS.label, Literal(genotype_string)))
                add_label(graph, genotype_instance, genotype_string)

            if genotype[0] != -1:
                first_sequence_genotype = get_genotype_sequence_instance(genotype[0], variant)
                graph.add((genotype_instance, HAS_FIRST_PART_PROP, first_sequence_genotype))
            if genotype[1] != -1:
                second_sequence_genotype = get_genotype_sequence_instance(genotype[1], variant)
                graph.add((genotype_instance, HAS_LAST_PART_PROP, second_sequence_genotype))

            if gt_genotype_quality is not None and gt_genotype_quality != -1:
                genotype_quality_instance = BNode()
                graph.add((genotype_quality_instance, RDF.type, CONDITIONAL_GENOTYPE_QUALITY_CLASS))
                graph.add((genotype_quality_instance, HAS_VALUE_PROP, Literal(gt_genotype_quality)))
                graph.add((genotype_instance, HAS_ATTRIBUTE_PROP, genotype_quality_instance))
            if gt_depth is not None and gt_depth != -1:
                coverage_instance = BNode()
                graph.add((coverage_instance, RDF.type, COVERAGE_CLASS))
                graph.add((coverage_instance, HAS_VALUE_PROP, Literal(gt_depth)))
                graph.add((genotype_instance, HAS_ATTRIBUTE_PROP, coverage_instance))
            if gt_alt_freq is not None and gt_alt_freq != -1:
                alt_freq_instance = BNode()
                graph.add((alt_freq_instance, RDF.type, ALLELE_FREQUENCY_CLASS))
                graph.add((alt_freq_instance, HAS_VALUE_PROP, Literal(gt_alt_freq)))
                graph.add((genotype_instance, HAS_ATTRIBUTE_PROP, alt_freq_instance))

        evidence = BNode()
        graph.add((evidence, HAS_SOURCE_PROP, sample_instance))
        if genotype_instance is not None:
            graph.add((evidence, HAS_ATTRIBUTE_PROP, genotype_instance))
        graph.add((variant_instance, HAS_EVIDENCE_PROP, evidence))


def include_prefixes(graph: Graph):
    nm = NamespaceManager(graph, bind_namespaces="rdflib")
    nm.bind(prefix='sample', namespace=Namespace(SAMPLE_NS), override=True,
            replace=True)
    nm.bind(prefix='landmark', namespace=Namespace(LANDMARK_NS), override=True,
            replace=True)
    nm.bind(prefix='var', namespace=Namespace(VARIANT_NS), override=True,
            replace=True)
    nm.bind(prefix='pos', namespace=Namespace(POSITION_NS), override=True,
            replace=True)
    nm.bind(prefix='alt', namespace=Namespace(SEQUENCE_VARIANT_NS), override=True,
            replace=True)
    nm.bind(prefix='ref', namespace=Namespace(SEQUENCE_REFERENCE_NS), override=True,
            replace=True)
    nm.bind(prefix='gt', namespace=Namespace(GENOTYPE_NS), override=True,
            replace=True)
    nm.bind(prefix='faldo', namespace=Namespace(FALDO_NS), override=True,
            replace=True)
    nm.bind(prefix='so', namespace=Namespace(SEQUENCE_ONTOLOGY_NS), override=True, replace=True)
    nm.bind(prefix='gfvo', namespace=Namespace(GFVO_NS), override=True, replace=True)
    nm.bind(prefix='up', namespace=Namespace(UNIPROT_NS), override=True, replace=True)

def generate_graph_for_vcf(vcf_file, threads: int) -> Graph:
    print(f"Processing {vcf_file}")
    vcf_graph: Graph = Graph()
    vcf_info = VCF(vcf_file, strict_gt=True, threads=threads)
    sample_list = include_samples(vcf_info, vcf_graph)
    include_ref_seqs(vcf_info, vcf_graph)
    for variant in vcf_info:
        include_variant(variant, sample_list, vcf_graph)
    return vcf_graph

def generate_rdf(vcf_files: list, output_rdf_file: str, threads: int):
    graph: Graph = Graph()
    #with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
    #with ThreadPoolExecutorStackTraced(max_workers=threads) as executor:
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        future_dict = {}
        for vcf_file in vcf_files:
            future_dict[vcf_file] = executor.submit(generate_graph_for_vcf, vcf_file, 1)

        for vcf_file, future in future_dict.items():
            try:
                vcf_graph = future.result()
            except Exception as exc:
                print(f'Error processing {vcf_file}: {str(exc)}')
                print(traceback.format_exc())
            else:
                print(f"{vcf_file} processed")
                graph = graph + vcf_graph

    include_prefixes(graph)
    print(f"Saving graph into {output_rdf_file}")
    graph.serialize(output_rdf_file, format="ttl")
    print("Done")

if __name__ == '__main__':
    arguments = create_arg_parser()
    generate_rdf(arguments.input, arguments.output, arguments.threads)