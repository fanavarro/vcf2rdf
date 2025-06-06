import concurrent.futures
import os
import traceback
from pathlib import Path

from rdflib import Graph, RDF, Literal, XSD, BNode, RDFS
from cyvcf2 import VCF
from rdflib.namespace import NamespaceManager, Namespace
from iri_utils import *
from namespaces import *
from gfvo_utils import add_identifier, add_label, add_location, add_alias
from info_ann_parser import parse_info_ann
from biotypes import TranscriptBiotypeLookup, GeneBiotypeLookup
from ontology_utils import OntologyLookUp

ANNOTATION_FIELD = 'ANN'
INFORMATION_CONTENT_ENTITY_CLASS = URIRef(GFVO_NS + 'InformationContentEntity')
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
TRANSCRIPT_CLASS = URIRef(SEQUENCE_ONTOLOGY_NS + '0000673')
GENE_CLASS = URIRef(SEQUENCE_ONTOLOGY_NS + '0000704')
FUNCTIONAL_ANNOTATION_CLASS = URIRef(GFVO_ANN_EXT_NS + 'FunctionalAnnotation')
IMPACT_CLASS = URIRef(GFVO_ANN_EXT_NS + 'Impact')
GFVO_FEATURE_CLASS = URIRef(GFVO_NS + 'Feature')

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
REFERS_TO_PROP = URIRef(GFVO_NS + 'refersTo')
HAS_ANNOTATION_PROP = URIRef(GFVO_NS + 'hasAnnotation')
ASSOCIATED_WITH_PROP = URIRef(SEQUENCE_ONTOLOGY_NS + 'associated_with')

transcript_biotype_lookup = TranscriptBiotypeLookup()
gene_biotype_lookup = GeneBiotypeLookup()
sequence_ontology_lookup = OntologyLookUp(Path(os.path.dirname(os.path.abspath(__file__)), 'resources/so.owl'))


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

def add_transcript_instance(transcript_id: str, biotype_name: str, graph: Graph) -> URIRef:
    transcript_instance = get_transcript_iri(transcript_id)
    graph.add((transcript_instance, RDF.type, TRANSCRIPT_CLASS))
    if transcript_biotype_lookup.contains_biotype(biotype_name):
        transcript_class = transcript_biotype_lookup.get_biotype_sequence_ontology(biotype_name)
        graph.add((transcript_instance, RDF.type, transcript_class))

    biotype_instance = URIRef(str(transcript_instance) + '/attribute/biotype')
    graph.add((biotype_instance, RDF.type, INFORMATION_CONTENT_ENTITY_CLASS))
    graph.add((biotype_instance, HAS_VALUE_PROP, Literal(biotype_name, datatype=XSD.string)))
    graph.add((transcript_instance, HAS_ATTRIBUTE_PROP, biotype_instance))

    return transcript_instance

def add_gene_instance(gene_id, biotype_name, graph):
    gene_instance = get_gene_iri(gene_id)
    graph.add((gene_instance, RDF.type, GENE_CLASS))

    if biotype_name is not None:
        if gene_biotype_lookup.contains_biotype(biotype_name):
            gene_class = gene_biotype_lookup.get_biotype_sequence_ontology(biotype_name)
            graph.add((gene_instance, RDF.type, gene_class))

        biotype_instance = URIRef(str(gene_instance) + '/attribute/biotype')
        graph.add((biotype_instance, RDF.type, INFORMATION_CONTENT_ENTITY_CLASS))
        graph.add((biotype_instance, HAS_VALUE_PROP, Literal(biotype_name, datatype=XSD.string)))
        graph.add((gene_instance, HAS_ATTRIBUTE_PROP, biotype_instance))
    return gene_instance

def get_ann_gene_instance(annotation, graph):
    if annotation.get('Gene_Name'):
        return get_gene_iri(annotation.get('Gene_Name'))
    else:
        return None

def get_ann_feature_instance(annotation, graph) -> URIRef:
    feature_instance = None
    feature_type = annotation.get('Feature_Type')
    # print(feature_type)
    feature_id = annotation.get('Feature_ID')

    if feature_type == 'transcript':
        biotype = annotation.get('Transcript_BioType')
        feature_instance = add_transcript_instance(transcript_id=feature_id, biotype_name=biotype, graph=graph)

    elif feature_type == 'intergenic_region':
        pass
    elif feature_type == 'gene':
        biotype = annotation.get('Gene_BioType')
        feature_instance = add_gene_instance(gene_id=feature_id, biotype_name=biotype, graph=graph)
    elif feature_type == 'gene_variant':
        # Gene variant has a notation like MIR3675.3, where the original gene is MIR3675
        # In the genome annotations, the "." is replaced by "-", but some gene variants are not included.
        biotype = annotation.get('Gene_BioType')
        feature_instance = add_gene_instance(gene_id=feature_id.replace('.', '-'), biotype_name=biotype, graph=graph)

    if feature_instance is not None:
        graph.add((feature_instance, RDF.type, GFVO_FEATURE_CLASS))
    #else:
        #print(f'Feature type "{feature_type}" not recognized. Feature ID = {feature_id}')

    return feature_instance

def add_annotation(variant: Variant, annotation: dict, graph: Graph):
    variant_instance = get_variant_iri(variant)
    functional_annotation = BNode()
    graph.add((functional_annotation, RDF.type, FUNCTIONAL_ANNOTATION_CLASS))

    # Add annotation as simple text by annotation properties
    for key, value in annotation.items():
        annotation_property_iri = get_annotation_property_iri(key)
        graph.add((functional_annotation, annotation_property_iri, Literal(value, datatype=XSD.string)))
        graph.add((annotation_property_iri, RDFS.label, Literal(key, datatype=XSD.string)))

    # Try to represent the annotation by linking transcript, genes, etc...
    # Get instance of the affected feature (transcript or gene or intergenic region...)
    feature_instance = get_ann_feature_instance(annotation, graph)
    if feature_instance is not None:
        graph.add((functional_annotation, REFERS_TO_PROP, feature_instance))


    gene_instance = get_ann_gene_instance(annotation, graph)
    if gene_instance is not None:
        graph.add((functional_annotation, REFERS_TO_PROP, gene_instance))

    # Get the allele of this effect
    allele_instance = get_allele_iri(variant, annotation.get('Allele'))
    if allele_instance is not None:
        graph.add((functional_annotation, REFERS_TO_PROP, allele_instance))

    # Include the annotation (effect or consequence)
    so_terms = sequence_ontology_lookup.get_iris_from_label(annotation.get('Annotation'))
    if so_terms is not None and len(so_terms) > 0:
        so_term = so_terms[0]
        graph.add((variant_instance, RDF.type, so_term))
        graph.add((functional_annotation, ASSOCIATED_WITH_PROP, so_term))

    # Include the impact
    impact = annotation.get('Annotation_Impact')
    if impact is not None:
        impact_node = BNode()
        graph.add((impact_node, RDF.type, IMPACT_CLASS))
        graph.add((impact_node, HAS_VALUE_PROP, Literal(impact)))
        graph.add((functional_annotation, HAS_ATTRIBUTE_PROP, impact_node))

    # Add annotation to variant
    graph.add((variant_instance, HAS_ANNOTATION_PROP, functional_annotation))


def add_annotations(variant: Variant, ann_field_info, graph: Graph):
    ann_string = variant.INFO.get(ANNOTATION_FIELD)
    annotations = parse_info_ann(ann_string, ann_field_info)
    for annotation in annotations:
        add_annotation(variant, annotation, graph)

def include_variant(variant: Variant, ann_field_info, sample_list: list, graph: Graph):
    variant_instance = get_variant_iri(variant)
    graph.add((variant_instance, RDF.type, SEQUENCE_ALTERATION_CLASS))
    add_identifier(graph, variant_instance, variant.ID)
    add_location(graph, variant)
    add_sv_type(graph, variant)
    add_ref_field(graph, variant)
    add_alt_field(graph, variant)
    add_sample_genotypes(variant, sample_list, graph)
    add_annotations(variant, ann_field_info, graph)


def add_alt_field(graph, variant):
    # Use of sequence variant to describe ALT field
    variant_instance = get_variant_iri(variant)
    for allele in variant.ALT:
        sequence_variant_instance = get_allele_iri(variant, allele)
        graph.add((sequence_variant_instance, RDF.type, SEQUENCE_VARIANT_CLASS))
        graph.add((sequence_variant_instance, HAS_VALUE_PROP, Literal(allele)))
        graph.add((variant_instance, HAS_ATTRIBUTE_PROP, sequence_variant_instance))


def add_ref_field(graph, variant):
    # Use of reference sequence to describe REF field
    variant_instance = get_variant_iri(variant)
    reference_sequence_instance = get_ref_sequence_iri(variant)
    graph.add((reference_sequence_instance, RDF.type, REFERENCE_SEQUENCE_CLASS))
    graph.add((reference_sequence_instance, HAS_VALUE_PROP, Literal(variant.REF)))
    graph.add((variant_instance, HAS_ATTRIBUTE_PROP, reference_sequence_instance))


def add_sv_type(graph, variant):
    variant_instance = get_variant_iri(variant)
    sv_type_class = get_sv_type_class(variant)
    if sv_type_class is not None:
        graph.add((variant_instance, RDF.type, sv_type_class))


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

def add_sample_genotypes(variant: Variant, sample_list: list, graph: Graph):
    # Include information about samples
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

def get_ann_field_header(vcf: VCF):
    try:
        ann_field_info = vcf.get_header_type(ANNOTATION_FIELD)
    except KeyError:
        ann_field_info = None

    return ann_field_info


def generate_graph_for_vcf(vcf_file, threads: int) -> Graph:
    print(f"Processing {vcf_file}")
    vcf_graph: Graph = Graph()
    vcf = VCF(vcf_file, strict_gt=True, threads=threads)
    ann_field_info = get_ann_field_header(vcf)
    sample_list = include_samples(vcf, vcf_graph)
    include_ref_seqs(vcf, vcf_graph)
    for variant in vcf:
        include_variant(variant, ann_field_info, sample_list, vcf_graph)
    return vcf_graph

def generate_rdf_file_for_vcf(vcf_file: Path, output_folder: Path, threads: int):
    output_file = Path(output_folder, vcf_file.stem + ".ttl")
    graph: Graph = generate_graph_for_vcf(vcf_file, threads = 1)
    include_prefixes(graph)
    graph.serialize(output_file, format="ttl")
    del graph
    return output_file


def generate_rdf_per_vcf(vcf_files: list, output_folder: Path, threads: int):
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        future_dict = {}
        for vcf_file in vcf_files:
            if output_folder is None:
                output_folder = Path(vcf_file).parent
            future_dict[vcf_file] = executor.submit(generate_rdf_file_for_vcf,Path(vcf_file), output_folder, 1)

        for vcf_file, future in future_dict.items():
            try:
                rdf_file = future.result()
            except Exception as exc:
                print(f'Error processing {vcf_file}: {str(exc)}')
                print(traceback.format_exc())
            else:
                print(f"{vcf_file} processed. Results saved to {rdf_file}")


def generate_rdf(vcf_files: list, output_rdf_file: Path, threads: int):
    graph: Graph = Graph()
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        future_dict = {}
        for vcf_file in vcf_files:
            future_dict[vcf_file] = executor.submit(generate_graph_for_vcf, Path(vcf_file), 1)

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
