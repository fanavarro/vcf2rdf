import argparse

from rdf_generator import generate_rdf, generate_rdf_per_vcf


def create_arg_parser():
    parser = argparse.ArgumentParser(prog='vcf2rdf')
    parser.add_argument('--input', '-i', required=True, nargs='+', help='Input VCF files.')
    parser.add_argument('--merged_output', '-o', required=False, help='Output RDF file. If set, all VCFs at input '
                                                               'will be merged into the same RDF file. If not, an RDF per '
                                                               'VCF will be created')
    parser.add_argument('--output_folder', '-d', required=False, help='Output folder for the RDF files.')
    parser.add_argument('--format', '-f', required=False, default='nt', help='Output format (nt, ttl, xml, pretty-xml, n3, trig, trix, json-ld)')
    parser.add_argument('--threads', '-t', required=False, default=1, type=int, help='Output RDF file.')
    return parser.parse_args()



if __name__ == '__main__':
    arguments = create_arg_parser()
    if arguments.merged_output:
        generate_rdf(arguments.input, arguments.merged_output, arguments.format, arguments.threads)
    else:
        generate_rdf_per_vcf(arguments.input, arguments.output_folder, arguments.format, arguments.threads)