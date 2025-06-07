import csv
import os
from pathlib import Path

gene_biotypes_dictionary_path = Path(os.path.dirname(os.path.abspath(__file__)), "resources/chromosome_identifiers.tsv")


class ChromosomeLookup:
    def __init__(self):
        self.dictionary = dict()
        with open(gene_biotypes_dictionary_path, newline='') as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')
            for row in reader:
                sequence_name = row["Sequence name"]
                genbank_name = row["GenBank seq accession"]
                ucsc_name = row["UCSC style name"]
                refseq_identifier = row["RefSeq seq accession"]

                self.dictionary[sequence_name] = refseq_identifier
                self.dictionary[genbank_name] = refseq_identifier
                self.dictionary[ucsc_name] = refseq_identifier


    def get_refseq_identifier(self, chromosome_name: str) -> str:
        return self.dictionary[chromosome_name]