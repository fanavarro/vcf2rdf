import csv
import os
from pathlib import Path

from rdflib import URIRef

from namespaces import SEQUENCE_ONTOLOGY_NS

transcript_biotypes_dictionary_path = Path(os.path.dirname(os.path.abspath(__file__)), "resources/transcript_biotypes.tsv")
gene_biotypes_dictionary_path = Path(os.path.dirname(os.path.abspath(__file__)), "resources/gene_biotypes.tsv")


class BiotypeLookup:
    def __init__(self, lookup_file: Path):
        self.dictionary = dict()
        with open(lookup_file, newline='') as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')
            for row in reader:
                biotype_name = row["name"]
                biotype_group = row["biotype_group"]
                biotype_so_acc = row["so_acc"]
                so_iri = URIRef(SEQUENCE_ONTOLOGY_NS + biotype_so_acc.replace('SO:', ''))
                self.dictionary[biotype_name] = so_iri

    def get_biotype_sequence_ontology(self, biotype_name):
        return self.dictionary.get(biotype_name)

    def contains_biotype(self, biotype_name: str):
        return biotype_name in self.dictionary

class TranscriptBiotypeLookup(BiotypeLookup):
    def __init__(self):
        super().__init__(transcript_biotypes_dictionary_path)

class GeneBiotypeLookup(BiotypeLookup):
    def __init__(self):
        super().__init__(gene_biotypes_dictionary_path)