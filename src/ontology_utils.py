from pathlib import Path

from rdflib import Graph, URIRef, RDFS, Literal


class OntologyLookUp:
    def __init__(self, ontology_path: Path):
        self.ontology = Graph()
        self.ontology.parse(ontology_path)

    def get_iris_from_label(self, label) -> URIRef:
        return list(self.ontology.subjects(predicate=RDFS.label, object=Literal(label), unique=True))