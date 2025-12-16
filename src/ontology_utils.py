from pathlib import Path

from rdflib import Graph, URIRef, RDFS, Literal, RDF


class OntologyLookUp:
    def __init__(self, ontology_path: Path):
        self.ontology = Graph()
        self.ontology.parse(ontology_path)

    def get_iris_from_label(self, label) -> URIRef:
        return list(self.ontology.subjects(predicate=RDFS.label, object=Literal(label), unique=True))

    def is_subclass_of(self, entity_iri: URIRef, parent_iri: URIRef) -> bool:
        query = "PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>\n"
        query += f'ASK WHERE {{ <{str(entity_iri)}> rdfs:subClassOf+ <{str(parent_iri)}> }}'
        return self.ontology.query(query).askAnswer