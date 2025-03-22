from rdflib import BNode
from rdflib.plugins.stores import sparqlstore
from rdflib.plugins.stores.sparqlstore import _node_to_sparql


def node_to_sparql_bnode_extension (node):
    if isinstance(node, BNode):
        return '<bnode:b%s>' % node
    return _node_to_sparql(node)

def get_sparql_update_store(query_endpoint, update_endpoint, auth = None):
    return sparqlstore.SPARQLUpdateStore(
        query_endpoint=query_endpoint,
        update_endpoint=update_endpoint,
        node_to_sparql=node_to_sparql_bnode_extension,
        auth=auth)