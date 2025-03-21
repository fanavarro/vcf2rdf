# vcf2rdf
Variant Call File to Resource Description Framework transformation.

This python scripts takes a set of VCF files and generates an RDF file describing the variations found in the VCF files. The RDF file generated is compliant with the Genomic Feature and Variation Ontology (GFVO), an ontology for describing genomic files such as VCF, GFF, GTF, etc. The folder [GFVO_docs](./GFVO_docs) contains documentation on this ontology provided by the authors, including the article, supplementary material, an RDF example and the ontology itself.

## Usage

```
python vcf2rdf [-h] --input INPUT [INPUT ...] --output OUTPUT
               [--threads THREADS]

options:
  -h, --help            show this help message and exit
  --input INPUT [INPUT ...], -i INPUT [INPUT ...]
                        Input VCF files.
  --output OUTPUT, -o OUTPUT
                        Output RDF file.
  --threads THREADS, -t THREADS
                        Output RDF file.
```

Examples of use:

  - `python vcf2rdf -i vcf_folder/*.vcf -o graph.ttl -t4`: Process all vcf files in the folder `vcf_folder` and store the resulting RDF graph into `graph.ttl`. Use 4 threads.
  - `python vcf2rdf -i file.vcf -o graph.ttl -t4`: Process the file `file.vcf` and store the resulting RDF graph into `graph.ttl`. Use 4 threads.
  - `python vcf2rdf -i file1.vcf file2.vcf -o graph.ttl -t4`: Process the files `file1.vcf` and `file2.vcf` and store the resulting RDF graph into `graph.ttl`. Use 4 threads.


## Example of queries of the resulting graph
The rdf graph resulting of translating VCF to RDF can be uploaded in a graph database, such as graphdb or virtuoso, and queried by using SPARQL query language. Some examples are shown next.


### Get the variant count per sample.
```
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX gfvo: <https://github.com/BioInterchange/Ontologies/gfvo#>
PREFIX so: <http://purl.obolibrary.org/obo/SO_>
select ?variant (count(?sample) as ?count)
where {
    ?variant a so:0001059 .
    ?sample a gfvo:Sample .
    ?variant gfvo:hasEvidence [
            gfvo:hasSource ?sample 
    ] .
} group by ?variant order by desc(?count)
```

### Get the variant, variant identifier, samples having the variant, genotype of the samples, and the sequences of the genotypes.
```
PREFIX dc: <http://purl.org/dc/elements/1.1/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX gfvo: <https://github.com/BioInterchange/Ontologies/gfvo#>
PREFIX so: <http://purl.obolibrary.org/obo/SO_>
select *
where {
    ?variant a so:0001059 ;
    	dc:identifier ?variantID .
    ?variant gfvo:hasEvidence [
        gfvo:hasSource ?sample ;
        gfvo:hasAttribute [
        	a gfvo:Genotype ;
    		rdfs:label ?genotypeLabel ;
    		gfvo:hasFirstPart [ gfvo:hasValue ?firstGenotypeSeq ] ;
    		gfvo:hasLastPart [ gfvo:hasValue ?secondGenotypeSeq ] 
        ]
    ] .

    ?sample a gfvo:Sample .
}
```

### Get the variant, variant identifier, samples having the variant, genotype of the samples, the sequences of the genotypes, and info for the sample (genotype quality, read depth, allele frequency and mutated reads depth)

```
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
PREFIX dc: <http://purl.org/dc/elements/1.1/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX gfvo: <https://github.com/BioInterchange/Ontologies/gfvo#>
PREFIX so: <http://purl.obolibrary.org/obo/SO_>
select *
where {
    ?variant a so:0001059 ;
    	dc:identifier ?variantID .
    ?variant gfvo:hasEvidence [
        gfvo:hasSource ?sample ;
        gfvo:hasAttribute [
        	a gfvo:Genotype ;
    		rdfs:label ?genotypeLabel ;
    		gfvo:hasFirstPart [ gfvo:hasValue ?firstGenotypeSeq ] ;
    		gfvo:hasLastPart [ gfvo:hasValue ?secondGenotypeSeq ] ;
            gfvo:hasAttribute [
                    a gfvo:ConditionalGenotypeQuality ;
            		gfvo:hasValue ?genotypeQuality
            ] ;
            gfvo:hasAttribute [
                    a gfvo:Coverage ;
             		gfvo:hasValue ?readDepth
             ] ;
            gfvo:hasAttribute [
                    a gfvo:AlleleFrequency ;
            		gfvo:hasValue ?alternativeFrequency
            ]
        ]
    ] .
    BIND(xsd:int(round(xsd:double(?readDepth) * ?alternativeFrequency)) AS ?variantReadDepth ).
    ?sample a gfvo:Sample .
}
```