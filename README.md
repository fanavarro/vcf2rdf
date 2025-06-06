# vcf2rdf
Variant Call File to Resource Description Framework transformation.

This python scripts takes a set of VCF files and generates an RDF file describing the variations found in the VCF files. The RDF file generated is compliant with the Genomic Feature and Variation Ontology (GFVO), an ontology for describing genomic files such as VCF, GFF, GTF, etc. The folder [GFVO_docs](./GFVO_docs) contains documentation on this ontology provided by the authors, including the article, supplementary material, an RDF example and the ontology itself.

## Usage

```
usage: vcf2rdf [-h] --input INPUT [INPUT ...] [--merged_output MERGED_OUTPUT]
               [--output_folder OUTPUT_FOLDER] [--threads THREADS]

options:
  -h, --help            show this help message and exit
  --input INPUT [INPUT ...], -i INPUT [INPUT ...]
                        Input VCF files.
  --merged_output MERGED_OUTPUT, -o MERGED_OUTPUT
                        Output RDF file. If set, all VCFs at input will be
                        merged into the same RDF file. If not, an RDF per VCF
                        will be created
  --output_folder OUTPUT_FOLDER, -d OUTPUT_FOLDER
                        Output folder for the RDF files.
  --threads THREADS, -t THREADS
                        Output RDF file.
```

Examples of use:

  - `python vcf2rdf -i vcf_folder/*.vcf -o graph.ttl -t4`: Process all vcf files in the folder `vcf_folder` and store the resulting RDF graph into `graph.ttl`. Use 4 threads.
  - `python vcf2rdf -i file.vcf -o graph.ttl -t4`: Process the file `file.vcf` and store the resulting RDF graph into `graph.ttl`. Use 4 threads.
  - `python vcf2rdf -i file1.vcf file2.vcf -o graph.ttl -t4`: Process the files `file1.vcf` and `file2.vcf` and store the resulting RDF graph into `graph.ttl`. Use 4 threads.
  - `python vcf2rdf -i file1.vcf file2.vcf -d variants_ttl -t4`: Process the files `file1.vcf` and `file2.vcf` and store the resulting RDF graphs into `variants_ttl/file1.ttl` and `variants_ttl/file1.ttl`. Use 4 threads.


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

### Get variant annotations

The annotations are stored as IRIs when possible, referring genes, transcripts and alleles:

```
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX gfvo: <https://github.com/BioInterchange/Ontologies/gfvo#>
PREFIX so: <http://purl.obolibrary.org/obo/SO_>
PREFIX fn: <https://github.com/BioInterchange/Ontologies/gfvo_ann_ext#>
select  *
where {
    ?variant a so:0001059 .
    ?variant gfvo:hasAnnotation [
        a fn:FunctionalAnnotation ;
    	so:associated_with [ rdfs:label ?effect ] ;
    	gfvo:refersTo ?allele ;
        gfvo:refersTo ?feature ;
    	gfvo:refersTo ?gene ;
        gfvo:hasAttribute [
            a fn:Impact ;
        	gfvo:hasValue ?impact 
        ] 
    ].
    ?allele a gfvo:SequenceVariant .
    ?feature a gfvo:Feature .
    ?gene a so:0000704 .      
} limit 500
```

Additionally, the annotations are stored through dynamically generated annotation properties to store the exact text of the annotation:

```
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX gfvo: <https://github.com/BioInterchange/Ontologies/gfvo#>
PREFIX so: <http://purl.obolibrary.org/obo/SO_>
PREFIX fn: <https://github.com/BioInterchange/Ontologies/gfvo_ann_ext#>
select  *
where {
    #values ?variant { <https://longseq.com/variant#000eb15a-5f08-5346-5751-143c79b6215f> } .
    ?variant a so:0001059 .
    ?variant gfvo:hasAnnotation [
        a fn:FunctionalAnnotation ;
    	?annotation ?annotationValue
    ].
    FILTER(regex(str(?annotation), str(fn:) ) )
} limit 500
```