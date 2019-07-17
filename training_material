Tapis app structure

```sh
package-name-version.dot.dot
|--system_name
|----test.sh
|----script.template
|----app.json
```

*Tapisapp development proceeds via the following steps:*

1. Test singularity module on the executionSystem
2. Describe the application using an Agave app description 
3. Create a shell template for running the app
4. Upload the application directory to a storageSystem
5. Post the app description to the Agave apps service
6. Debug your app by running jobs and updating the app until it works as intended

1. Test singularity module on the executionSystem

Download test file

```sh
files-get -S data.iplantcollaborative.org /shared/iplantcollaborative/example_data/Samtools_mpileup/ex1.bam
```

```
module load biocontainers
module load samtools/ctr-1.9--h91753b0_5
```

```
samtools sort -o ex1_sort.bam -@ 12 ex1.bam
```

2. Describe the application using an Agave app description (i.e. json)

```
{"name":"jawon-samtools-sort",
 "parallelism":"SERIAL",
 "version":"1.9",
 "helpURI":"http://samtools.sourceforge.net/samtools.shtml",
 "label":"",
 "shortDescription":"Sort SAM file",
 "longDescription":"Perform sort operation on a SAM file",
 "author":"Jawon Song",
 "datePublished":"07/17/2019",
 "tags":["aligner","NGS","SAM"],
 "ontology":["http://sswapmeet.sswap.info/agave/apps/Application", "http://sswapmeet.sswap.info/sequenceServices/SequenceServices"],
 "executionHost":"lonestar4.tacc.teragrid.org",
 "executionType":"HPC",
 "deploymentPath":"/iplant/home/jawon/applications/samtools-0.1.8/lonestar",
 "templatePath":"samtools-sort.sge",
 "testPath":"test-sort.sh ",
 "checkpointable":"false",
 "modules":["purge","load TACC","load iRODS"],
 "inputs":[
    {"id":"inputBam",
     "value":
        {"default":"",
         "validator":"",
         "visible":true,
         "required":true},
     "details":
        {"label":"The BAM file to sort",
         "description":""},
     "semantics":
        {"ontology":["http://sswapmeet.sswap.info/mime/application/X-bam"],
         "minCardinality":1,
         "maxCardinality":1,
         "fileTypes":["raw-0"]}}],
 "parameters":[
    {"id":"outPrefix",
     "value":
        {"default":"sorted",
         "type":"string",
         "validator":"",
         "required":false,
         "visible":true},
     "details":
        {"label":"Prefix for the output file",
         "visible":true},
     "semantics":
        {"ontology":["xs:string"]}},
    {"id":"maxSortMemory",
     "value":
        {"default":"500000000",
         "type":"number",
         "validator":"",
         "required":true,
         "visible":true},
     "details":
        {"label":"Maxiumum memory in bytes, used for sorting",
         "visible":true},
     "semantics":
        {"ontology":["xs:integer"]}},
    {"id":"sortByName",
     "value":
        {"default":false,
         "type":"bool",
         "validator":"",
         "required":false,
         "visible":true},
     "details":
        {"label":"Sort by name rather than coordinate",
         "visible":true},
     "semantics":
        {"ontology":["xs:boolean"]}}],
"outputs":[{"id":"bam",
     "value":
        {"default":"",
         "validator":"",
         "required":true,
         "visible":true},
     "details":
        {"label":"BAM file",
         "description":""},
     "semantics":
        {"ontology":["http://sswapmeet.sswap.info/mime/application/X-bam"],
         "minCardinality":1,
         "fileTypes":["raw-0"]}}]}
```
