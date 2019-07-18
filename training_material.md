Tapis app structure

```sh
package-name-version.dot.dot
|--system_name
|----test.sh
|----script.template
|----app.json
```

*Prerequesites 

1. Completed `tacc-systems-create`
2. Completed `clients-create`
3. Completed `auth-tokens-create`
4. Installed `agave CLI`

*Tapisapp development proceeds via the following steps:*

1. Test singularity module on the executionSystem
2. Describe the application using an Agave app description 
3. Create a shell template for running the app
4. Upload the application directory to a storageSystem
5. Post the app description to the Agave apps service
6. Debug your app by running jobs and updating the app until it works as intended

1. Test singularity module on the executionSystem

 *Download test file

 ```sh
 $files-get -S data.iplantcollaborative.org /shared/iplantcollaborative/example_data/Samtools_mpileup/ex1.bam
 ```

 *Set up $PATH using singularity images on biocontainers repository

 ```
 $module load biocontainers
 $module load samtools/ctr-1.9--h91753b0_5
 ```

 *Test out commands on the command line prompt
 ```
 $samtools sort -o ex1_sort.bam -@ 12 ex1.bam
 ```

2. Describe the application using an Agave app description (i.e. json)

```
{"available": true,
 "checkpointable": false,
 "defaultMemoryPerNode": 64,
 "defaultProcessorsPerNode": 12,
 "defaultMaxRunTime": "02:00:00",
 "defaultNodeCount": 1,
 "defaultQueue": "normal",
 "deploymentPath":"jawon/applications/samtools-1.9/stampede2",
 "templatePath":"test.template",
 "testPath":"test.sh ",
 "deploymentSystem": "data.iplantcollaborative.org",
 "executionSystem": "tacc-stampede2-wonaya",
 "executionType": "HPC",
 "helpURI": "http://samtools.github.io",
 "label": "samtools_sort",
 "longDescription": "",
 "name": "jawon-samsort-s2",
 "parallelism": "SERIAL",
 "shortDescription": "sort sam/bam file",
 "templatePath": "test.template",
 "testPath": "test.sh",
 "version": "1.9",
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
    {"id":"outname",
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
        {"default":500000000,
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

3. Create a shell template for running the app

```
#input 
inputBam=${inputBam}
#param
outname=${outname}
maxSortMemory=${maxSortMemory}
sortByName=${sortByName}

module purge
module load TACC
module load biocontainers
module load samtools/ctr-1.9--h91753b0_5

ARGS=""
ARGS=" -o ${outname}  "

# Run the actual program
echo "samtools sort -@ 12 ${ARGS} ${inputBam}"
```

4. Upload the application directory to a storageSystem

```
$imkdir -p /iplant/home/jawon/applications/samtools-1.9/stampede2/
$iput -r * /iplant/home/jawon/applications/samtools-1.9/stampede2/.
```

5. Post the app description to the Agave apps service

```
$apps-addupdate -F samsort.json
```

6. Debug your app by running jobs and updating the app until it works as intended

```
{
  "name": "jawon-samsort-s2 test-1563469468",
  "appId": "jawon-samsort-s2-1.9",
  "archive": true,
  "inputs": {
    "inputBam": "agave://data.iplantcollaborative.org/shared/iplantcollaborative/example_data/Samtools_sort_BAM/sample_1.bam"
  },
  "parameters": {
    "maxSortMemory": 500000000
  }
}
```

Lastly, 

```
$jobs-submit -F job.json
```

You can monitor the status of your job

```
$jobs-status $JOB-ID
```







![Alt Text](https://gifimage.net/wp-content/uploads/2018/06/what-have-you-done-gif-1.gif)
