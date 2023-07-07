---
title: "Week_6.26"
output:
  pdf_document: default
  html_document: default
---
****************************************
# 2023.6.28
## Goal: move the "practice" from scratch1 to aealmada_561 
*operated on Discovery cluster to avoid data loss*

```{bash}
# find the practice folder: /scratch1/difeizhu/practice
cd
cds
ls
cd practice
pwd

# move practice folder to aealmada_561 directory
mv /scratch1/difeizhu/practice /project/aealmada_561
cd /project/aealmada_561
ls
cd practice
pwd

# move practice folder to judyz directory
mv /project/aealmada_561/practice /project/aealmada_561/judyz
```



****************************************
# 2023.6.29
## Goal: debugging nextflow ChIPseq pipeline by updating the nf_core version

Errors: ERROR ~ No such variable: ch_design_reads_csv

* error occured in the main.nf file
* guessed could be the nf_core version problem

### Find the files to run ChIPseq pipeline
```{bash}
cd /project/aealmada_561
ls

cd judyz
ls

cd practice
ls
```

### Deleted the output files from last run
*same procedure before submitting the new job*

* deleted from GUI
* deleted directories: work, results
* deleted files: nextflow.err, nextflow.out
* checked status by command: ls


### New job 1
* modified run.sh only

By updating the nf_core version, several additional parameters are required
[https://nf-co.re/launch?id=1687564431_64f34655737a#input_output_options]
```{bash}
nano run.sh

# file displayed:

#!/bin/bash
#SBATCH --account=aealmada_561
#SBATCH --partition=main
#SBATCH --job-name=PolII_test
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=difeizhu@usc.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=6gb
#SBATCH --time=14:00:00
#SBATCH --output=nextflow.out
#SBATCH --error=nextflow.err
pwd; hostname; date
echo "Lets go"
nextflow run nf-core/chipseq -r 1.2.2 \  # modify this line
-profile singularity \
--single_end \
--input design.csv \
# need to add --outdir parameters
# need to add --genome parameters
--fasta /scratch1/difeizhu/practice/gencode/GRCh38.p13.genome.fa \
--gtf /scratch1/difeizhu/practice/gencode/gencode.v32.annotation.gtf \
--macs_gsize 3.2e9 \
--blacklist /scratch1/difeizhu/practice/gencode/hg38-blacklist.v2.bed \
--email difeizhu@usc.edu \
-resume \
-c nextflow.config
date
```


Modified version of run.sh

* added --outdir parameters
* added --genome parameters
```{bash}
nano run.sh

#!/bin/bash
#SBATCH --account=aealmada_561
#SBATCH --partition=main
#SBATCH --job-name=PolII_test
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=difeizhu@usc.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=6gb
#SBATCH --time=14:00:00
#SBATCH --output=nextflow.out
#SBATCH --error=nextflow.err
pwd; hostname; date
echo "Lets go"
nextflow run nf-core/chipseq -r 2.0.0 \
-profile singularity \
--single_end \
--input design.csv \
--outdir /project/aealmada_561/judyz/practice \
--genome GRCh38 \
--fasta /scratch1/difeizhu/practice/gencode/GRCh38.p13.genome.fa \
--gtf /scratch1/difeizhu/practice/gencode/gencode.v32.annotation.gtf \
--macs_gsize 3.2e9 \
--blacklist /scratch1/difeizhu/practice/gencode/hg38-blacklist.v2.bed \
--email difeizhu@usc.edu \
-resume \
-c nextflow.config
date
```


Submitted job 1
```{bash}
sbatch run.sh
```

**Results**: failed :(

* new folder generated: pipeline_info
* "results" folder disappeared
* checked nextflow.out file
  + still launched the original version 1.2.2



### New job 2
*might need to call the 2.0.0 version first
```{bash}
nextflow pull -r 2.0.0 nf-core/chipseq
## done - revision: 51eba00b32 [2.0.0]
```


Deleted files and submitted job 2
```{bash}
sbatch run.sh
```


**Results**: failed :(

* checked nextflow.out file

  + launched the new version
  + ERROR ~ * --macs_gsize: expected type: Number, found: String (3.2e9)
  + WARN: Found unexpected parameters:* --single_end: true - Ignore this warning: params.schema_ignore_params = "single_end" 


Checked .nextflow.log file for more details
```{bash}
cat .nextflow.log
## basically same as the summary
```



### New job 3
*for macs_gsize, convert 3.2e9 to numeric version
*ignore the warning first

Modified version of run.sh

* --macs_gsize 3200000000 \
```{bash}
nano run.sh

#!/bin/bash
#SBATCH --account=aealmada_561
#SBATCH --partition=main
#SBATCH --job-name=PolII_test
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=difeizhu@usc.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=6gb
#SBATCH --time=14:00:00
#SBATCH --output=nextflow.out
#SBATCH --error=nextflow.err
pwd; hostname; date
echo "Lets go"
nextflow run nf-core/chipseq -r 2.0.0 \
-profile singularity \
--single_end \
--input design.csv \
--outdir /project/aealmada_561/judyz/practice \
--genome GRCh38 \
--fasta /scratch1/difeizhu/practice/gencode/GRCh38.p13.genome.fa \
--gtf /scratch1/difeizhu/practice/gencode/gencode.v32.annotation.gtf \
--macs_gsize 3200000000 \
--blacklist /scratch1/difeizhu/practice/gencode/hg38-blacklist.v2.bed \
--email difeizhu@usc.edu \
-resume \
-c nextflow.config
date
```


Deleted files and submitted job 3
```{bash}
sbatch run.sh
```


**Results**: failed :(

* checked nextflow.out file

  + ERROR ~ Unable to read script: '/home1/difeizhu/.nextflow/assets/nf-core/chipseq/./workflows/chipseq.nf' -- cause: /scratch1/difeizhu/practice/gencode/GRCh38.p13.genome.fa
*DARNNN i forgot to modify the file paths since i moved the folders yesterday!!*



### New job 4
* modified all the file paths

Modified version of run.sh

* changed file paths
```{bash}
#!/bin/bash
#SBATCH --account=aealmada_561
#SBATCH --partition=main
#SBATCH --job-name=PolII_test
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=difeizhu@usc.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=6gb
#SBATCH --time=14:00:00
#SBATCH --output=nextflow.out
#SBATCH --error=nextflow.err
pwd; hostname; date
echo "Lets go"
nextflow run nf-core/chipseq -r 2.0.0 \
-profile singularity \
--single_end \
--input design.csv \
--outdir /project/aealmada_561/judyz/practice \
--genome GRCh38 \
--fasta /project/aealmada_561/judyz/practice/gencode/GRCh38.p13.genome.fa \
--gtf /project/aealmada_561/judyz/practice/gencode/gencode.v32.annotation.gtf \
--macs_gsize 3200000000 \
--blacklist /project/aealmada_561/judyz/practice/gencode/hg38-blacklist.v2.bed \
--email difeizhu@usc.edu \
-resume \
-c nextflow.config
date
```


Original version of design.csv
```{bash}
group,replicate,fastq_1,fastq_2,antibody,control
POLR2A,3,/scratch1/difeizhu/practice/fastq/ENCFF000XWJ.fastq.gz,,POLR2A,ENCSR000EEN
POLR2A,2,/scratch1/difeizhu/practice/fastq/ENCFF000PNT.fastq.gz,,POLR2A,ENCSR000BLH
POLR2A,4,/scratch1/difeizhu/practice/fastq/ENCFF000XWR.fastq.gz,,POLR2A,ENCSR000EEN
POLR2A,1,/scratch1/difeizhu/practice/fastq/ENCFF000PNM.fastq.gz,,POLR2A,ENCSR000BLH
ENCSR000EEN,1,/scratch1/difeizhu/practice/fastq/ENCFF000XTF.fastq.gz,,,
ENCSR000BLH,1,/scratch1/difeizhu/practice/fastq/ENCFF000POQ.fastq.gz,,,
```

File path:
/project/aealmada_561/judyz/practice/fastq

Modified version of design.csv
```{bash}
group,replicate,fastq_1,fastq_2,antibody,control
POLR2A,3,/project/aealmada_561/judyz/practice/fastq/ENCFF000XWJ.fastq.gz,,POLR2A,ENCSR000EEN
POLR2A,2,/project/aealmada_561/judyz/practice/fastq/ENCFF000PNT.fastq.gz,,POLR2A,ENCSR000BLH
POLR2A,4,/project/aealmada_561/judyz/practice/fastq/ENCFF000XWR.fastq.gz,,POLR2A,ENCSR000EEN
POLR2A,1,/project/aealmada_561/judyz/practice/fastq/ENCFF000PNM.fastq.gz,,POLR2A,ENCSR000BLH
ENCSR000EEN,1,/project/aealmada_561/judyz/practice/fastq/ENCFF000XTF.fastq.gz,,,
ENCSR000BLH,1,/project/aealmada_561/judyz/practice/fastq/ENCFF000POQ.fastq.gz,,,
```


Deleted files and submitted job 4
```{bash}
sbatch run.sh
```


**Results**: failed :(

* observation: longer runtime, made some progress??
* new folder created: genome
* checked nextflow.out file *very long*
  + ERROR ~ Error executing process > 'NFCORE_CHIPSEQ:CHIPSEQ:PREPARE_GENOME:GTF2BED (gencode.v32.annotation.gtf)' Caused by:
  Failed to submit process to grid scheduler for execution


**Analysis of the errors**
ERROR ~ Error executing process > 'NFCORE_CHIPSEQ:CHIPSEQ:PREPARE_GENOME:GTF2BED (gencode.v32.annotation.gtf)'
*got the genome folder, maybe the fasta worked!!!*

Caused by: Failed to submit process to grid scheduler for execution
*googled this one*

Command output: sbatch: error: invalid partition specified: short,  sbatch: error: Batch job submission failed: Invalid partition name specified
*maybe change partition name to short, currently is main*


*What I've learned*
In the nextflow.config file
```{bash}
nextflow.enable.dsl=2
process {
  executor='slurm'
  queue='short' # same as setting the partition to short
  memory='32 GB'
  maxForks=10
}
```



### New job 5
* modified nextflow.config file

Modified version of nextflow.config

* deleted queue='short' because error showed invalid partition specified: short
```{bash}
nextflow.enable.dsl=2
process {
  executor='slurm'
  memory='32 GB'
  maxForks=10
}
```
*met weird error messages: disk quota exceeded. Didn't do anything, error disappeared.*


Deleted files and submitted job 5
```{bash}
sbatch run.sh
```


**Results**: failed :(

* observation: longer runtime, made more progress
* previous error fixed!!
* checked nextflow.out file
  + ERROR: Please check samplesheet header -> group,replicate,fastq_1,fastq_2,antibody,control != sample,fastq_1,fastq_2,antibody,control




****************************************
# 2023.6.30
## Goal: keep debugging chipseq pipeline. Modify design.csv file


### New job 1
Current version of design.csv
```{bash}
group,replicate,fastq_1,fastq_2,antibody,control
POLR2A,3,/project/aealmada_561/judyz/practice/fastq/ENCFF000XWJ.fastq.gz,,POLR2A,ENCSR000EEN
POLR2A,2,/project/aealmada_561/judyz/practice/fastq/ENCFF000PNT.fastq.gz,,POLR2A,ENCSR000BLH
POLR2A,4,/project/aealmada_561/judyz/practice/fastq/ENCFF000XWR.fastq.gz,,POLR2A,ENCSR000EEN
POLR2A,1,/project/aealmada_561/judyz/practice/fastq/ENCFF000PNM.fastq.gz,,POLR2A,ENCSR000BLH
ENCSR000EEN,1,/project/aealmada_561/judyz/practice/fastq/ENCFF000XTF.fastq.gz,,,
ENCSR000BLH,1,/project/aealmada_561/judyz/practice/fastq/ENCFF000POQ.fastq.gz,,,
```


Modified version of design.csv
Followed this template [https://nf-co.re/chipseq/2.0.0/usage]
*The group and replicate columns were replaced with a single sample column as of v2.0 of the pipeline. The sample column is essentially a concatenation of the group and replicate columns.*
```{bash}
sample,fastq_1,fastq_2,antibody,control
POLR2A_REP3,/project/aealmada_561/judyz/practice/fastq/ENCFF000XWJ.fastq.gz,,POLR2A,ENCSR000EEN
POLR2A_REP2,/project/aealmada_561/judyz/practice/fastq/ENCFF000PNT.fastq.gz,,POLR2A,ENCSR000BLH
POLR2A_REP4,/project/aealmada_561/judyz/practice/fastq/ENCFF000XWR.fastq.gz,,POLR2A,ENCSR000EEN
POLR2A_REP1,/project/aealmada_561/judyz/practice/fastq/ENCFF000PNM.fastq.gz,,POLR2A,ENCSR000BLH
ENCSR000EEN_REP1,/project/aealmada_561/judyz/practice/fastq/ENCFF000XTF.fastq.gz,,,
ENCSR000BLH_REP1,/project/aealmada_561/judyz/practice/fastq/ENCFF000POQ.fastq.gz,,,
```

While modifying the file: Error writing design.csv: Disk quota exceeded 
*Don't understand the error, it also appeared yesterday.*

* i restarted the cluster terminal several time and the problem disappeared
* future strategy: try not to modify a lot in the file


Deleted files and submitted job 1
```{bash}
sbatch run.sh
```


**Results**: failed :(

* checked nextflow.out file *very long*
  + ERROR: Please check samplesheet -> Control identifier has to match does a provided sample identifier!
  + Control: 'ENCSR000BLH'

*need to understand the parameters*




