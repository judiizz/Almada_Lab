
# step 1: mapping reads
# file: Snakefile_1a.smk

# specify .smk file by option -s
# only showed execution plan instead of actually performing the steps.
snakemake -s Snakefile_1.smk -np mapped_reads/A.bam 

# Specify the number of cores used.
snakemake -s Snakefile_1.smk --core 1 mapped_reads/A.bam



# step 2: generalizing 
# generate the target file by replacing the wildcard {sample} with the value B
snakemake -s Snakefile_2.smk -np mapped_reads/B.bam

# specify multiple targets
snakemake -s Snakefile_2.smk -np mapped_reads/A.bam mapped_reads/B.bam

# Bash magic: put multiple targets in a single pass {A,B}
# creating the given path for each element and separating the resulting paths by a whitespace
# not snakemake property
snakemake -s Snakefile_2.smk -np mapped_reads/{A,B}.bam

# try to execute by specifying the core
# only create the output file mapped_reads/B.bam, because already executed the workflow before
# no input file is newer than the output file mapped_reads/A.bam
snakemake -s Snakefile_2.smk --core 1 mapped_reads/{A,B}.bam

# update the file modification date of the input file data/samples/A.fastq
touch data/samples/A.fastq

# rerun the job
# report showed wildcard: sample=A
snakemake -s Snakefile_2.smk --core 1 mapped_reads/{A,B}.bam




# step 3: sorting read alignments
# overview of this rule
snakemake -s Snakefile_2.smk -np sorted_reads/B.bam

# try to execute by specifying the core
snakemake -s Snakefile_2.smk --core 1 sorted_reads/B.bam

# try to update the file modification of the first run, and see how the both runs work
touch data/samples/B.fastq

# rerun the job
# report showed both jobs being executed :) 
# the second job will be executed again, because we updated the file from the begining!
snakemake -s Snakefile_2.smk --core 1 sorted_reads/B.bam




# step 4: indexing read alignments and visualizing the DAG of jobs
# create a visualization of the DAG using the dot command by Graphviz
# snakemake specifies the DAG in dot language, pipes into dot, .svg format
snakemake -s Snakefile_2.smk --dag sorted_reads/{A,B}.bam.bai | dot -Tsvg > dag.svg

# create the indexed version A.bam.bai
# report showed for A, rules samtools_sort and samtools_index were executed
# for B, only samtools_index was executed
# because we did not execute the previous run for A, snakemake automatically detected this!!
snakemake -s Snakefile_2.smk --core 1 sorted_reads/{A,B}.bam.bai



# step 5: calling genomic variants
# execute the bcftools_call
snakemake -s Snakefile_2.smk --core 1 calls/all.vcf

# obtain the DAG for the file calls/all.vcf
snakemake -s Snakefile_2.smk --dag calls/all.vcf | dot -Tsvg > dag_2.svg

git remote set-url origin git@github.com:judiizz/Almada_Lab.git
git init

git remote add origin https://github.com/judiizz/Almada_Lab.git

git push origin master