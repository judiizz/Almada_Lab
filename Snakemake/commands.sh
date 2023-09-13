
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

Judyzdf20030516#

# step 6: calling genomic variants
snakemake -s Snakefile_2.smk --core 1 plots/quals.svg


# step 7: adding a target rule
# report shows the target rule triggered the job plot_quals
# report shows two jobs: plot_quals and all, plot_quals is triggered by all
snakemake -s Snakefile_2.smk -n



# EXERCISE

# Create the DAG for the complete workflow
snakemake -s Snakefile_2.smk --dag plots/quals.svg | dot -Tsvg > dag_comp.svg


# execute the complete workflow and have a look at the resulting plot/quals.svg
snakemake -s Snakefile_2.smk --core 1 plots/quals.svg


# Snakemake provides handy flags for forcing re-execution of parts of the workflow
# Have a look at the command line help with snakemake --help and search for the flag --forcerun
    #  Force the re-execution or creation of the given rules or files. Use this option
    #  if you changed a rule and want to have all its output in your workflow updated.
# use this flag to re-execute the rule samtools_sort and see what happens
    # report shows forced execution

snakemake -s Snakefile_2.smk --forcerun sorted_reads/B.bam  --core 1 sorted_reads/B.bam 



# Snakemake display the reason for each job (under reason:)
# perform a dry-run -np that forces some rules to be reexecuted
# report shows forced execution for rules bwa_map and
snakemake -s Snakefile_2.smk -np mapped_reads/{A,B}.bam calls/all.vcf --forcerun mapped_reads/{A,B}.bam calls/all.vcf





################################################
# ADVANCED features

# STEP 1: specifying the threads
snakemake --cores 10

# use --forceall to enforce a complete re-execution of the workflow
    # execution of the selected (or the first) rule and all rules it is dependent on regardless of already created output.
# coombine with different values for --cores
    # the threads = 8 is specified in the rule bwa_map
snakemake -s Snakemake_adv.smk --core 10 --forceall plots/quals.svg 
# bwa_map: 0.504, no B


# more cores
snakemake -s Snakemake_adv.smk --core 50 --forceall plots/quals.svg 
# bwa_map: A and B 0.864

# more cores
snakemake -s Snakemake_adv.smk --core 100 --forceall plots/quals.svg 
# bwa_map: A and B 0.764

# result: adding more cores than threads could perform more tasks?



# STEP 2: config files
# functioned well
snakemake -s Snakemake_adv.smk --core 50 --forceall plots/quals.svg 