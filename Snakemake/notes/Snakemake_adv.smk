
# STEP 1: SPECIFYING THE NUMBER OF USED THREADS
# use more than 1 thread could speed up the computation
# could specify the threads a rule needs
# by default, the rule is assumed to need 1 thread
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    
    # add threads: 
    threads: 8

    # propagated to the shell command
    shell:
        "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"

# execution: the number of threads the jobs need is considered by the Snakemake scheduler
# the sum of the threads of all jobs running at the same time DOES NOT EXCEED a given number of CPU cores
# when core > threads: snakemake will saturate the remaining cores with other jobs
# when core < threads: number of threads a rule uses will be reduced to the number of given cores
# if --core without number, all available cores are used




# STEP 2: CONFIG FILES
# customize your workflow: congfig file mechanism: JSON or YAML

# add this to the top of the snakefile
configfile: "config.yaml"
# snakemake load the config file and store the contents into a globally available dictionary called config
# create a new file called: config.yaml
# add this into the config.yaml
samples:
    A: data/samples/A.fastq
    B: data/samples/B.fastq

# could now remove the statement define SAMPLES and change here
rule bcftools_call:

    # propagate every element under samples
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
    output:
        "calls/all.vcf"
    shell:
        "bcftools mpileup -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"