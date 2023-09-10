

# STEP 1: MAPPING READS
# snakemake rule name
rule bwa_map:

    # used by the rule
    input:
        "data/genome.fa",
        "data/samples/A.fastq"
    
    # created by the rule: in this case, bwa and samtools creates a .bam file
    output:
        "mapped_reads/A.bam"
    
    # shell commands to execute
    shell:
        # snakemake automatically replace {input} woth the input files before executing the command
        # concatenate them with whitespace: data/genome.fa data/samples/A.fastq
        # same as {output}  

        # baw mem pipes into samtools to generate compressed BAM file
        # output of samtools is redirected to the output file: > {output}
        "bwa mem {input} | samtools view -Sb - > {output}"


# STEP 2: GENERALIZING THE READ MAPPING RULE
