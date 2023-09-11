

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
    # must start with the same directory!!

#*****************************************************************************************************


# STEP 2: GENERALIZING THE READ MAPPING RULE
# by using wildcards: replace the A in "data/samples/A.fastq" with wildcards {sample}

### NOTES:

# When Snakemake determines that this rule can be applied to generate a target file by replacing 
# the wildcard {sample} in the output file with an appropriate value, it will propagate that value 
# to all occurrences of {sample} in the input files and thereby determine the necessary input for the resulting job.

rule bwa_map:
    input:
        "data/genome.fa",
        # wildcard
        "data/samples/{sample}.fastq"
    output:

        # wildcard
        "mapped_reads/{sample}.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"

# all output files of a rule have to contain exactly the SAME wildcards


#*****************************************************************************************************

# STEP 3: SORTING READ ALIGNMENTS
# to sort the read alignments in the BAM files generated above
# samtools sort command
# add this new rule following the baw_map

rule samtools_sort:

    # take input file from mapped_reads directory
    # the input here is the output from the previous rule
    input:
        "mapped_reads/{sample}.bam"

    # store a sorted version in the sorted_reads directory
    output:
        "sorted_reads/{sample}.bam"

    # samtools requires a prefix specified with the flag -T
    # check the samtools commands
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"


#*****************************************************************************************************

# STEP 4: INDEXING READ ALIGNMENTS AND VISUALIZING THE DAG OF JOBS
# use samtools to index the sorted read alignments
# so we can quickly access reads by the genomic location they were mapped to
rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    
    # index bam with format .bam.bai
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"

# DAG: node is each job (rule), edge is the dependency that connect each job



#*****************************************************************************************************

# STEP 5: CALLING GENOMIC VARIANTS
# aggregate the mapped reads from all samples and jointly call genomic variants on them
# combine two denpendencies: samtools and bcftools

# snakemake helper function for collecting input files
expand("sorted_reads/{sample}.bam", sample=SAMPLES)
# SAMPLES is the list
["sorted_reads/A.bam", "sorted_reads/B.bam"]

# function is useful when the pattern contains multiple wildcards: sample name and replicates
expand("sorted_reads/{sample}.{replicate}.bam", sample=SAMPLES, replicate=[0, 1])
# SAMPLES is the list
["sorted_reads/A.0.bam", "sorted_reads/A.1.bam", "sorted_reads/B.0.bam", "sorted_reads/B.1.bam"]

#NOTES:
# first let snakemake kno which samples we want to consider
    # snakemake works backwards from requested output, not from available input

# Python syntax, so we define a list of samples at the top of the snakefile.smk
SAMPLES = ["A", "B"]

# rule
rule bcftools_call:

    # specifying names for input or output files
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
    
    output:
        "calls/all.vcf"

    # with the specified names, the files can be referred here by input.fa or input.bam
    # split the string over multiple indented lines
    shell:
        "bcftools mpileup -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"
