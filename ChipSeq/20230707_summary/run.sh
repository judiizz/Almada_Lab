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
