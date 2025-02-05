#!/bin/bash
# sh De_Novo_Assembly_SKESA.sh
# use SKESA to assemble raw reads
# run from your directory with fastq.gz files
# ly276@cornell.edu 
# 08/29/2023

cd "$1" || exit 1 # navigate to your working directory

# function to run skesa for a pair of files
run_skesa() {
  local file="$1"
  skesa --reads "$file","${file%_1.fastq.gz}_2.fastq.gz" --cores 20 --memory 64 --use_paired_ends --contigs_out "${file%_1.fastq.gz}_skesa.fasta"
} #change the extension accordingly; for fna.gz files, change "*_1.fastq.gz" to "*_1.fna.gz"

# export the function so that parallel can access it
export -f run_skesa

# run skesa in parallel for up to 10 jobs
find . -maxdepth 1 -name "*_1.fastq.gz" | parallel -j 10 run_skesa
