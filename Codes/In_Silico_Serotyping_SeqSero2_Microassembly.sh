#!/bin/bash
# sh In_Silico_Serotyping_SeqSero2_Microassembly.sh 
# run from your directory with fasta files; deposits results in SeqSero2_Output directory
# the script only works for paired-end reads
# Ruixi Chen rc836@cornell.edu
# Oct 26, 2021


mkdir SeqSero2_Output

for i in *_1.fastq.gz #change the extension accordingly; for fna.gz files, change "*_1.fastq.gz" to "*_1.fna.gz"
do
SAMPLE=$(echo ${i} | sed "s/_1.fastq\.gz//") 
SeqSero2_package.py -p 10 -t 2 -i ${SAMPLE}_1.fastq.gz ${SAMPLE}_2.fastq.gz
done
