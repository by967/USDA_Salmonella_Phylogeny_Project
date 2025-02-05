#!/bin/bash
# sh In_Silico_Serotyping_SeqSero2_kmer.sh 
# run from your directory with fna files; deposits results in SeqSero2_Output directory
# Ruixi Chen rc836@cornell.edu
# Oct 26, 2021


mkdir SeqSero2_Output

for f in *.fna #change the extension accordingly; for FASTA files, change "*.fna" to "*.fasta"
do
SeqSero2_package.py -m k -t 4 -i ./$f -d SeqSero2_Output/${f%.fna}
done
