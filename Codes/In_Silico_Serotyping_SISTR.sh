#!/bin/bash
# sh In_Silico_Serotyping_SISTR.sh 
# run from your directory with fasta files; deposits results in sistr_out directory
# Laura Carroll lmc297@cornell.edu
# Jan 25, 2019

mkdir sistr_out

for f in *.fna #change the extension accordingly; for FASTA files, change "*.fna" to "*.fasta"
do
sistr -i ./$f $f -o sistr_out/${f%.fas} -f csv --qc -t 1 --verbose 
done
