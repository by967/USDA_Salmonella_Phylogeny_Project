#!/bin/bash
# sh Quality_Control_QUAST.sh
# QUAST - Quality Assessment Tool for Genome Assemblies
# jk2739@cornell.edu
# all assemblies should be stored within the working directory

mkdir quast_results

for f in *.fna #change the extension accordingly; for FASTA files, change "*.fna" to "*.fasta"
do
python quast.py -o ./quast_results/quast_${f%_contigs.fna} --min-contig 1 $f
done

#collect report txt files

mkdir quast_reports
for f in *.fna
do cd quast_results/quast_${f%_contigs.fna}
cat report.txt > ${f%_contigs.fna}_report.txt
cp ${f%_contigs.fna}_report.txt ../../quast_reports
cd ../..
done


