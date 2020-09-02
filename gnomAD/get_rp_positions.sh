#!/bin/bash


# extracted RP positions
cat trimed_gencode.v31.gtf | awk '$3 ~ /gene/ && $14 ~ /^"RP[SL][0-9]+"/ {print $14, $1, $4, $5}' | sed 's/"//g'| sed 's/;//g' | sed 's/ /,/g' > RP_positions.csv

# extract RP exon positions
echo "gene_name,chr,start,end,transcript_name,exon_number" > RP_exon_positions.csv; cat trimed_gencode.v31.gtf | awk '$3 ~ /exon/ && $16 ~ /^"RP[SL][0-9P]+A?"/ {print $16, $1, $4, $5, $20, $22}' | sed 's/"//g'| sed 's/;//g' | sed 's/ /,/g' >> RP_exon_positions.csv


