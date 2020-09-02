#!/bin/bash

# run this in a directory with all the BED files generated from IGV export features

# makes a directory to store junctions that start at exon 3 or 4 and end at exon 5
mkdir exon4to5

# get junctionBED that connects to exon 5 (by position)
# ls *.bed lists all the files with the .bed extension in the directory
# parallel can be thought as looping through the bed files
# in awk, $2~/2369[3-4]/ filters for junctions that start in the 23693000s or 23694000s (which includes junctions that start at exon 3 and 4)
# and $3~/23695/ filters for junctions that end in the 23695000s (which only includes junctions that go to exon 5)
ls *.bed | parallel " echo track graphType=junctions > exon4to5/{} ; awk ' \$2~/2369[3-4]/ && \$3~/23695/ ' {} >> exon4to5/{} "

# goes into the directory that was made earlier
cd exon4to5

# saves the .bed file names in the exon4to5 directory into a file
ls *.bed > names.txt

# makes a new directory to store tidy versions of each bed files
mkdir tidy

# make files with tidy version of junction reads
# the loop takes the contents from names.txt and loops  through them
# sed 1d deletes the first line of each file, which is "track graphType=junctions"
# awk -F'[\t,]' sets the field separator as tab and comma
# $2 + $11 and $3 - $12 gives the atual start/end positions of each junction
# $5 gives how many reads mapped to that junction, and $6 gives which strand that junction is on
# then, the contents are sorted by start position, and written to a txt file
while read line; do sed 1d $line | awk -F'[\t,]' ' { print $2 + $11 , "\t" , $3 - $12 , "\t" , $5 , "   \t" , $6} ' | sort > tidy/${line%.*}.txt; done < names.txt
