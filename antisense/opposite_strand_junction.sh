# split the gtf file into + and - 
# 7th column in the gtf file has information about which strand a gene is on, so awk selects for + or - and makes a subseted gtf file
awk '$7 ~ /\+/ {print $0}' trimed_gencode.v31.gtf > pos_gencode.v31.gtf
awk '$7 ~ /\-/ {print $0}' trimed_gencode.v31.gtf > neg_gencode.v31.gtf

# make a variable with paths to all the bam files
allfiles=$(ls */*/*.bam)

# make directory to store sam files with the positive genes
mkdir pos_genes

# the code below makes sam files with junction reads mapping to the opposite strand of the gene
# the variable "thisbase" will contain the file name of each of the files being processed
# samtools view -H returns the header of the bed file and writes it out to a sam file in the directory pos_genes
# bedtools intersect looks at which reads in the bed file overlaps with genes in the gtf file 
# from the reads that overlapped with the genes in the gtf file:
## samtools view -q 10 filters for reads that have a MAPQ greater than 10, 
## awk command filters for reads that:
##### 1. have FLAGs where the first alignment in pair is on the negative strand
##### 2. have a "N" in the CIGAR string, as this indicates a junction
##### 3. don't have "A:+" in the XS tag, as this indicates which strand the read mapped to. The reason "A:-" is not used because some reads do not have the XS tag at all, and the reason the command looks at columns 20, 21, and 22 is because the XS tag could be on one of the three columns
for f in $allfiles; do thisbase=$(basename $f | sed "s/\.bam//"); samtools view -H $f > pos_genes/${thisbase}.sam; bedtools intersect -abam $f -b pos_gencode.v31.gtf | samtools view -q 10 - | awk  '$2~/81|161|163|83|89/ && $6~/N/ && ($20 !~ /A:+/ && $21 !~ /A:+/ && $22 !~ /A:+/)' >> pos_genes/${thisbase}.sam; done

# make directory to store sam files with the negative genes
mkdir neg_genes

# same as the for loop above
for f in $allfiles; do thisbase=$(basename $f | sed "s/\.bam//"); samtools view -H $f > neg_genes/${thisbase}.sam; bedtools intersect -abam $f -b neg_gencode.v31.gtf | samtools view -q 10 - | awk  '$2!~/81|161|163|83|89/ && $6~/N/ && ($20 !~ /A:-/ && $21 !~ /A:-/ && $22 !~ /A:-/)' >> neg_genes/${thisbase}.sam; done
