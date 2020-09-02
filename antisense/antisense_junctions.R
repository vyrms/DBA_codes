library(Rsubread)
library(dplyr)
library(tibble)
library(biomaRt)


# move the positive/ negative gtf files to the respective pos_genes/ neg_genes directory

# go to the directory pos_genes
setwd('pos_genes/')
# list the sam files
pos_files <- list.files(pattern =".sam$")

# makes a list from featureCounts
pos_genes <- featureCounts(pos_files, annot.ext = 'pos_gencode.v31.gtf', 
                            isGTFAnnotationFile = TRUE,
                            isPairedEnd=TRUE)
# $counts contains the number of reads per gene
pos_counts <- data.frame(pos_genes$counts)

#gene name conversions
ids <- rownames(pos_counts) # gene names from counts
ids <- gsub("\\.\\d*$", "", ids) # remove version tag
ensembl_database <- useEnsembl(biomart="ensembl", 
                               dataset="hsapiens_gene_ensembl") # connect to ensemble database for gene names
results <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 filter = "ensembl_gene_id", 
                 values = ids, mart = ensembl_database) # use genes to obtain hgnc_symbol

#for each entry in results find corresponding tag
for(i in 1:nrow(results)){
  
  idx <- grep(paste(results$ensembl_gene_id[i], 
                    "\\.\\d*$", sep= ""), 
              rownames(pos_counts)) # find corresponding tag
  
  if (nchar(results$hgnc_symbol[i]) == 0){ # if gene result is bank replace with id
    rownames(pos_counts)[idx] <- results$ensembl_gene_id[i]}
  else if (results$hgnc_symbol[i] %in% rownames(pos_counts)){ # if symbol already used replaced with id
    rownames(pos_counts)[idx] <- results$ensembl_gene_id[i] 
  }else{ # else replace with gene symbol
    rownames(pos_counts)[idx] <- results$hgnc_symbol[i]
  }
}


# filter for genes that have more than n reads in m samples
n = 5 # number of reads in a sample
m = 2 # number of samples that have more than n reads
i = 1
pos_filtered <- pos_counts
while(i <= nrow(pos_filtered)){
  keep <- c()
  
  pass <- sum(pos_filtered[i,-1] > n)
  if(pass < m){
    pos_filtered <- pos_filtered[-i,]
    i = i - 1
  }
  i = i + 1
}







# negative strand genes

setwd('D:/hosei_to_be_deleted/neg_genes')
neg_files <- list.files(pattern =".sam$")

neg_genes <- featureCounts(neg_files, annot.ext = 'neg_gencode.v31.gtf', 
                           isGTFAnnotationFile = TRUE,
                           isPairedEnd=TRUE)
neg_counts <- data.frame(neg_genes$counts)

#gene name conversions
ids <- rownames(neg_counts) # gene names from counts
ids <- gsub("\\.\\d*$", "", ids) # remove version tag
ensembl_database <- useEnsembl(biomart="ensembl", 
                               dataset="hsapiens_gene_ensembl") # connect to ensemble database for gene names
results <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 filter = "ensembl_gene_id", 
                 values = ids, mart = ensembl_database) # use genes to obtain hgnc_symbol

#for each entry in results find corresponding tag
for(i in 1:nrow(results)){
  
  idx <- grep(paste(results$ensembl_gene_id[i], 
                    "\\.\\d*$", sep= ""), 
              rownames(neg_counts)) # find corresponding tag
  
  if (nchar(results$hgnc_symbol[i]) == 0){ # if gene result is bank replace with id
    rownames(neg_counts)[idx] <- results$ensembl_gene_id[i]}
  else if (results$hgnc_symbol[i] %in% rownames(neg_counts)){ # if symbol already used replaced with id
    rownames(neg_counts)[idx] <- results$ensembl_gene_id[i] 
  }else{ # else replace with gene symbol
    rownames(neg_counts)[idx] <- results$hgnc_symbol[i]
  }
}

# filter for genes that have more than n reads in m samples (n = 8, m = 2)
i = 1
neg_filtered <- neg_counts
while(i <= nrow(neg_filtered)){
  keep <- c()
  
  pass <- sum(neg_filtered[i,-1] > n)
  if(pass < m){
    neg_filtered <- neg_filtered[-i,]
    i = i - 1
  }
  i = i + 1
}


# switch wd to where files should be saved
setwd('D:/hosei_to_be_deleted/opposite')


# export
write.table(pos_filtered, 
            "opposite_pos_genes_5_2.csv", 
            sep=",", quote = FALSE, row.names = TRUE, col.names = TRUE) 
write.table(neg_filtered, 
            "opposite_neg_genes_5_2.csv", 
            sep=",", quote = FALSE, row.names = TRUE, col.names = TRUE) 



# row sums
rowSums(pos_filtered)


