library(readr)

# read in files
setwd("input_files/")
pheno.data <- read_tsv("name_conversion.txt")
pheno.data <- as.data.frame(pheno.data)

junction_table <- read_tsv("all_junction_table.txt")

# change col names to match ones in name conversion table
colnames(junction_table) <- gsub("^121317.", "", colnames(junction_table))
colnames(junction_table) <- gsub("_", ".", colnames(junction_table))
colnames(junction_table) <- gsub("-", ".", colnames(junction_table))

# change col names to match rest of paper
matched <- match(colnames(junction_table[,6:ncol(junction_table)]), pheno.data$index_sequence)

colnames(junction_table) <- c(colnames(junction_table[,1:5]), pheno.data$individual_number[matched])

# write the final table into a file
write.table(junction_table, 
            "all_junction_table_converted.txt", 
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
