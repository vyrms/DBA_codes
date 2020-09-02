library(dplyr)

# import lof metrics data
test <- read.table("D:/input_files/gnomad.v2.1.1.lof_metrics.by_gene.txt", sep="\t", header=TRUE)

# select RPs
test <- test[grepl("^RP[LS][0-9P]+A?$", test$gene),]

# add new column for lof difference
test <- test %>% mutate(oe_lof_diff = oe_lof_upper - oe_lof_lower)

# filter out genes with lof difference greater than equal to 1
high_oe_lof_variance <- test %>% filter(oe_lof_diff >= 1)
