library(ggplot2)
library(dplyr)

# look at allele frequency per exon
alldata = read.csv("D:/input_files/gnomAD_exon_essential_rp_low_pli_stranded.csv")


## look at splice donor sites (end of exon)
# remove duplicated donor sites
all_donor <- alldata[!duplicated(
  alldata[, !(names(alldata) %in% c("start", 
                                    "start_count",
                                    "start_freq",
                                    "transcript_name", 
                                    "exon_number"))]), ]


## remove RPs with high oe lof variance
# aggregate allele grequency greater than 0.001 to be 0.001
all_donor <- all_donor %>% mutate(end_freq_binned = case_when(end_freq >= 0.001 ~ 0.001,
                                                              TRUE ~ end_freq))
# RP genes related to DBA
dba_causing_gnames <- c("RPS19", "RPL5", "RPS26", "RPL11", "RPL35A", "RPS10", "RPS24", "RPS17", "RPS7",
                        "RPL26", "RPL15", "RPS29", "RPS28", "RPL31", "RPS27", "RPL27", "RPL35", "RPL18", "RPS15A")

# select RPs that don't have high difference in lof metric
rp_donor <- all_donor[grepl("^RP[LS][0-9P]+A?$", all_donor$gene_name) & !(all_donor$gene_name %in% high_oe_lof_variance$gene),]
rp_donor <- rp_donor %>% mutate(dba = case_when(gene_name %in% dba_causing_gnames ~ "DBA",
                                                TRUE ~ "non-DBA"))

# hist of pLI distribution in RPs
rp_donor %>% group_by(gene_name) %>% slice(1) %>% ggplot() +
  geom_histogram(aes(x=pLI)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(title="Histogram of pLI in DBA Related RPs",
       x="pLI",
       y="Count") +
  facet_wrap(. ~ dba)


# plots with 19 DBA related RPs
# make subset for DBA causing RPs
dba_causing <- all_donor[all_donor$gene_name %in% dba_causing_gnames & !(all_donor$gene_name %in% high_oe_lof_variance$gene),]
# make label for DBA causing RPs
dba_causing_label <- paste("DBA Causing RPs (", 
                           length(unique(dba_causing$gene_name)),
                           " genes, ",
                           nrow(dba_causing),
                           " Exons)",
                           sep="")
# make subset for non-DBA causing RPs
not_dba_causing <- all_donor[grepl("^RP[LS][0-9P]+A?$", all_donor$gene_name) & !(all_donor$gene_name %in% dba_causing_gnames) & !(all_donor$gene_name %in% high_oe_lof_variance$gene),]
# make label for non-DBA causing RPs
not_dba_causing_label <- paste("Not DBA Causing RPs (", 
                               length(unique(not_dba_causing$gene_name)),
                               " genes, ",
                               nrow(not_dba_causing),
                               " Exons)",
                               sep="")
# make subset for high pLI genes
high_pli <- all_donor[!grepl("^RP[LS][0-9P]+A?$", all_donor$gene_name) & all_donor$pLI > 0.96,]
# make label for high pLI genes
high_pli_label <- paste("pLI > 0.96 genes (", 
                        length(unique(high_pli$gene_name)),
                        " genes, ",
                        nrow(high_pli),
                        " Exons)",
                        sep="")
# make subset for low pLI genes
low_pli <- all_donor[!grepl("^RP[LS][0-9P]+A?$", all_donor$gene_name) & all_donor$pLI < 0.04,]
# make label for low pLI genes
low_pli_label <- paste("pLI < 0.04 genes (", 
                       length(unique(low_pli$gene_name)),
                       " genes, ",
                       nrow(low_pli),
                       " Exons)",
                       sep="")

# make plots with 19 DBA related RPs
ggplot() +
  geom_line(data=dba_causing, 
            aes(x=end_freq, 
                y=cumsum(..count../nrow(dba_causing)),
                color=dba_causing_label), 
            stat="bin", bins=100000) +
  geom_line(data=not_dba_causing, 
            aes(x=end_freq, 
                y=cumsum(..count../nrow(not_dba_causing)),
                color=not_dba_causing_label), 
            stat="bin", bins=100000) +
  geom_line(data=high_pli, 
            aes(x=end_freq, 
                y=cumsum(..count../nrow(high_pli)),
                color=high_pli_label), 
            stat="bin", bins=100000) +
  geom_line(data=low_pli, 
            aes(x=end_freq, 
                y=cumsum(..count../nrow(low_pli)),
                color=low_pli_label), 
            stat="bin", bins=100000) +
  coord_cartesian(xlim = c(0, 0.8), ylim = c(0.6, 1))+
  scale_color_manual(values = c("red", "black", "grey", "orange")) +
  guides(color = guide_legend("RP Genes by DBA relatedness and other genes by pLI")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = c(0.7, 0.4)) +
  labs(title="Binned Cumulative Line Graph of DBA Related/Non-Related RPs (100000 bins)", 
       x="Sum of Allele Frequency Across Positions in Donor Sites", 
       y="Cumulative Sum of Exons over Total Number of Exons")

# zoom in 0.001
ggplot() +
  geom_line(data=dba_causing, 
            aes(x=end_freq_binned, 
                y=cumsum(..count../nrow(dba_causing)),
                color=dba_causing_label), 
            stat="bin", bins=100) +
  geom_line(data=not_dba_causing, 
            aes(x=end_freq_binned, 
                y=cumsum(..count../nrow(not_dba_causing)),
                color=not_dba_causing_label), 
            stat="bin", bins=100) +
  geom_line(data=high_pli, 
            aes(x=end_freq_binned, 
                y=cumsum(..count../nrow(high_pli)),
                color=high_pli_label), 
            stat="bin", bins=100) +
  geom_line(data=low_pli, 
            aes(x=end_freq_binned, 
                y=cumsum(..count../nrow(low_pli)),
                color=low_pli_label), 
            stat="bin", bins=100) +
  coord_cartesian(ylim = c(0.6, 1))+
  scale_color_manual(values = c("red", "black", "grey", "orange")) +
  guides(color = guide_legend("RP Genes by DBA relatedness and other genes by pLI")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = c(0.7, 0.4)) +
  labs(title="Binned Cumulative Line Graph of DBA Related/Non-Related RPs zoomed on 0.01(100 bins)", 
       x="Sum of Allele Frequency Across Positions in Donor Sites", 
       y="Cumulative Sum of Exons over Total Number of Exons")

