# used to download gnomAD data genewise
# also used to analyze genewise


import re
import os
import random

import gnomAD_download


# filters for genes with pli greater than 0.96 and for RPs
def pli_filter():
    with open("input_files/gnomad.v2.1.1.lof_metrics.by_gene.txt", "r") as infile:
        header = infile.readline().strip("\n").split("\t")
        print(header.index("pLI"))
        print(header[header.index("pLI")])
        outdict = {}
        for line in infile:
            line = line.strip("\n").split("\t")
            name = line[0]
            if line[header.index("pLI")] == "NA":
                continue
            pli = float(line[header.index("pLI")])
            if pli > 0.96 or re.search(r"^RP[SL][0-9P]+A?$", name):
                outdict[name] = pli
    return outdict


# get information about selected genes
def get_gene_info(high_pli={}):
    with open("input_files/trimed_gencode.v31.gtf", "r") as infile:
        outdict = {}
        for line in infile:
            line = line.replace(";", "").replace(" ", "\t").replace("\"", "").strip("\n").split("\t")
            # find gene line in our high_pli list
            if line[2] == 'gene' and line[13] in high_pli.keys():
                # gname, chr, start, end, gene_len, pLI
                gene_out = [line[13], line[0], line[3], line[4], int(line[4]) - int(line[3]), high_pli[line[13]]]
                outdict[line[13]] = gene_out
    return outdict


# aggregates gnomAD data into a csv file (genewise)
def combine(gene_infos={}):
    outlist = [["gene", "chr", "start", "end", "gene_len", "pLI"]]
    for gene in gene_infos.keys():
        # add the gene info to outlist
        outlist.append(gene_infos[gene])
        zeros_to_add = len(outlist[0]) - 7
        for n in range(zeros_to_add):
            outlist[len(outlist) - 1].append(0)

        # look in gnomAD data
        filename = "input_files/gnomAD_data/" + gene + ".csv"
        with open(filename, "r") as infile:
            infile.readline()
            for line in infile:
                # read in a line
                line = line.strip("\n").split(",")

                # add column for that variant type if that type isn't in the list
                # and add 0s to columns that dont have the type
                if line[0] not in outlist[0]:
                    outlist[0].append(line[0])
                    for i in range(1, len(outlist)):
                        if len(outlist[i]) < len(outlist[0]):
                            outlist[i].append(0)

                # add count to the last row aka the current row
                if line[5] == "None":
                    line[5] = 0
                outlist[len(outlist) - 1][outlist[0].index(line[0])] += float(line[5])

    # write result to csv file
    with open("input_files/gnomAD_genefreq_sum.csv", "w") as outfile:
        for line in outlist:
            outfile.write(",".join(map(str, line)) + "\n")


# filters for genes below pLI 0.04
def pli_filter_lower():
    # number of nonRP high pLI genes in pLI data file= 2483
    high_pli_count = 2483
    with open("input_files/gnomad.v2.1.1.lof_metrics.by_gene.txt", "r") as infile:
        header = infile.readline().strip("\n").split("\t")
        tempdict = {}
        for line in infile:
            line = line.strip("\n").split("\t")
            name = line[0]
            if line[header.index("pLI")] == "NA":
                continue
            pli = float(line[header.index("pLI")])
            if pli < 0.04 and not re.search(r"^RP[SL][0-9P]+A?$", name):
                tempdict[name] = pli

    # see which genes are in the gtf file as well
    outdict = {}
    with open("input_files/trimed_gencode.v31.gtf", "r") as gtf:
        for line in gtf:
            line = line.replace(";", "").replace(" ", "\t").replace("\"", "").strip("\n").split("\t")
            if line[2] == 'gene' and line[13] in tempdict.keys():
                outdict[line[13]] = tempdict[line[13]]

    # randomly choose 2483 genes
    random.seed(12345)
    out_gnames = random.sample(list(outdict), high_pli_count)
    outdict = {gname: outdict[gname] for gname in out_gnames}

    return outdict


def main():
    high_pli = pli_filter()
    gene_infos_high = get_gene_info(high_pli)

    low_pli = pli_filter_lower()
    gene_infos_low = get_gene_info(low_pli)

    gene_infos = dict(gene_infos_low) + dict(gene_infos_high)
    gene_infos_copy = dict(gene_infos_low) + dict(gene_infos_high)

    already_downloaded = os.listdir("input_files/gnomAD_data")
    already_downloaded = [genename.split(".", 1)[0] for genename in already_downloaded]
    # remove genes that are already downloaded
    for genename in already_downloaded:
        try:
            del gene_infos_copy[genename]
        except KeyError:
            pass
    print(len(gene_infos_copy))
    gnomAD_download.download([gene[0:4] for gene in gene_infos_copy.values()])
    # combine(gene_infos)


main()
