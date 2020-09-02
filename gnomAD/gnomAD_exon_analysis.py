import os
import re
import random


# reads in and makes csv into a list
def load_csv(path=""):
    output = []
    with open(path, "r") as infile:
        infile.readline()
        for line in infile:
            row = line.strip("\n").split(",")
            for val in range(len(row)):
                try:
                    row[val] = int(row[val])
                except ValueError:
                    pass
            output.append(row)
    return output


# get positions for all genes in gnomAD_data folder
def get_positions():
    # get names of genes to extract positions for
    gnomad_files = os.listdir("input_files/gnomAD_data")
    gnames = [file.split(".", 1)[0] for file in gnomad_files]

    # write the header for output file
    with open("input_files/essential_rp_low_positions.csv", "w") as outfile:
        outfile.write("gene_name,chr,start,end,transcript_name,exon_number\n")

    # open the gtf file and extract position data
    with open("input_files/trimed_gencode.v31.gtf", "r") as gtf:
        # open the output file
        with open("input_files/essential_rp_low_positions.csv", "a") as outfile:
            # loop in the gtf file
            for line in gtf:
                line = line.strip(";\n").replace("\"", "")
                line = re.split("\t|; | ", line)
                # write positions to output file
                if line[2] == "exon" and line[15] in gnames:
                    outfile.write(f"{line[15]},{line[0]},{line[3]},{line[4]},{line[19]},{line[21]}\n")


# use all files in gnomAD_data folder to make one file with all the splice region mutations
def map_splice():
    gnomad_files = os.listdir("input_files/gnomAD_data")

    # define how many base pairs from exon are considered splice region
    inexon = 3
    inintron = (3, 8)

    # load positions of genes
    exon_data = load_csv("input_files/essential_rp_low_positions.csv")
    for i in range(len(exon_data)):
        exon_data[i].extend([0, 0, 0, 0])
    print(exon_data)

    for file in gnomad_files:
        with open(f"input_files/gnomAD_data/{file}", "r") as variant_data:
            gname = file.split(".", 1)[0]
            variant_data.readline()
            # for each variant in the gnomAD data
            for variant in variant_data:
                variant = variant.strip("\n").split(",")
                if variant[0] == "splice_region_variant":
                    # loop in exon position data to see which exon the variant goes to
                    atgene = False  # see if the loop is at the desired gene (used to exit loop)
                    for i in range(len(exon_data)):
                        if exon_data[i] is None:
                            continue
                        if variant[4] == "None":
                            variant[4] = 0
                        elif variant[5] == "None":
                            variant[5] = 0
                        # look at start of exon
                        if exon_data[i][2] - inintron[1] <= int(variant[3]) <= exon_data[i][2] - inintron[0] or \
                                exon_data[i][2] <= int(variant[3]) <= exon_data[i][2] + inexon:
                            exon_data[i][6] += int(variant[4])  # add count
                            exon_data[i][7] += float(variant[5])  # add frequency
                            atgene = True
                        # look at end of exon
                        elif exon_data[i][3] - inexon <= int(variant[3]) <= exon_data[i][3] or \
                                exon_data[i][3] + inintron[0] <= int(variant[3]) <= exon_data[i][3] + inintron[1]:
                            exon_data[i][8] += int(variant[4])  # add count
                            exon_data[i][9] += float(variant[5])  # add frequency
                            atgene = True
                        # if the loop went through the desired gene, break
                        elif atgene and exon_data[i][0] != gname:
                            break

    # write output
    with open("input_files/gnomAD_exon_essential_rp_low.csv", "w") as outfile:
        outfile.write("gene_name,chr,start,end,transcript_name,exon_number,"
                      "start_count,start_freq,end_count,end_freq\n")
        print(exon_data)
        for exon in exon_data:
            if exon is None:
                continue
            outfile.write(",".join(map(str, exon)) + "\n")


# filters for genes with pli greater than 0.96 and for RPs (helper for add_pli)
def pli_filter():
    with open("input_files/gnomad.v2.1.1.lof_metrics.by_gene.txt", "r") as infile:
        header = infile.readline().strip("\n").split("\t")
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


# adds pLI info to csv table
def add_pli():
    with open("input_files/gnomAD_exon_essential_rp_low.csv", "r") as exon_data:
        exon_data.readline()

        # get pLI info
        high_rp = pli_filter()
        low = pli_filter_lower()
        plis = {**high_rp, **low}

        with open("input_files/gnomAD_exon_essential_rp_low_pli.csv", "w") as outfile:
            outfile.write("gene_name,pLI,chr,start,end,transcript_name,exon_number,"
                          "start_count,start_freq,end_count,end_freq\n")
            for line in exon_data:
                try:
                    line = line.strip("\n").split(",")
                    line.insert(1, str(plis[line[0]]))
                    outfile.write(",".join(line) + "\n")
                except KeyError:
                    print(line[0])


# filters for genes with pLI less than 0.04
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


# add strand info on csv
def add_strandedness():
    data = load_csv("input_files/gnomAD_exon_essential_rp_low_pli.csv")

    strand = {}
    with open("input_files/trimed_gencode.v31.gtf", "r") as gtf:
        for line in gtf:
            line = line.replace(";", "").replace(" ", "\t").replace("\"", "").strip("\n").split("\t")
            if line[2] == 'gene' and line[13] not in strand.keys():
                strand[line[13]] = line[6]

    for i in range(len(data)):
        if strand[data[i][0]] == "-":
            endtemp = data[i][4]
            ectemp = data[i][9]
            eftemp = data[i][10]
            data[i][4] = data[i][3]
            data[i][9] = data[i][7]
            data[i][10] = data[i][8]
            data[i][3] = endtemp
            data[i][7] = ectemp
            data[i][8] = eftemp

    with open("input_files/gnomAD_exon_essential_rp_low_pli_stranded.csv", "w") as outfile:
        outfile.write("gene_name,pLI,chr,start,end,transcript_name,exon_number,"
                      "start_count,start_freq,end_count,end_freq\n")
        for line in data:
            try:
                outfile.write(",".join([str(x) for x in line]) + "\n")
            except KeyError:
                print(line[0])


def main():
    # get position of gene/exon in gnomAD_data folder
    get_positions()

    # map variants in gnomAD_data folder to exons
    map_splice()

    # add column of pLI
    add_pli()

    # switch start and end of exon if gene in on negative strand
    add_strandedness()


main()
