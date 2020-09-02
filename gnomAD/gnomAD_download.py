# downloads data from gnomAD and aggregates necessary data


import requests
import json


# load all RP names and positions
def load_csv():
    rp_list = []
    with open("input_files/RP_positions.csv", "r") as infile:
        for line in infile:
            rp = line.strip("\n").split(",")
            rp_list.append(rp)
    print(rp_list)
    return rp_list


# if the region has too many variants for gnomAD to show
def breakdown(gname, query, ch, start, end, url, headers):
    print("region too big")
    past_midpt = int(start)
    change = (int(end) - int(start)) // 100
    for midpt in range(int(start) + change, int(end) + 1, change):
        json_input = {
            "query": query,
            "variables": {"chrom": f"{ch}",
                          "datasetId": "gnomad_r3",
                          "referenceGenome": "GRCh38",
                          "start": past_midpt,
                          "stop": midpt
                          }
        }
        response = requests.post(url, json=json_input, headers=headers)
        r_json = response.json()
        variants = r_json['data']['region']['variants']
        if variants is None:
            print(r_json)
            print(response.status_code)
            breakdown(gname=gname, query=query, ch=ch, start=past_midpt, end=midpt, url=url, headers=headers)
        past_midpt = midpt

        # saves the downloaded content as csv
        with open(f"input_files/gnomAD_data/{gname}.csv", "a") as outfile:
            for var in variants:
                if var['consequence'] is None:
                    var['consequence'] = "N/A"
                changes = var['variant_id'].split("-")
                outfile.write(var['consequence'] + ',' +
                              changes[2] + ',' +
                              changes[3] + ',' +
                              str(var['pos']) + ',' +
                              str(var['genome']['ac']) + ',' +
                              str(var['genome']['af']) + "\n")


# download data in csv format
def download(rp_list=[]):
    # set url for requests to be sent to
    url = "https://gnomad.broadinstitute.org/api"

    # http request header (constant)
    headers = {"Content-Type": "application/json"}
    # part of the request body (constant)
    query = """
        query VariantInRegion($chrom: String!, $start: Int!, $stop: Int!, $datasetId: DatasetId!, $referenceGenome: ReferenceGenomeId!) {
  region(start: $start, stop: $stop, chrom: $chrom, reference_genome: $referenceGenome) {
    clinvar_variants {
      clinical_significance
      clinvar_variation_id
      gold_stars
      major_consequence
      pos
      variant_id
    }
    variants(dataset: $datasetId) {
      consequence
      flags
      gene_id
      gene_symbol
      hgvs
      hgvsc
      hgvsp
      lof
      lof_filter
      lof_flags
      pos
      rsid
      variant_id: variantId
      exome {
        ac
        ac_hemi
        ac_hom
        an
        af
        filters
        populations {
          id
          ac
          an
          ac_hemi
          ac_hom
        }
      }
      genome {
        ac
        ac_hemi
        ac_hom
        an
        af
        filters
        populations {
          id
          ac
          an
          ac_hemi
          ac_hom
        }
      }
    }
  }
}
        """

    # memory errors :
    # PTPRD, ZNF609, AFF2, RERE
    for rp in rp_list:
        # set which part of the genome we want to download
        print(rp)
        gname = rp[0]
        ch = rp[1]
        start = rp[2]
        end = rp[3]

        # make the complete request body
        json_input = {
            "query": query,
            "variables": {"chrom": f"{ch}",
                          "datasetId": "gnomad_r3",
                          "referenceGenome": "GRCh38",
                          "start": int(start),
                          "stop": int(end)
                          }
        }

        response = requests.post(url, json=json_input, headers=headers)
        print(response.status_code)
        r_json = response.json()

        # saves the downloaded content as csv
        data = r_json['data']['region']['variants']
        if data is None:
            with open(f"input_files/gnomAD_data/{gname}.csv", "w") as outfile:
                outfile.write("Annotation,Before,After,Position,Allele.Count,Allele.Frequency\n")
            breakdown(gname, query, ch, start, end, url, headers)
            continue

        # saves the downloaded content as csv
        with open(f"input_files/gnomAD_data/{gname}.csv", "w") as outfile:
            outfile.write("Annotation,Before,After,Position,Allele.Count,Allele.Frequency\n")
            for var in data:
                if var['consequence'] is None:
                    var['consequence'] = "N/A"
                changes = var['variant_id'].split("-")
                outfile.write(var['consequence'] + ',' +
                              changes[2] + ',' +
                              changes[3] + ',' +
                              str(var['pos']) + ',' +
                              str(var['genome']['ac']) + ',' +
                              str(var['genome']['af']) + "\n")


# use to download just the RP genes from gnomAD
def main():
    rp_list = load_csv()
    download(rp_list)
