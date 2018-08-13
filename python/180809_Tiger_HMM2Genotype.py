import argparse
import pandas as pd
import os
import re

parser = argparse.ArgumentParser(description='From hmm out put to  genotypes')
parser.add_argument('--h', help='The hmm output folder', type=str)
parser.add_argument('--i', help='The raw input folder contain all the markers', type=str)
parser.add_argument('--o', help='A folder to put the outputs', type=str)
args = parser.parse_args()

hmm_folder = args.h
input_folder = args.i
Hmm_file = re.match(".*/(\\d+.*)\\.vcf/", hmm_folder).group(1)
Output_file = args.o + "genotype_" + Hmm_file + ".txt"
Output_file2 = args.o + "genotype2_" + Hmm_file + ".txt"


#input_folder = "/Users/yanjunzan/Documents/impute/results/Tiger_input/999.vcf_DIR/"
#hmm_folder = "/Users/yanjunzan/Documents/impute/results/Tiger_input/999.vcf/"

def match_in_hmm(input_folder, hmm_folder):
    inputs = os.listdir(input_folder)
    hmms = os.listdir(hmm_folder)
    # find inputs
    input_vcf = [i for i in inputs if re.match(".*vcf.*short$",i)]
    hmms_hmm = [i for i in hmms if re.match(".*hmm\\.out\\.txt$",i)]
    id = re.match("(.*)\\.vcf.*",hmms_hmm[1]).group(1)
    input_vcf_chr = [re.match(".*vcf\\.(\\d+)\\.short",i).group(1) for i in input_vcf]
    hmms_hmm_chr = [re.match(".*vcf\\.(\\d+)\\.hmm.*", i).group(1) for i in hmms_hmm]
    chr = set(input_vcf_chr).intersection(hmms_hmm_chr)
    chr = sorted([int(i) for i in chr])
    vcfs = [id + ".vcf." + str(i) + ".short" for i in chr]
    hmms = [id + ".vcf." + str(i) + ".hmm.out.txt" for i in chr]
    d = {key: value for key, value in zip(vcfs, hmms)}
    return d


def comb_hmm_input(Hmm_file, Input_file):
    # read in the raw file
    #Input_file = "/Users/yanjunzan/Documents/impute/results/Tiger_input/999.vcf_DIR/999.vcf.9.short"
    Input = pd.read_csv(Input_file, header=None, sep="\t")
    Input.columns = ["chr", "pos", "H", "allele.h", "L", "allele.l"]
    Nrow = Input.shape[0]

    # Read in the hmm
    #Hmm_file = "/Users/yanjunzan/Documents/impute/results/Tiger_input/999.vcf/999.vcf.9.hmm.out.txt"
    hmm = open(Hmm_file)
    lines = hmm.readlines()
    genotype_raw = lines[2].strip("\n").split(sep=" ")
    genotype = genotype_raw[0:Nrow]
    Input.loc[:, "genotype"] = genotype
    return Input.sort_values(['chr', 'pos'], ascending=[True, True])


output_file = match_in_hmm(input_folder,  hmm_folder)

counter = 0
for k, v in output_file.items():
    vcf = input_folder + k
    hmm = hmm_folder + v
    if counter == 0:
        geno = comb_hmm_input(Hmm_file=hmm, Input_file=vcf)
    else:
        geno = geno.append(comb_hmm_input(Hmm_file=hmm, Input_file=vcf))
    counter = counter + 1
geno.sort_values(['chr', 'pos'], ascending=[True, True]).to_csv(path_or_buf=Output_file, sep="\t", index=False)


## transfer it to a t c g

# geno2 = list()
#
# for i in list(range(0, geno.shape[0])):
#     High = geno.iloc[i, 2] + geno.iloc[i, 2]
#     Low = geno.iloc[i, 4] + geno.iloc[i, 4]
#     HL = geno.iloc[i, 2] + geno.iloc[i, 4]
#     if geno.iloc[i, 6] == "CC":
#         geno2.append(High)
#     elif geno.iloc[i, 6] == "LL":
#         geno2.append(Low)
#     else:
#         geno2.append(HL)
#
# geno.loc[:, "SNP_call"] = geno2
#
# geno.sort_values(['chr', 'pos'], ascending=[True, True]).to_csv(path_or_buf=Output_file2, sep="\t", index=False)
