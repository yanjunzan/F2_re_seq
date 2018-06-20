import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description='Create a file containing the parental and offspring genotypes')
parser.add_argument('--c', help='The refined crossover file', type=str)
parser.add_argument('--i', help='The raw input file contain all the markers', type=str)
parser.add_argument('--o', help='A folder to put the outputs', type=str)
args = parser.parse_args()

CO_file = args.c
Input_file = args.i
Output_file = args.o + "genotype_" + CO_file

#Output_file = "/Users/yanjunzan/Documents/impute/HMM/TIGER/TIGER_Scripts-for-distribution/yanjun_test/"+"genotype_" + "test.txt"
# read in CO file

#CO_file = "/Users/yanjunzan/Documents/impute/HMM/TIGER/TIGER_Scripts-for-distribution/yanjun_test/rough_COs_1325.refined.breaks.txt"
CO = pd.read_csv(CO_file, header=None, sep="\t")
CO.columns = ["sample", "chr", "start", "end", "genotype"]

# read in position file
#Input_file = "/Users/yanjunzan/Documents/impute/HMM/TIGER/TIGER_Scripts-for-distribution/input/1325_chr1.parts.txt"
Input = pd.read_csv(Input_file, header=None, sep="\t")
Input.columns = ["chr", "pos", "H", "allele.h", "L", "allele.l"]

# generate holder
geno = ["NA"] * Input.shape[0]
Input.loc[:, "genotype"] = geno
# fix first line
chr_all = pd.Series(CO["chr"], dtype="category")
chromosomes = chr_all.cat.categories

for j in chromosomes:
    C_chr = CO["chr"] == j
    I_chr = Input["chr"] == j
    for i in np.arange(0, CO.loc[C_chr].shape[0]):
        A = Input.loc[I_chr, "pos"] >= CO.loc[C_chr].iloc[i, ]["start"]
        B = Input.loc[I_chr, "pos"] <= CO.loc[C_chr].iloc[i, ]["end"]
        Input.loc[A & B & I_chr, "genotype"] = str(CO.loc[C_chr].iloc[i, ]["genotype"])

Input.sort_values(['chr', 'pos'], ascending=[True, True]).to_csv(path_or_buf=Output_file, sep="\t", index=False)

