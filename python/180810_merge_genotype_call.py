import argparse
import pandas as pd
import os
import re
from itertools import chain
parser = argparse.ArgumentParser(description='Create a file containing the parental and offspring genotypes')
#parser.add_argument('--ID_of', help='A list of offsprings names', type=str, nargs='+')
parser.add_argument('--d', help='directory for all the individuals', type=str)
parser.add_argument('--f', help='location file', type=str)
parser.add_argument('--o', help='A folder to put the outputs', type=str)

args = parser.parse_args()

Folder = args.d
File = args.f
Output_file = args.o + "genotype_matrix_.txt"

#Folder = "/Users/yanjunzan/Documents/impute/results/Tiger_input/"
#File = "/Users/yanjunzan/Documents/impute/results/Within_fam_all_mrk3.txt"

## find all files first

Files_all = [i for i in os.listdir(Folder) if re.match("genotype_.*", i)]

# creat a dic

all_pos = {}

with open(File) as POS:
    for line in POS:
        chr, pos = line.strip("\n").split("\t")
        key = chr + "_" + pos
        all_pos[key] = [chr, pos]


# cread a pd dataframe for the same size
Row = len(all_pos)
Col = 3 + len(Files_all)
key_all= list(all_pos.keys())

Holder = pd.DataFrame(index=list(all_pos.keys()), columns=range(Col))
a = ["key","chr", "pos"]
a.extend(Files_all)
Holder.columns = a

Holder.loc[:, "key"] = key_all
# open and read in the first file

# for i in range(len(Files_all)):
#
#     with open(Folder + Files_all[i]) as IN:
#         for line in IN:
#             chr, pos, H, HL, L, LL, g = line.strip("\n").split("\t")
#             key_now = chr + "_" + pos
#             #print(key_now)
#             Holder.loc[key_now, "chr"] = chr
#             Holder.loc[key_now, "pos"] = pos
#             Holder.loc[key_now, Files_all[i]] = g
#
# Holder.sort_values(['chr', 'pos'], ascending=[True, True]).to_csv(path_or_buf=Output_file, sep="\t", index=False)
#


for i in range(len(Files_all)):
    key_now = list()
    with open(Folder + Files_all[i]) as IN:
        for line in IN:
            chr, pos, H, HL, L, LL, g = line.strip("\n").split("\t")
            key_now.append( chr + "_" + pos)
            #print(key_now)
#             Holder.loc[key_now, "chr"] = chr
#             Holder.loc[key_now, "pos"] = pos
#             Holder.loc[key_now, Files_all[i]] = g
#
# Holder.sort_values(['chr', 'pos'], ascending=[True, True]).to_csv(path_or_buf=Output_file, sep="\t", index=False)
#


st = [i for i, e in enumerate(key_all) if e in key_now]

