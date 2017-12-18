#!/usr/bin/python
import re  # regular expression
import sys
import gzip
#  import pandas

ID_of = re.compile(str(sys.argv[1]))
Fid = re.compile(str(sys.argv[2]))
Mid = re.compile(sys.argv[3])
vcf_file = str(sys.argv[4])
out_file = str(sys.argv[5])
# ID_of = re.compile("175")
# Fid = re.compile(".*1997.*")
# Mid = re.compile(".*1812.*")
# vcf_file = "/Users/yanjunzan/Documents/impute/data/testplate.F1/chr1_CM000093.4.recode.vcf.gz"
sys.stdout = open(out_file, "wt")


def find_index(name, lines):
    indices = []
    for i, elem in enumerate(lines):
        if re.match(name, elem):
            indices.append(i)
    return indices


def get_geno(x):
    if re.match("^\.", x):
        return [5, 5]
    else:
        out = re.match("(\d{1}).{1}(\d{1}).*", x)
        return out.group(1, 2)


with gzip.open(vcf_file, "rb") as vcf:
    for line in vcf:
        if re.match(b"^##", line):
            pass
        elif re.match(b"^#CHROM", line):
            line = line.decode().rstrip("\n")
            header = re.split("\t", line)
            index_off = find_index(ID_of, header)
            index_F = find_index(Fid, header)
            index_M = find_index(Mid, header)
            #  print(index_off, index_F, index_M)
        else:
            line = line.decode().rstrip("\n")
            con = re.split("\t", line)
            Chr = con[0]
            Pos = con[1]
            ID = con[2]
            REF = con[3]
            ALT = con[4]
            if len(REF) > 1 or len(ALT) > 1:
                pass
            elif re.match("\.", con[index_off[0]]):
                pass
            elif re.match("\.", con[index_F[0]]):
                pass
            elif re.match("\.", con[index_M[0]]):
                pass
            else:
                geno_of1 = int(get_geno(con[index_off[0]])[0])
                geno_of2 = int(get_geno(con[index_off[0]])[1])
                geno_F1 = int(get_geno(con[index_F[0]])[0])
                geno_F2 = int(get_geno(con[index_F[0]])[1])
                geno_M1 = int(get_geno(con[index_M[0]])[0])
                geno_M2 = int(get_geno(con[index_M[0]])[1])
                if not (sum([geno_F1, geno_F2, geno_M1, geno_M2]) == 0 or sum([geno_F1, geno_F2, geno_M1, geno_M1]) == 4):
                    print(Chr, Pos, REF, ALT, geno_F1, geno_F2, geno_M1, geno_M2, geno_of1, geno_of2, sep="\t")
