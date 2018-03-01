#!/usr/bin/python
import re  # regular expression
import sys
import gzip
#  import pandas

ID_of = re.compile(str(sys.argv[1]))
Fid_h = re.compile(str(sys.argv[2]))
Mid_h = re.compile(str(sys.argv[3]))
Fid_l = re.compile(str(sys.argv[4]))
Mid_l = re.compile(sys.argv[5])
vcf_file = str(sys.argv[6])
out_file = str(sys.argv[7])
#ID_of = re.compile(str(".*1325.*"))
#Fid_h = re.compile(str(".*1690.*"))
#Mid_h = re.compile(str(".*1812.*"))
#Fid_l = re.compile(str(".*1997.*"))
#Mid_l = re.compile(str(".*1925.*"))

# ID_of = re.compile("175")
# Fid = re.compile(".*1997.*")
# Mid = re.compile(".*1812.*")
#vcf_file = "/Users/yanjunzan/Documents/impute/git/data/180208.all.223+700.f2.P60.vcf.gz"
#out_file = "/Users/yanjunzan/Documents/impute/git/data/test.180208.vcf"
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
            #print(line)
        elif re.match(b"^#CHROM", line):
            line = line.decode().rstrip("\n")
            header = re.split("\t", line)
            index_off = find_index(ID_of, header)
            index_F_h = find_index(Fid_h, header)
            index_M_h = find_index(Mid_h, header)
            index_F_l = find_index(Fid_l, header)
            index_M_l = find_index(Mid_l, header)
            #  print(index_off, index_F, index_M)
        elif re.match(b"\n",line):
            pass
        else:
            #print(line,"is here")
            line = line.decode().rstrip("\n")
            #print(line)
            con = re.split("\t", line)
            Chr = con[0]
            Pos = con[1]
            ID = con[2]
            REF = con[3]
            ALT = con[4]
            if len(REF) > 1 or len(ALT) > 1 or re.match("\.", con[index_off[0]]):
                continue
            else:
                geno_of1 = int(get_geno(con[index_off[0]])[0])
                geno_of2 = int(get_geno(con[index_off[0]])[1])
                geno_F1_h = int(get_geno(con[index_F_h[0]])[0])
                geno_F2_h = int(get_geno(con[index_F_h[0]])[1])
                geno_M1_h = int(get_geno(con[index_M_h[0]])[0])
                geno_M2_h = int(get_geno(con[index_M_h[0]])[1])
                geno_F1_l = int(get_geno(con[index_F_l[0]])[0])
                geno_F2_l = int(get_geno(con[index_F_l[0]])[1])
                geno_M1_l = int(get_geno(con[index_M_l[0]])[0])
                geno_M2_l = int(get_geno(con[index_M_l[0]])[1])
                if sum([geno_F1_h, geno_F2_h, geno_M1_h, geno_M2_h, geno_F1_l, geno_F2_l, geno_M1_l, geno_M2_l]) == 4 and (sum([geno_F1_h, geno_F2_h, geno_M1_h, geno_M2_h]) == 0 or sum([geno_F1_h, geno_F2_h, geno_M1_h, geno_M2_h]) == 4 ):
                    print(Chr, Pos, REF, ALT, geno_F1_h, geno_M2_h, geno_F1_l, geno_M2_l, geno_of1, geno_of2, sep="\t")

