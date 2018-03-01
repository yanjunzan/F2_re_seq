##################################################################################
####### Script for reading in input files from QTLFast/Haley-Knott ###############
####### Author: Örjan Carlborg, Date: 180226                       ###############
####### 180226_ReadQTLFastFiles.R                                  ###############
##################################################################################

# File with genotypes from Per Wahlberg for HWS x LWS intercross at generation 41 
RawGenotypesF2 <- read.table(file="~/Dropbox/Projects/Ongoing/HL_F2_SeqMrk/F2_re_seq/data/F2/hl061106.mrkdata",header=FALSE,sep=" ")

# File with information on where the markers on the different chromosomes are located in the marker data file
Mrkinfo <- read.table(file="~/Dropbox/Projects/Ongoing/HL_F2_SeqMrk/F2_re_seq/data/F2/hl061106_mrkinfo_forR.csv",header=FALSE,sep=";")

# File with location of genotyped markers on the linkage map published in Wahlberg, 2009 in BMC Genomics
Mapinfo <- read.table(file="~/Dropbox/Projects/Ongoing/HL_F2_SeqMrk/F2_re_seq/data/F2/WahlbergLinkageMap2009.txt",header=TRUE,sep="\t")
Mapinfo <-  Mapinfo[,1:4]
# Note: Discrepancy between files as mrk number 47 on chromosome 1 in the map used for QTL mapping is not 
# included in the  Wahlberg 2009 linkage map

# How many chromosomes & markers have data?
nochr <- Mrkinfo[1,1]
nomrk <- Mrkinfo[1,2]

# Get the information on which columns contain the genotypes for the respective markers
## How many markers on each chromosome
chrom_mrks <- Mrkinfo[2:(nochr+1),1]
## Read in the markerlocations
mrks_chrom <- matrix(nrow = nochr, ncol = max(chrom_mrks))
mrks_chrom <- Mrkinfo[2:(nochr+1),2:(max(chrom_mrks)+1)]
## Rownames with chromsome numbers
rownames(mrks_chrom) <- c(1:nochr)
