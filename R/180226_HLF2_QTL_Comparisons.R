#################### Author: Örjan Carlborg, IMBIM, Uppsala University ##############################
# Date of last revision: 180227

######################################## Description: ###############################################
# R-code for comparing results from analysis of the F2 intercross population between HWS and LWS
# birds from generation 40. The original results from these analyses are published in two report
# first based on a sparse Micro-satellite based linkage map in:
# Jacobsson, L., Park, H.-B., Wahlberg, P., Fredriksson, R., Perez-Enciso, M., Siegel, P. B., & Andersson, L. (2005). Many QTLs with minor additive effects are associated with a large difference in growth between two selection lines in chickens. Genetical Research, 86(02), 115. doi:10.1017/S0016672305007767
# and later after adding 300+ SNP markers in
# Wahlberg, P., Carlborg, Ö., Foglio, M., Tordoir, X., Syvänen, A.-C., Lathrop, M., et al. (2009). Genetic analysis of an F(2) intercross between two chicken lines divergently selected for body-weight. BMC Genomics, 10(1), 248. doi:10.1186/1471-2164-10-248
#
# Analyses here are based on the data/qtl probabilities from Wahlberg, 2009

########### The following processing ofdata and analyses are performed in this script:
# Section A. Reading in the data
# Section B. Further processing of the raw data
# Section C. Comparing genotypes (GBSA vs SNP/MS) & QTL genotype probabilities

#########################################################################################
##################### Section A. Reading in data ########################################

####### A.1 Reading in phenotypes for F2 population  ########
# Data used in Wahlberg, 2009 & Jacobsson, 2005  studies    #
# Input-file from Wahlberg, 2009 study                      #

Phenotypes <- read.table(file="~/Dropbox/Projects/Ongoing/HL_F2_SeqMrk/F2_re_seq/data/F2/hl071022.phen",header=TRUE,sep="\t")
tst <- Phenotypes<(-100000)
Phenotypes[tst]<-NA
rownames(Phenotypes)<-Phenotypes[,1]
Phenotypes <- Phenotypes[,2:ncol(Phenotypes)]

######## A.2. Reading in fixed effects for F2 population ########
# Data  used Wahlberg, 2009 & Jacobsson, 2005 studies           #
# Input-file from Wahlberg, 2009 study                          #

Fixed <- read.table(file="~/Dropbox/Projects/Ongoing/HL_F2_SeqMrk/F2_re_seq/data/F2/hl071022.fix",header=TRUE,sep=" ")
rownames(Fixed) <- Fixed[,1]
Fixed <- Fixed[,2:ncol(Fixed)]

####### A.3. Reading in genome-wide genotype-probabilities                 ####
# F2 QTL genotype probabilities, calculated using Haley et al, 2004 algorithm #
# Input-file from Wahlberg, 2009 study                                        #
# Genotype probabilities were then calculated using the QTLExpress software   #

Geno <- read.table(file="~/Dropbox/Projects/Ongoing/HL_F2_SeqMrk/F2_re_seq/data/F2/hl070831.khcoeffqe",header=FALSE,sep="\t")
Geno[,8] <- Geno[,4]-Geno[,7] # Calculate the indicator regression variable as in Haley et al, 2004 (P(QQ)-p(qq), where Q=HWS allele, q=LWS allele)

###### A.4. Reading in Genetic- and Physical map (Wahlberg, 2009)       ########## 
# Map is available in the Supplements of the paper Wahlberg, 2009                #
# Input-file created by Örjan Carlborg based on the published genetic map        #
# Physical locations in this file are given on the GalGal3 assembly              #
# One missmatch on Chromosome 27 detected -> marker moved 0.1 cM in marker map   #

MrksF2 <- read.table(file="~/Dropbox/Projects/Ongoing/HL_F2_SeqMrk/F2_re_seq/data/F2/150202_OCG_hlF2markermapWahlberg2009.txt",header=TRUE,sep="\t")
MrksF2[,1] <- as.numeric(as.character(MrksF2[,1]))
MrksF2[,2] <- as.numeric(as.character(MrksF2[,2]))
colnames(MrksF2)[2]<-c("GalGal3Mb")
MrksF2[,3] <- as.numeric(as.character(MrksF2[,3]))
MrksF2[,4] <- as.character(MrksF2[,4])
# Chr 27, 59 cM does not exist in mrk probabilities -> move to 58.9 cM
MrksF2[which(MrksF2[,1]==27&MrksF2[,3]==59),3]<-58.9
MrksF2[,5] <- as.integer(floor(MrksF2[,3]))
colnames(MrksF2)[5]<-c("cMInt")
MrksF2[,6] <- MrksF2[,2]*1000000 # Convert Mb locations to bp
MrksF2[,6] <- as.numeric(as.character(MrksF2[,6]))
colnames(MrksF2)[6]<-c("GalGal3bp")
MrksF2[,4] <- gsub('-', '_', MrksF2[,4])
MrksF2<-MrksF2[!is.na(MrksF2[,1]),]

#########################################################################################
################### Section B. Further processing of the raw data               ######### 

################### B1. Liftover of marker locations GalGal3 to GalGal5         #########
########## SNP/MS in Wahlberg, 2009 on GalGal3, GBSA markers on GalGal5         #########
########## To integrate data, SNP/MS marker locations of Wahlberg et al 2009    ######### 
########## are transferred to GalGal5                                           #########

library(rtracklayer)
liftdata <- MrksF2[!is.na(MrksF2[,6]),c(1,4,6)]
liftdata[,1] <- paste('chr',liftdata[,1],sep="")
liftnames <- data.frame(liftdata[,2])
liftnames[,1]<-as.character(liftnames[,1])
colnames(liftnames)<-c("Marker")
GG3 <- GRanges(seqnames=Rle(c(liftdata[,1])), ranges=IRanges(start=c(liftdata[,3]),end=c(liftdata[,3]),width=NULL))
mcols(GG3,use.names=TRUE) <- liftnames
chain <- import.chain("~/Dropbox/Projects/Ongoing/HL_F2_SeqMrk/F2_re_seq/data//liftover/galGal3ToGalGal5.over.chain")
GG5 <- liftOver(GG3, chain)
GG5 = unlist(GG5)
length(GG3)-length(GG5)
setdiff(GG3$Marker, GG5$Marker)
tmp<-c()
tmp<-cbind(as.character(GG5$Marker),as.character(seqnames(GG5)),as.numeric(round(start(ranges(GG5))/1000000,1)))
tmp2<-c(tmp)
MrksF2[MrksF2[,4]%in%tmp[,1],7]<-tmp[,3]
colnames(MrksF2)[7]<-c("GalGal5Mb")
LocMrksF2<-MrksF2[,c(1,3,7,2,4)]

####### B.2. Manual assignment of GalGal5 location for markers where liftover failed #######
####### Locations in file found by manually searching the NCBI database in February, 2018  #

LocManMrksF2 <- read.table(file="~/Dropbox/Projects/Ongoing/HL_F2_SeqMrk/F2_re_seq/data/F2/180226_OCG_hlF2markersWahlberg2009ManualPositionNCBI.txt",header=TRUE,sep="\t")
LocManMrksF2[,1] <- as.character(LocManMrksF2[,1])
LocManMrksF2[,4] <- as.numeric(LocManMrksF2[,4])
for (i in 1:nrow(LocManMrksF2)){
  LocMrksF2[which(LocMrksF2[,5]==LocManMrksF2[i,1]),3] <- LocManMrksF2[i,4]
  MrksF2[which(MrksF2[,4]==LocManMrksF2[i,1]),7] <- LocManMrksF2[i,4]
}

# For some markers, the location on the GG4 assembly could not be found, neither through
# automatic liftover or manual searches. These are removed from subsequent analyses.
MrksF2 <- MrksF2[which(!is.na(MrksF2[,7])),]
MrksF2[,7]<-as.numeric(MrksF2[,7])
MrksF2[,8]<-paste("Geno_a$",MrksF2[,4],sep="")
MrksF2[,9]<-NA
LocMrksF2 <- LocMrksF2[which(!is.na(LocMrksF2[,3])),]

##### B.3. Build data-structure with F2 QTL-genotype probabilities for regression analysis ####
##### Probabilities from Wahlberg SNP/MS dataset (2009) read in above

Geno_a <- matrix(data=NA,nrow=length(unique(Geno[,2])))
Geno_a[,1]<-unique(Geno[,2])
rownames(Geno_a) <- Geno_a[,1]
for (i in 1:nrow(MrksF2)){
  Geno_a <- cbind(Geno_a,Geno[((Geno[,1]==MrksF2[i,1]) & (Geno[,3]==MrksF2[i,5])),8])
}
Geno_a <- Geno_a[,2:ncol(Geno_a)]
mrks <- MrksF2[,4]
colnames(Geno_a)<-mrks
Geno_a <- as.data.frame(Geno_a)
colnames(Geno_a) <- gsub('-', '_', colnames(Geno_a))

##### B.4. Read in datastructure with the raw SNP/MS genotypes from Per Wahlberg  ###########
## This is done by running the script 180226_ReadQTLFastFiles.R
source(file ="~/Dropbox/Projects/Ongoing/HL_F2_SeqMrk/F2_re_seq/R/180226_ReadQTLFastFiles.R")

##### B.5. Load the Rqtl object mapthis including genotype, phenotype and linkage map #######
## data for the GBSA marker data

load(file='~/Dropbox/Projects/Ongoing/HL_F2_SeqMrk/F2_re_seq/data/HL_F2_mapthis.Rdata')

############################################################################################################
########### Section C. Comparing GBSA & SNP/MS marker genotypes / QTL genotype probabilities     ###########

############################################################################################################
########## C.1. Extracting genotypes for specific markers from the SNP/MS data as reference      ###########

#### Example, the nearest markers to Growth1 QTL on GGA1 is rs15502284
qtl_chromosome <- 1
qtl_marker <- 'rs15502284'

#### Find the location for this in the genotype object containing the raw SNP/MS genotype information
if (qtl_chromosome > 1) {
  qtl_marker_location <- mrks_chrom[qtl_chromosome,which(Mapinfo$Marker==qtl_marker)]
} else if (qtl_chromosome == 1) {
  qtl_marker_location <- mrks_chrom[qtl_chromosome,which(Mapinfo$Marker==qtl_marker)]
  if (qtl_marker_location > 46) {
    # On GGA1 '+1' is used for markers with ID >=47 as mrk 47 on the genetic map 
    # was not included in the linkage map published in Wahlberg 2009
    qtl_marker_location <- mrks_chrom[qtl_chromosome,(which(Mapinfo$Marker==qtl_marker)+1)]
  }
}

#### Extract the SNP/MS genotypes for all individuals at this location
## Genotypes given as: 'ID, allele 1, allele 2'
qtl_mrk_genotypes <- RawGenotypesF2[c(1,2*qtl_marker_location,2*qtl_marker_location+1)]

#### Extracting SNP/MS genotypes only for a subset of the individuals (e.g. a FS family)
## Example, FS family including F2's 353 & 358
#1997 0 0 1 2
#1812 0 0 0 1
#1690 0 0 1 1
#1925 0 0 0 2
#175 1997 1812 1 
#294 1690 1925 0 
#353 175 294 1 
#358 175 294 1 

family  <- c(1997,1812,1690,1925,175,294,353,358)
qtl_mrk_family <- qtl_mrk_genotypes[qtl_mrk_genotypes[,1]%in%family,1:3]

############################################################################################################
########### C.2. Extract genotype probabilities for an individual across an entire chromosome ######## 

##### C.2.1. Function to prepare data for plotting                                                        ########

prepare_plot_prob <- function (plot_chromosome,plot_name) {
  ## First get GBSA probabilities
  tmp1 <- eval(parse(text = paste("mapthis$geno$",plot_chromosome,"$prob[,,1]",sep="'")))
  tmp2 <- eval(parse(text = paste("mapthis$geno$",plot_chromosome,"$prob[,,3]",sep="'")))
  rqtlprob_gbsa<-tmp1-tmp2
  rqtlprob_gbsa<-as.data.frame(rqtlprob_gbsa)
  rownames(rqtlprob_gbsa)<-mapthis$pheno[,"F2_ID"]

  ## then locations of GBSA markers
  gbsa_cms <- eval(parse(text = paste("as.vector(unname(mapthis$geno$",plot_chromosome,"$map))",sep="'")))
  gbsa_locations<-eval(parse(text = paste('(as.numeric(gsub(',plot_name,',"",names(mapthis$geno$',plot_chromosome,'$map))))',sep="'")))
  colnames(rqtlprob_gbsa)<-gbsa_locations

  ## Second get locations for SNP/MS markers 
  chrom_mrknames <- Mapinfo[which(Mapinfo$Chromosome==plot_chromosome),1]
  which(colnames(Geno_a) %in% chrom_mrknames)
  Chrom_geno <- Geno_a[which(as.numeric(rownames(Geno_a)) %in% mapthis$pheno[,"F2_ID"]),which(colnames(Geno_a) %in% chrom_mrknames)]
  snpms_locations<-(as.numeric(LocMrksF2$GalGal5Mb[which(LocMrksF2$Chromosome==plot_chromosome)]))
  colnames(Chrom_geno)<-snpms_locations

  ## Prepare for plotting the probabilities across the chromosome
  # Create variables that contain one column for every scored marker 
  tmp1 <- rqtlprob_gbsa
  tmp1[] <- 'NA'
  tmp2 <- Chrom_geno
  tmp2[] <- 'NA'

  # One for the GBSA and one for the MS/SNP markers, columns empty if marker of other kind
  plot_gbsa_prob <- cbind(rqtlprob_gbsa,tmp2)
  plot_qtl_prob <- cbind(Chrom_geno,tmp1)
  # Sort them to plot markers in correct order
  plot_gbsa_prob <- plot_gbsa_prob[ , order(as.numeric(colnames(plot_gbsa_prob)))]
  plot_qtl_prob <- plot_qtl_prob[ , order(as.numeric(colnames(plot_qtl_prob)))] 
}

##### C.2.2 Function for plotting genotype probabilities across an entire chromosome for a given individual
## The individual to plot is provided in the 'check_this_id' variable and the chromosome in 'check_this_chrom'
## Output are plots with genotype-probabilities from SNP/MS data (qtlprobs) and GBSA data (rqtlprobs)
checkid <- function (check_this_id,check_this_chrom)  {
  qtlprobs<-Geno$V8[which(Geno$V2==check_this_id & Geno$V1==check_this_chrom)]
  tmp1 <- eval(parse(text = paste("mapthis$geno$",check_this_chrom,"$prob[which(mapthis$pheno[1]==check_this_id),,1]",sep="'")))
  tmp2 <- eval(parse(text = paste("mapthis$geno$",check_this_chrom,"$prob[which(mapthis$pheno[1]==check_this_id),,3]",sep="'")))
  rqtlprob<-as.numeric(tmp1-tmp2)
  plot(rqtlprob)
  plot(qtlprobs)
}

#### Provide ID & Chromosome & run 'checkid' function 
## this will output SNP/MS & GBSA qtl probabilities for individual 
ID_tocheck <- 353
Chrom_tocheck <- 1
checkid(ID_tocheck,Chrom_tocheck)

### ToDo: plot probabilities at scored markers  & corresponding Mb locations

##### Checking locations across 'plot_chromosome' for individuals

chromosome_to_plot <- 2

Chrom_names <- read.table(file="~/Dropbox/Projects/Ongoing/HL_F2_SeqMrk/F2_re_seq/data/180228_Contig_names.txt",sep =" ",stringsAsFactors = FALSE)
Chrom_name <- Chrom_names$V1[chromosome_to_plot]

prepare_plot_prob(chromosome_to_plot,Chrom_name)

### And now plot
# Whole chromosome
dev.off()
plot_individual <- 4
plot(as.numeric(colnames(plot_qtl_prob)),as.numeric(plot_qtl_prob[plot_individual,]), col = "blue", pch = 2, cex = 1)
points(as.numeric(colnames(plot_gbsa_prob)),as.numeric(plot_gbsa_prob[plot_individual,]),col = "red", pch = 1, cex = 0.5)

# Region on chromosome
plot_region <- c(165,175) #Start & stop of region in Mb (GalGal5)

mrk_to_plot <- which((as.numeric(colnames(plot_qtl_prob))>plot_region[1])&(as.numeric(colnames(plot_qtl_prob))<plot_region[2]))
plot(as.numeric(colnames(plot_qtl_prob[mrk_to_plot])),as.numeric(plot_qtl_prob[plot_individual,mrk_to_plot]), col = "blue", pch = 2, cex = 1)
points(as.numeric(colnames(plot_gbsa_prob[mrk_to_plot])),as.numeric(plot_gbsa_prob[plot_individual,mrk_to_plot]),col = "red", pch = 1, cex = 0.5)

#### Compare genotype probabilities across chromosome #####
max(as.numeric(colnames(plot_qtl_prob)),na.rm=TRUE)

head(rqtlprob_gbsa)
head(Chrom_geno)



############################################################################################################
########### C.3. Study genotype probabilities at particular locations, e.g. QTL                     ######## 

## Identify individuals with deviating SNP/MS  & GBSA genotypes at a particular location
deviants <- QTL1_SNP_gbsa[which(abs(QTL1_SNP_gbsa[,2]-QTL1_SNP_gbsa[,3])>0.5),1]
ID_tocheck <- deviants[81]

## Fewer individuals have GBSA genotypes. Find which these are 
# and which of them have phenotypes allowing them to be used for QTL mapping
used_ids_gbsa <- fake.f2$pheno$F2_ID[which(!is.na(fake.f2$pheno[2]))]

# Extract the genotype probabilities
QTL1_GBSA<-qtl$prob[[1]][which(!is.na(fake.f2$pheno[2])),1]-qtl$prob[[1]][which(!is.na(fake.f2$pheno[2])),3]
#
QTL1_SNP_gbsa <- cbind(QTL1_SNP[which(QTL1_SNP[,1]%in%used_ids_gbsa),],QTL1_GBSA)
tst<-lm(fake.f2$pheno$BW8[which(!is.na(fake.f2$pheno[2]))]~QTL1_SNP_gbsa[,2])
summary(tst)
cor(QTL1_SNP_gbsa[,2],QTL1_SNP_gbsa[,3])
plot(QTL1_SNP_gbsa[,2]-QTL1_SNP_gbsa[,3])
QTL1_SNP_gbsa[which(abs(QTL1_SNP_gbsa[,2]-QTL1_SNP_gbsa[,3])>0.5),1]


#Pull out probabilities for Growth1
QTL1_SNP<-cbind(Geno$V2[which(Geno$V3=='448' & Geno$V1==1)],Geno$V8[which(Geno$V3=='448' & Geno$V1==1)])

QTL1_SNP_gbsa[QTL1_SNP_gbsa[,1]==353]
QTL1_SNP_gbsa[QTL1_SNP_gbsa[,1]==358]


## For example id 353
ID353_c1<-Geno$V8[which(Geno$V2=='353' & Geno$V1==1)]
id353_rqtlprob<-mapthis$geno$`1`$prob[4,,1]-mapthis$geno$`1`$prob[4,,3]
id353_rqtlprob<-as.numeric(id353_rqtlprob)
mapthis$geno$`1`$map[1]


plot(id353_rqtlprob)
plot(ID353_c1[435:485]) # 167.6-185.1
plot(id353_rqtlprob[121:157]) #Peak 136
ID353_c1[400:498]
ID353_c1[448]
plot(id353_rqtlprob[136:146])
mapthis$geno$`1`$map[136:146]  #175-185 Mb




mapthis$geno$`1`$map[[3]]

##### Checking raw GBSA genotypes for the markers surrounding the Growth 1 peak ########
mapthis$geno$`1`$data[mapthis$pheno[1]==353,135]
mapthis$geno$`1`$data[mapthis$pheno[1]==353,136]
mapthis$geno$`1`$data[mapthis$pheno[1]==353,137]

mapthis$geno$`1`$data[mapthis$pheno[1]==358,135]
mapthis$geno$`1`$data[mapthis$pheno[1]==358,136]
mapthis$geno$`1`$data[mapthis$pheno[1]==358,137]

##### Checking locations across Chromosome 1 for individuals that deviate in Growth 1

gbsa_cms <- unname(mapthis$geno$`1`$map)
gbsa_locations<-trunc(as.numeric(gsub("CM000093.4.","",names(mapthis$geno$`1`$map))))
rqtlprob_gbsa<-mapthis$geno$`1`$prob[,,1]-mapthis$geno$`1`$prob[,,3]

snpms_locations<-trunc(as.numeric(LocMrksF2$GalGal5Mb[which(LocMrksF2$Chromosome==1)]))
length(snpms_locations)
hkqtlprob_snpms <- Geno_a[]
which(rownames(Geno_a) %in% mapthis$pheno[1]))

dim(rqtlprob_gbsa)[2]

###### Genetic and physical locations for QTL ######
# The genetic map used in Jacobsson, 2005 is published in the following paper:
# Jacobsson, L., Park, H.-B., Wahlberg, P., Jiang, S., Siegel, P. B., & Andersson, L. (2004). Assignment of fourteen microsatellite markers to the chicken linkage map. Poultry Science, 83(11), 1825–1831.
# Here, we use the marker locations in the genetic map of Wahlberg, 2009 to assign the physical locations
# of the mapped QTL in the GalGal assembly
#
# First, read in the information on the QTL mapped in Jacobsson, 2005 together
# with the markers-based confidence intervals. 

QTLJacobsson2005 <- read.table(file="~/Dropbox/Projects/Ongoing/HL_F2_SeqMrk/F2_re_seq/data/F2/150203_OCG_QTLJacobsson2005.txt",header=TRUE,sep="\t")
QTLJacobsson2005[,4]<-as.character(QTLJacobsson2005[,4])
QTLJacobsson2005[,4] <- gsub(' ', '', QTLJacobsson2005[,4])
QTLJacobsson2005[,5]<-as.character(QTLJacobsson2005[,5])
QTLJacobsson2005[,5] <- gsub(' ', '', QTLJacobsson2005[,5])

# The linkage map on Chromosome 7 is flipped - swap places for flanking markers
tmp <- QTLJacobsson2005[QTLJacobsson2005[,2]=="7",4]
QTLJacobsson2005[QTLJacobsson2005[,2]=="7",4] <- QTLJacobsson2005[QTLJacobsson2005[,2]=="7",5]
QTLJacobsson2005[QTLJacobsson2005[,2]=="7",5] <- tmp

# Then, read in the locations for the QTL reported in Wahlberg, 2009. In this paper, no confidence
# interval for the QTL are given, only the peak location

QTLWahlberg2009 <- read.table(file="~/Dropbox/Projects/Ongoing/HL_F2_SeqMrk/F2_re_seq/data/F2/150203_OCG_QTLWahlberg2009.txt",header=TRUE,sep="\t")
QTLWahlberg2009[,2]<-as.numeric(QTLWahlberg2009[,2])
QTLWahlberg2009[,3]<-as.numeric(QTLWahlberg2009[,3])
QTLWahlberg2009[,4]<-as.character(QTLWahlberg2009[,4])

###### Physical map-locations for QTL in Jacobsson, 2005 #########
# Determine Mb locations on GalGal5 for the confidence-interval markers QTL in Jacobsson, 2005.

# Two markers flanking QTL (MCW130 & BMP7) are missing in Wahlberg, 2009
# Instead use the neareast marker in Jacobsson, 2005 that is present in Wahlberg, 2009
# These are MCW0293 (MCW130) & MC3R (BMP7). We use the MC3R_PYRO from Wahlberg, 2009 for MC3R
QTLJacobsson2005[,5] <- gsub('MCW130', 'MCW0293', QTLJacobsson2005[,5])
QTLJacobsson2005[,5] <- gsub('BMP7', 'MC3R_PYRO', QTLJacobsson2005[,5])
for (i in 1:nrow(QTLJacobsson2005)){
  #Left marker
  QTLJacobsson2005[i,6] <- LocMrksF2[which((substr(LocMrksF2[,5],1,3)%in%substr(QTLJacobsson2005[i,4],1,3)) & (substr(LocMrksF2[,5],nchar(LocMrksF2[,5])-2,nchar(LocMrksF2[,5]))%in%substr(QTLJacobsson2005[i,4],nchar(QTLJacobsson2005[i,4])-2,nchar(QTLJacobsson2005[i,4])))),5]
  QTLJacobsson2005[i,7] <- LocMrksF2[which((substr(LocMrksF2[,5],1,3)%in%substr(QTLJacobsson2005[i,4],1,3)) & (substr(LocMrksF2[,5],nchar(LocMrksF2[,5])-2,nchar(LocMrksF2[,5]))%in%substr(QTLJacobsson2005[i,4],nchar(QTLJacobsson2005[i,4])-2,nchar(QTLJacobsson2005[i,4])))),3]  
  QTLJacobsson2005[i,8] <- LocMrksF2[which((substr(LocMrksF2[,5],1,3)%in%substr(QTLJacobsson2005[i,4],1,3)) & (substr(LocMrksF2[,5],nchar(LocMrksF2[,5])-2,nchar(LocMrksF2[,5]))%in%substr(QTLJacobsson2005[i,4],nchar(QTLJacobsson2005[i,4])-2,nchar(QTLJacobsson2005[i,4])))),2]  
  #Right marker
  QTLJacobsson2005[i,9] <- LocMrksF2[which((substr(LocMrksF2[,5],1,3)%in%substr(QTLJacobsson2005[i,5],1,3)) & (substr(LocMrksF2[,5],nchar(LocMrksF2[,5])-2,nchar(LocMrksF2[,5]))%in%substr(QTLJacobsson2005[i,5],nchar(QTLJacobsson2005[i,5])-2,nchar(QTLJacobsson2005[i,5])))),5]
  QTLJacobsson2005[i,10] <- LocMrksF2[which((substr(LocMrksF2[,5],1,3)%in%substr(QTLJacobsson2005[i,5],1,3)) & (substr(LocMrksF2[,5],nchar(LocMrksF2[,5])-2,nchar(LocMrksF2[,5]))%in%substr(QTLJacobsson2005[i,5],nchar(QTLJacobsson2005[i,5])-2,nchar(QTLJacobsson2005[i,5])))),3]  
  QTLJacobsson2005[i,11] <- LocMrksF2[which((substr(LocMrksF2[,5],1,3)%in%substr(QTLJacobsson2005[i,5],1,3)) & (substr(LocMrksF2[,5],nchar(LocMrksF2[,5])-2,nchar(LocMrksF2[,5]))%in%substr(QTLJacobsson2005[i,5],nchar(QTLJacobsson2005[i,5])-2,nchar(QTLJacobsson2005[i,5])))),2]  
}
QTLJacobsson2005<-QTLJacobsson2005[,c(1:4,7:8,5,10:11)]
colnames(QTLJacobsson2005)[c(2,5,6,8,9)]<-c("Chromosome","MbGG4L","cML","MbGG4R","cMR")
QTLJacobsson2005[,2]<-as.numeric(QTLJacobsson2005[,2])
QTLJacobsson2005[,5]<-as.numeric(QTLJacobsson2005[,5])
QTLJacobsson2005[,6]<-as.numeric(QTLJacobsson2005[,6])
QTLJacobsson2005[,8]<-as.numeric(QTLJacobsson2005[,8])
QTLJacobsson2005[,9]<-as.numeric(QTLJacobsson2005[,9])

###### Physical map-locations for QTL in Wahlberg, 2009 #########
# Identify the markers that flank the QTL-peak 
for (i in 1:nrow(QTLWahlberg2009)){
  #Left marker
  LChr <- (LocMrksF2[,1]==QTLWahlberg2009[i,2])
  LMrks <- (LocMrksF2[,2] <= QTLWahlberg2009[i,3])
  QTLWahlberg2009[i,5] <- LocMrksF2[which(LChr&LMrks)[length(which(LChr&LMrks))],5]
  QTLWahlberg2009[i,6] <- LocMrksF2[which(LChr&LMrks)[length(which(LChr&LMrks))],3]  
  QTLWahlberg2009[i,7] <- LocMrksF2[which(LChr&LMrks)[length(which(LChr&LMrks))],2]  
  #Right marker
  UMrks <- (LocMrksF2[,2] >= QTLWahlberg2009[i,3])
  QTLWahlberg2009[i,8] <- LocMrksF2[which(LChr&UMrks)[1],5]
  QTLWahlberg2009[i,9] <- LocMrksF2[which(LChr&UMrks)[1],3]  
  QTLWahlberg2009[i,10] <- LocMrksF2[which(LChr&UMrks)[1],2]
}
QTLWahlberg2009<-QTLWahlberg2009[,c(1:2,4:10)]
colnames(QTLWahlberg2009)[c(4,5,6,7,8,9)]<-c("MrkL","MbGG4L","cML","MrkR","MbGG4R","cMR")
QTLWahlberg2009[,2]<-as.numeric(QTLWahlberg2009[,2])
QTLWahlberg2009[,5]<-as.numeric(QTLWahlberg2009[,5])
QTLWahlberg2009[,6]<-as.numeric(QTLWahlberg2009[,6])
QTLWahlberg2009[,8]<-as.numeric(QTLWahlberg2009[,8])
QTLWahlberg2009[,9]<-as.numeric(QTLWahlberg2009[,9])

######### Overlap between significant F15 sweeps and F2 QTL ##########
# Read information on significant sweeps in F15 - these are obtained using script
# 150123_OCG_BackwardEliminationFinalProcedure.R

#load(file="~/Google Drive/projects/HL_AIL_F15/results/OCG_Final_Rev1/150209_OCG_BE_Bagging_FDR20.Rdata")
#load(file="~/Google Drive/projects/HL_AIL_F15/results/OCG_Final_Rev1/150209_OCG_BE_Bagging_FDR5.Rdata")
#load(file="~/Google Drive/projects/HL_AIL_F15/results/OCG_Final_Rev1/150209_OCG_BE_original_FDR5.Rdata")
#load(file="~/Google Drive/projects/HL_AIL_F15/results/OCG_Final_Rev1/150209_OCG_BE_original_FDR20.Rdata")
#load(file="~/Google Drive/projects/HL_AIL_F15/results/OCG_Final_Rev1/150209_OCG_BE_all_sweeps_F15")

SweepsF15 <- all_sweeps
SweepsF15 <- as.data.frame(SweepsF15)
rownames(SweepsF15) <- c(seq(1,nrow(SweepsF15),1))
SweepsF15[,1] <- as.character(SweepsF15[,1])
SweepsF15[,1] <- gsub('-', '_', SweepsF15[,1])
SweepsF15[,2] <- as.numeric(as.character(SweepsF15[,2]))
SweepsF15[,3] <- as.numeric(as.character(SweepsF15[,3]))

# Screen over F15-sweeps to identify overlaps with QTL
######## Sweep-overlaps with QTL from Jacobsson, 2005 #######
# Overlaps not only to BW56 QTL, but for all traits reported in the study
SweepsF15[,c(4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)]<-NA
k <- nrow(SweepsF15)
rowsSweepsF15 <- nrow(SweepsF15)
for (i in 1:rowsSweepsF15){
  tst1 <- SweepsF15[i,2]==QTLJacobsson2005[,2]
  tst2 <- SweepsF15[i,3]>QTLJacobsson2005[,5]
  tst3 <- SweepsF15[i,3]<QTLJacobsson2005[,8]  
  if (length(which(tst1 & tst2 & tst3))>0){
    SweepsF15[i,4] <- as.character(QTLJacobsson2005[which(tst1 & tst2 & tst3)[1],1])
    SweepsF15[i,20] <- as.character(QTLJacobsson2005[which(tst1 & tst2 & tst3)[1],3])
    SweepsF15[i,5:10] <- QTLJacobsson2005[which(tst1 & tst2 & tst3)[1],4:9]
    if (length(which(tst1 & tst2 & tst3))>1){
      for (j in (2:length(which(tst1 & tst2 & tst3)))) {
        k <- k+1
        SweepsF15[k,1:3] <- SweepsF15[i,1:3]
        SweepsF15[k,4] <- as.character(QTLJacobsson2005[which(tst1 & tst2 & tst3)[1],1])
        SweepsF15[k,5:10] <- QTLJacobsson2005[which(tst1 & tst2 & tst3)[j],4:9]
        SweepsF15[k,20] <- as.character(QTLJacobsson2005[which(tst1 & tst2 & tst3)[j],3])  
      }
    }
  }
}

######## Sweep-overlaps with QTL from Wahlberg, 2009 #######
# Overlaps not only to BW56 QTL, but for all traits reported in the study
# Wahlberg, 2009 gives no CI for the QTL
# Here, we use a threshold of +/- 5 Mb from the middle of the marker-interval where the
# QTL was mapped to determine overlaps

SweepsF15[,c(11,12,13,14,15,16,17)]<-NA
for (i in 1:rowsSweepsF15){
  midpoint <- (QTLWahlberg2009[,8]+QTLWahlberg2009[,5])/2
  QTL_L <- midpoint - 5
  QTL_R <- midpoint + 5
  tst1 <- SweepsF15[i,2]==QTLWahlberg2009[,2]
  tst2 <- SweepsF15[i,3]>QTL_L
  tst3 <- SweepsF15[i,3]<QTL_R
  if (length(which(tst1 & tst2 & tst3))>0){
    SweepsF15[i,11] <- as.character(QTLWahlberg2009[which(tst1 & tst2 & tst3)[1],1])
    SweepsF15[i,12:17] <- QTLWahlberg2009[which(tst1 & tst2 & tst3)[1],4:9] #If multiple matches, take the first                  
    SweepsF15[i,20] <- as.character(QTLWahlberg2009[which(tst1 & tst2 & tst3)[1],3])
    if (length(which(tst1 & tst2 & tst3))>1) {
      for (j in (2:length(which(tst1 & tst2 & tst3)))) {
        k <- k+1
        SweepsF15[k,1:3] <- SweepsF15[i,1:3]      
        SweepsF15[k,11] <- as.character(QTLWahlberg2009[which(tst1 & tst2 & tst3)[1],1])
        SweepsF15[k,12:17] <- QTLWahlberg2009[which(tst1 & tst2 & tst3)[j],4:9] #If multiple matches, take the first                    
        SweepsF15[k,20] <- as.character(QTLWahlberg2009[which(tst1 & tst2 & tst3)[j],3])  
      }
    }
  }
}

# Identify in which regions the selected sweep-markers are located
for (i in 1:nrow(SweepsF15)){
  SweepsF15[i,18]<-MrksF15[which(MrksF15[,4]%in%SweepsF15[i,1]),5]
}
colnames(SweepsF15)[c(4,11,18)]<-c("QTL2005","QTL2009","Region")

####### Output of overlap F15 sweeps and QTL #########
# First sort the dataframe and then output R-object needed to produce Table 2
SweepsF15 <- SweepsF15[order(SweepsF15[,2],SweepsF15[,3],SweepsF15[,4]),]
#save(SweepsF15, file="~/Google Drive/projects/HL_AIL_F15/results/OCG_Final_Rev1/150210_OCG_SweepsQTL.Rdata")

####### Single-locus F2 analysis ##########
####### Estimate effects at locations of F15 sweeps in F2 population ######
# First, identify markers in F2 within a defined physical distance from F15 Sweep ######
# Here, set to be sweepdist = 5 Mb
sweepdist <- 5 # Distance from peak F15 sweep marker; change to use different interval
F2SweepMrknames <- c()
F2SweepMrkNo <- c()
F15QTLSweeps <- unique(SweepsF15[,18])
for (j in 1:length(F15QTLSweeps)){
  i <- which(SweepsF15[,18]%in%F15QTLSweeps[j])[1]
  F2SweepMrks <- MrksF2[((MrksF2[,1]==SweepsF15[i,2]) & (MrksF2[,2] > (SweepsF15[i,3]-sweepdist)) & (MrksF2[,2] < (SweepsF15[i,3]+sweepdist))),]
  F2SweepMrknames <- c(F2SweepMrknames,as.character(t(F2SweepMrks[!is.na(F2SweepMrks[,1]),4])))
  F2SweepMrkNo <- c(F2SweepMrkNo,length(t(F2SweepMrks[!is.na(F2SweepMrks[,1]),4])))
}

# Estimate genetic effects of all F2 markers in the interval around the F15 sweep as calculated above
F2SweepMrkEstimates <- c()
for (i in 1:length(F2SweepMrknames)){
  xnam_F2 <- paste("Geno_a$", F2SweepMrknames[i], sep="")
  fmla_F2 <- as.formula(paste("Phenotypes$BW8 ~ Fixed$SEX +", paste(xnam_F2, collapse= "+",sep="")))
  F2.lm <- lm(fmla_F2, y=TRUE)
  F2SweepMrkEstimates <- c(F2SweepMrkEstimates,F2.lm$coefficients[3])
}

# Find the marker in each interval with the largest additive effect in the F2 population
F2MaxMrks <- c()
F2MaxMrks <- c(F2MaxMrks,F2SweepMrkEstimates[which(abs(F2SweepMrkEstimates[1:F2SweepMrkNo[1]]) == max(abs(F2SweepMrkEstimates[1:F2SweepMrkNo[1]])))])

for (i in 2:length(F2SweepMrkNo)){
  abs_a <- abs(F2SweepMrkEstimates[(sum(F2SweepMrkNo[1:(i-1)])+1):sum(F2SweepMrkNo[1:i])])
  max_a <- max(abs_a)
  F2MaxMrks <- c(F2MaxMrks, abs_a[which(abs_a == max_a)])
}

# The microsatellite marker MCW0234 is removed as it has identical genotypes to the SNP rs13796585 in this data
# If not removed, the model-matrix becomes singular
F2MaxMrks <- F2MaxMrks[names(F2MaxMrks)!="Geno_a$MCW0234"]

######## Multi-locus F2 analysis ########
# The estimates of the genetic effects for the sweeps in the F15 is done using a multi-locus model
# Therefore, a similar multilocusmodel is fitted including one marker from each F15 determined sweep in the F2 data

# Identify markers that are within mrkdist=x Mb of marker in model. Mrks in this interval
# are not to be included simultaneously in the model due to potential co-linearities
# in the statistical model. Here set mrkdist conservatively to 10 Mb

mrkdist<-10
for (i in 1:nrow(MrksF2)){
  if (any(MrksF2[i,8]%in%names(F2MaxMrks))){
    MrksF2[i,9]<-which(names(F2MaxMrks)%in%MrksF2[i,8])
    k<-i-1
    # Only do this for markers that have not already been decided to be downstream of sweep
    if (k>0) {
      while ((MrksF2[i,1]==MrksF2[k,1])&((MrksF2[i,7]-MrksF2[k,7])<mrkdist)) {
        MrksF2[k,9]<-MrksF2[i,9]
        k<-k-1
      }
    }
    k<-i+1
    if (k<(nrow(MrksF2)+1)) {
      while ((MrksF2[i,1]==MrksF2[k,1])&((MrksF2[k,7]-MrksF2[i,7])<mrkdist)) {
        MrksF2[k,9]<-MrksF2[i,9]
        k<-k+1
      }
    }
  } else{
    if (is.na(MrksF2[i,9])) {
      MrksF2[i,9]<-0      
    }
  }
}

# Identify the F2 markers that are selected to represent the sweeps
SweepmrksF2 <-matrix(data=NA,nrow=length(F2MaxMrks),ncol=2)
SweepmrksF2[,1]<-which(MrksF2[,8]%in%names(F2MaxMrks)) #Marker number
SweepmrksF2[,2]<-MrksF2[which(MrksF2[,8]%in%names(F2MaxMrks)),9] #Sweep ID

# If the markers selected to represent the sweeps are closer than (mrkdist) Mb as defined above, 
# keep only the one with the largest estimate as co-factor
Remove_mrks <-c()
if (length(unique(MrksF2[SweepmrksF2[,1],9])) < length(F2MaxMrks)) {
  for (i in unique(MrksF2[SweepmrksF2[,1],9])){
    if (length(which(MrksF2[SweepmrksF2[,1],9]%in%i))>1){
      same_sweep <- which(names(F2MaxMrks)%in%paste("Geno_a$",MrksF2[SweepmrksF2[which(MrksF2[SweepmrksF2[,1],9]%in%i),1],4],sep=""))
      Remove_mrks <- c(Remove_mrks,MrksF2[SweepmrksF2[same_sweep[which(abs(as.numeric(F2MaxMrks[same_sweep]))!=max(abs(as.numeric(F2MaxMrks[same_sweep]))))],1],4])
    }
  }
}
F2MaxMrks <- F2MaxMrks[!names(F2MaxMrks)%in%paste("Geno_a$",Remove_mrks,sep="")]

# First, estimate the effects of the F2-sweep markers in the multi-locus model
# Remove MCW0234 marker as its its genotype overlaps with other in same sweep
fmla_F2_all <- as.formula(paste("Phenotypes$BW8 ~ Fixed$SEX +", paste(names(F2MaxMrks),"-Geno_a$MCW0234", collapse= "+",sep="")))
F2_all.lm <- lm(fmla_F2_all, y=TRUE)
summary(F2_all.lm)

# Then, scan across all genotyped markers in the F2 by regressing on the QTL genotype probabilities
# at the location of the marker. The sweep-markers, that are not part of the tested sweep, 
# are included as co-factors in the model

for (i in 1:nrow(MrksF2)){
  # If the marker is included in the full model or near the marker that is
  if (MrksF2[i,9]>0) {
    # If it is the marker that is in the model - just use full model
    if (MrksF2[i,8]%in%names(F2MaxMrks)){
        fmla_F2scan <- as.formula(paste("Phenotypes$BW8 ~ Fixed$SEX +", paste(names(F2MaxMrks), collapse= "+",sep="")))      
        F2_scan.lm <- lm(fmla_F2scan, y=TRUE)      
      # If it is not the actual sweep marker, replace it with the nearby marker:
      } else {
        #Identify which marker to be replaced & replace it in model
        replmrk<-c()
        replmrk<-SweepmrksF2[which(SweepmrksF2[,2]%in%MrksF2[i,9]),1]
        fmla_F2scan <- as.formula(paste("Phenotypes$BW8 ~ Fixed$SEX +",MrksF2[i,8],"+", paste(names(F2MaxMrks[which(!names(F2MaxMrks)%in%MrksF2[replmrk,8])]), collapse= "+",sep="")))
        F2_scan.lm <- lm(fmla_F2scan, y=TRUE)
      }
    }
  # If the marker is not in, or near, the sweep marker - just add it to the model
  else{
    fmla_F2scan <- as.formula(paste("Phenotypes$BW8 ~ Fixed$SEX +",MrksF2[i,8],"+", paste(names(F2MaxMrks), collapse= "+",sep="")))      
    F2_scan.lm <- lm(fmla_F2scan, y=TRUE)      
  }
  m<-summary(F2_scan.lm)
  est<-m$coefficients[which(rownames(m$coefficients)==MrksF2[i,8]),1]
  SE<-m$coefficients[which(rownames(m$coefficients)==MrksF2[i,8]),2]
  sign<-m$coefficients[which(rownames(m$coefficients)==MrksF2[i,8]),4]  
  MrksF2[i,10]<-est
  MrksF2[i,11]<-SE
  MrksF2[i,12]<-sign
}
colnames(MrksF2)[9:12]<-c("Sweep","MeanMM","SEMM","p-valueMM")

# Then also perform a simple single marker scan without co-factors in the model
for (i in 1:nrow(MrksF2)){
  fmla_F2scan <- as.formula(paste("Phenotypes$BW8 ~ Fixed$SEX +",MrksF2[i,8]))      
  F2_scan.lm <- lm(fmla_F2scan, y=TRUE)      
  m<-summary(F2_scan.lm)
  est<-m$coefficients[which(rownames(m$coefficients)==MrksF2[i,8]),1]
  SE<-m$coefficients[which(rownames(m$coefficients)==MrksF2[i,8]),2]
  sign<-m$coefficients[which(rownames(m$coefficients)==MrksF2[i,8]),4]  
  MrksF2[i,13]<-est
  MrksF2[i,14]<-SE
  MrksF2[i,15]<-sign
}
colnames(MrksF2)[13:15]<-c("MeanSM","SESM","p-valueSM")

######## Multi-locus sweep-scan in the F15 population ######
# Here, estimate the effects of all sweeps in the F15 population using the
# significant sweeps as co-factors in the analysis

SweepsF15[,19] <- paste("geno_a$", SweepsF15[,1], sep="")  
MrksF15[,7] <- paste("geno_a$", MrksF15[,4], sep="")

# Different phenotypes
#phe<-as.numeric(fix_data$pheno$BW4)
phe<-as.numeric(fix_data$pheno$BW8)
min.lm <- lm(phe ~ fx1)
phe<-rstandard(min.lm)^2



# Do not fit markers that are closer than (too_close) Mb; here set to 2 Mb
too_close <- 2
for (i in 1:nrow(MrksF15)){
  MrksF15[i,]
  # Check if any marker in full model within (too_close) Mb of the tested marker
  if (any((SweepsF15[which(SweepsF15[,2]%in%MrksF15[i,2]),3]-MrksF15[i,3])<too_close)) {
  # If it is in the same sweep
    if (MrksF15[i,5]%in%SweepsF15[,18]){
      if (SweepsF15[which(SweepsF15[,18]%in%MrksF15[i,5]),1]==MrksF15[i,4]){
        fmla_full <- as.formula(paste("phe ~ fx1 +", paste(SweepsF15[,19], collapse= "+",sep="")))
        F15_sweep.lm <- lm(fmla_full)
      } else {
        fmla_full_1rep <- as.formula(paste("phe ~ fx1 +",MrksF15[i,7] ,"+", paste(SweepsF15[which(!SweepsF15[,18]%in%MrksF15[i,5]),19], collapse= "+",sep="")))
        F15_sweep.lm <- lm(fmla_full_1rep, y=TRUE)
      }
    } 
    # If it is in another sweep on the chromosome
    else {
      too_near <- which((SweepsF15[which(SweepsF15[,2]%in%MrksF15[i,2]),3]-MrksF15[i,3])<too_close)
      ToRemove <- SweepsF15[which(SweepsF15[,2]%in%MrksF15[i,2])[too_near],19]
      fmla_full_1rep <- as.formula(paste("phe ~ fx1 +",MrksF15[i,7] ,"+", paste(SweepsF15[which(SweepsF15[,19]!=ToRemove),19], collapse= "+",sep="")))
      F15_sweep.lm <- lm(fmla_full_1rep, y=TRUE)    
    }  
  }
  else {
    fmla_full <- as.formula(paste("phe ~ fx1 +",MrksF15[i,7],"+", paste(SweepsF15[,19], collapse= "+",sep="")))
    F15_sweep.lm <- lm(fmla_full)  
  }
  m<-summary(F15_sweep.lm)
  est<-m$coefficients[which(rownames(m$coefficients)==paste("geno_a$",MrksF15[i,4],sep="")),1]
  SE<-m$coefficients[which(rownames(m$coefficients)==paste("geno_a$",MrksF15[i,4],sep="")),2]
  sign<-m$coefficients[which(rownames(m$coefficients)==paste("geno_a$",MrksF15[i,4],sep="")),4]  
  MrksF15[i,8]<-est
  MrksF15[i,9]<-SE
  MrksF15[i,10]<-sign
}
colnames(MrksF15)[7:10]<-c("Model","Mean","SE","p-value")

vGWASmrks <- which(MrksF15[,10]<0.12)

########## Collect the results from the F2 and F15 analyses ############
# Results are stored in a dataframe *Results* and then saved to file for plotting

F15_SweepScan <- MrksF15[,c(2,3,4,8,9,10)]
F2_Scan <- MrksF2[,c(1,7,4,10,11,12,13,14,15)]
head(F15_SweepScan)
head(F2_Scan)
nrow(F15_SweepScan)
nrow(F2_Scan)
Results <- matrix(data=NA,nrow=(nrow(F15_SweepScan)+nrow(F2_Scan)),ncol=14)
Results <- as.data.frame(Results)
Results[1:nrow(F2_Scan),1:9]<-F2_Scan[,1:9]  
Results[(nrow(F2_Scan)+1):nrow(Results),1:3]<-F15_SweepScan[,1:3]
Results[(nrow(F2_Scan)+1):nrow(Results),10:12]<-F15_SweepScan[,4:6]
Results[1:nrow(F2_Scan),13]<-LocMrksF2[,2] #cM F2
Results[,1]<-as.numeric(Results[,1])
Results[,2]<-as.numeric(Results[,2])
Results[,4]<-as.numeric(Results[,4])
Results[,5]<-as.numeric(Results[,5])
Results[,6]<-as.numeric(Results[,6])
Results[,7]<-as.numeric(Results[,7])
Results[,8]<-as.numeric(Results[,8])
Results[,9]<-as.numeric(Results[,9])
Results[,10]<-as.numeric(Results[,10])
Results[,11]<-as.numeric(Results[,11])
Results[,12]<-as.numeric(Results[,12])
Results[,13]<-as.numeric(Results[,13])
Results[,14]<-Results[,13]

# Sort the Results dataframe by Chromosome and Physical location (Mb) in the GalGal4 assembly
Results <- Results[order(Results[,1],Results[,2]),]
colnames(Results)[1:14]<-c("Chromosome","GG4Mb","Marker","MeanF2MM","SEF2MM","P-valueF2MM","MeanF2SM","SEF2SM","P-valueF2SM","MeanF15MM","SEF15MM","p-valueF15MM","cM","cMIP")

# Physical location for marker rs13834451 does not fit the linkage map
# Set its location as missing here
Results[which(Results[,3]=="rs13834451"),13]<-NA

#Estimate location on Linkage-map for F15 markers via linear interpolation
for (i in which(is.na(Results[,13]))){
  #Check for markers flanking segments to be interpolated
  k<-i-1
  while (is.na(Results[k,13])){
    k<-k-1      
  }
  flankingL<-k
  k<-i+1
  while (is.na(Results[k,13])){
    k<-k+1      
  }
  flankingU<-k  
  intMb <- (Results[i,2]-Results[flankingU,2])/(Results[flankingU,2]-Results[flankingL,2])
  distcM <- Results[flankingU,13]-Results[flankingL,13]
  Results[i,14] <- Results[flankingU,13]+distcM*intMb
}

# Calculate the distribution of Sweep-lengths in cM (F2) & coverage of genome

Results2 <- Results[!is.na(Results$MeanF15MM),]
for (i in 1:nrow(Results2)) {
  Results2[i,15] <- regs$Reg_new[regs[,3]==Results2[i,3]]  
}
Results2[,15]
maxcM <- do.call(rbind,lapply(split(Results2,Results2$V15),function(chunk) chunk[which.max(chunk$cMIP),]))
mincM <- do.call(rbind,lapply(split(Results2,Results2$V15),function(chunk) chunk[which.min(chunk$cMIP),]))
lengths <-(maxcM[,14]-mincM[,14])
F2lengths <- 0.5*(1-(exp(1)^(-2*(lengths/100))))
mean(lengths)
range(lengths)
sum(lengths)
lengths2 <- lengths[lengths > 0]
length(lengths2)
mean(lengths2)
range(lengths2)
sum(lengths2)

regs<-read.delim("~/Google Drive/projects/HL_AIL_F15/data/141217_ZS_cluster_with_effects.txt",sep="\t",stringsAsFactors=F)
F15Gmap <- regs[,c(3,1,7)]
F15maxcM <- do.call(rbind,lapply(split(F15Gmap,F15Gmap$Reg_new),function(chunk) chunk[which.max(chunk$G_map),]))
F15mincM <- do.call(rbind,lapply(split(F15Gmap,F15Gmap$Reg_new),function(chunk) chunk[which.min(chunk$G_map),]))
F15lengths <-(F15maxcM[,3]-F15mincM[,3])
F15lengths <- 0.5*(1-(exp(1)^(-2*(F15lengths/100))))
mean(F15lengths)
range(F15lengths)
sum(F15lengths)
F15lengths2 <- F15lengths[F15lengths > 0]
length(F15lengths2)
mean(F15lengths2)
range(F15lengths2)
sum(F15lengths2)

maxMb <- do.call(rbind,lapply(split(Results2,Results2$V15),function(chunk) chunk[which.max(chunk$GG4Mb),]))
minMb <- do.call(rbind,lapply(split(Results2,Results2$V15),function(chunk) chunk[which.min(chunk$GG4Mb),]))
Mblengths <-(maxMb[,2]-minMb[,2])
mean(Mblengths)
range(Mblengths)
sum(Mblengths)

# cM/Mb
F15cMMb <- F15lengths/Mblengths
F15cMMb <- F15cMMb[F15cMMb!="NaN"]
mean(F15cMMb)
hist(F15cMMb)

F2F15mapinflation <- F15lengths/F2lengths
F2F15mapinflation<-F2F15mapinflation[F2F15mapinflation!="NaN"]
mean(F2F15mapinflation)
range(F2F15mapinflation)
hist(F2F15mapinflation)

r=0.05*6.5
-log(1-2*r)/2

#Write to file for test-plotting to visualize overlap between F2 and F15 results
#save(Results,file="~/Google Drive/projects/HL_AIL_F15/results//OCG_Final_Rev1/150210_OCG_Additive_effects_markers.Rdata")
Resultsmat <- as.matrix(Results)
write(t(Resultsmat), file = "~/Google Drive/projects/HL_AIL_F15/results/OCG_150206/150206_OCG_ScanF2F15.txt", ncolumns=ncol(Resultsmat),sep = "\t")