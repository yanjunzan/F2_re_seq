#########################################################################################################
# Script for building a linkage map in the HWS x LWS F2 using GBSA scored genotypes
#
# Author Örjan Carlborg 180201
#########################################################################################################

library(qtl)
setwd(dir = '~/Dropbox/Projects/Ongoing/HL_F2_SeqMrk/F2_re_seq/')

############# 1. Read in the data from file & check it ############
## File with the GBSA type I genotypes, i.e. GBSA scored markers with alleles completely fixed in founder lines
#f2cross <- read.cross(format="csv",dir="/Users/orjanc/Dropbox/Projects/Ongoing/HL_F2_SeqMrk/F2_re_seq/data/",file = "F2_n707Rqtl_ID.csv",na.strings = "-",genotypes = c("A", "H", "B", "C", "D"), estimate.map = FALSE, map.function = "haldane", sep = ";")

## File with the GBSA type II genotypes, i.e. GBSA scored markers with alleles fixed between founders within F2 family
f2cross <- read.cross(format="csv",dir="/Users/orjanc/Dropbox/Projects/Ongoing/HL_F2_SeqMrk/F2_re_seq/data/",file = "F2_827.within.fam.mat.bin50.cut.bugcorrected.180223.csv",na.strings = "NA",genotypes = c("A", "H", "B", "C", "D"), estimate.map = FALSE, map.function = "haldane", sep = ";")
############# 2. Have a look at the imported data #############

## First an overall summary
summary(f2cross)

## Generate a plot to visualize missing genotypes across the markers 
plotMissing(f2cross)

## Visual check for missing genotypes per individual & marker
par(mfrow=c(1,2), las=1)
plot(ntyped(f2cross), ylab="No. typed markers", main="No. genotypes by individual") 
plot(ntyped(f2cross, "mar"), ylab="No. typed individuals",main="No. genotypes by marker")

############# 3. Do a first cleaning of the data ##########################

######## 3.A. First drop the markers that with genotype scores in few individuals as they complicates map estimation #######
## Here use a threshold of 50 markers based on the visual inspection aboveas this is enough to 
## estimate recombination with "reasonable" accuracy

nt.bymar <- ntyped(f2cross, "mar")
todrop <- names(nt.bymar[nt.bymar < 50])
mapthis <- drop.markers(f2cross, todrop)

#### Double-check that the markers have been dropped
par(mfrow=c(1,2), las=1)
plot(ntyped(mapthis), ylab="No. typed markers", main="No. genotypes by individual") 
plot(ntyped(mapthis, "mar"), ylab="No. typed individuals",main="No. genotypes by marker")

############# 3.B. Identify individuals that have GBSA scores for few markers ################
#### Here, this is an indication that these individuals had low marker coverage
#### The individual-ID's are given in the phenotype-column of the F2-cross infile

id_lt200mrk <- mapthis$pheno[which(ntyped(mapthis)<300),1] 

## Check the sequence-coverage for the selected individuals
#### First read in the sequence coverage file
coverage_F2 <- read.table(file = './data/180208_F2_coverage.txt', row.names = 1)
#### Print out the coverages
coverage_F2[row.names(coverage_F2) %in% id_lt200mrk,1]
#### As expected, the scoring of few markers was due to the low coverage for these individuals.

## The individuals with few genotyped markers are dropped from the linkage map construction
mapthis <- subset(mapthis, ind=(ntyped(mapthis)>200))
#### Note that this does not indicate that the genotypes are of low quality - can thus be included in later mapping!

#### Check dataset again after the low-coverage individuals have been dropped
par(mfrow=c(1,2), las=1)
plot(ntyped(mapthis), ylab="No. typed markers", main="No. genotypes by individual") 
plot(ntyped(mapthis, "mar"), ylab="No. typed individuals",main="No. genotypes by marker")

############ 3.C. Check whether any individuals are duplicates that needs to be removed #############

## Visually inspect the relationships between the individuals in the population
cg <- comparegeno(mapthis)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes") 
rug(cg[lower.tri(cg)])

##  Identify the individuals that are highly related 
#### Here, a threshold of 0.9 is used
wh <- which(cg > 0.9, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]
wh

## Inspect their genotype similarity
#### This is done by checking the pairwise genotype tables; no relatedness higher so below for illustration only
#g <- pull.geno(mapthis)
#table(g[223,], g[238,])

########### 3.D. Check whether any markers are duplicates ############
## No duplicates are expected given the way the genotypes are scored & the size of the population
print(dup <- findDupMarkers(mapthis, exact.only=FALSE))
# Duplicates here due to using markers with few genotyped individuals - removed if threshold is higher

########### 3.E. Check for distorted segregation patterns ###########
## This is an important quality control for this type of markers
## This as it could indicate an error in GBSA genotype assignment despite most individuals having reads for them
gt <- geno.table(mapthis)
## A Bonferroni corrected threshold for the number of markers is used to identify and remove markers
gt[gt$P.value < 0.05/totmar(mapthis),]
todrop <- rownames(gt[gt$P.value < 0.05/totmar(mapthis),])
mapthis <- drop.markers(mapthis, todrop)
## Note1: More stringent filtering could, however, be used to obtain an even higher quality dataset for map construction
## Note2: These markers do not necessarily have to be removed for mapping as a shortage of heterozygotes is expected for markers/bins with few reads

########### 3.F. Evaluate the recombination frequencies for the markers ##########
## This is a way to identify individuals that have skewed genotypes
#### Could be an indication of DNA contamination / pedigree errors
g <- pull.geno(mapthis)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq))
par(mfrow=c(1,3), las=1)
for(i in 1:3) plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "AB", "BB")[i], ylim=c(0,1))

## Drop these individuals due to the skewed genotype frequencies
#### Individuals to keep
id_gtfreq <- which(gfreq[1,]<0.5 & gfreq[2,]<0.65 & gfreq[2,]>0.35 & gfreq[3,]<0.5)
#### Drop
mapthis <- subset(mapthis, ind=(id_gtfreq))

########## 4. Estimate a first rough draft of the linkage-map ############

####  4.A. Estimate the recombination-fractions between pairs of markers #####
mapthis <- est.rf(mapthis)
mapthis_tst <- markerlrt(mapthis)

## Plot the obtained recombination fractions against LOD-scores
#### Do not want any high recombination-fractions together with high LOD-scores
rf <- pull.rf(mapthis)
lod <- pull.rf(mapthis, what="lod")
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")

#### Or weird linkages indicating that markers are wrongly ordered in the input data
#### such as between chromosomes or orders within chromsomes

## Evaluate by plotting pairwise recombination-fractions and LOD scores across the genome
dev.off()
plotRF(mapthis, alternate.chrid=TRUE)

## Looks pretty OK from the graph - possibly a couple of linkage-groups that could be explored further 
#### later to see if there might be some structural rearrangements on some segregating haplotypes
#### or problematic genotypes in some areas

######## 4.B. Get the complete genetic linkage map  #########
genmap <- est.map(mapthis, error.prob=0.005)
## As numerical summary
summaryMap(genmap)
## As a plot
plotMap(genmap)
mapthis <- replace.map(genmap, newmap)

########### 5. Improve quality of linkage map ###########
## The raw linkage map is too long, suggesting that there might be too high recombination-rates scored in some individuals
## Therefore try to identify these

########### 5.A. Find & remove individuals with excessive recombination
## These could be ones with problems in assigning GBSA scores or sample mixups
plot(countXO(mapthis), ylab="Number of crossovers")
hist(countXO(mapthis), ylab="Number of crossovers")
mean(countXO(mapthis), ylab="Number of crossovers")

## Bi-modal distribution of crossovers identifies individuals with approximately double the
## amount of recominations than the majority. Mixup of F2 samples or bug in code??

## Remove these problematic individuals
#### Threshold of 70 recombinations picked arbitrarily from visual inspection.
mapthis2 <- subset(mapthis, ind=(countXO(mapthis) < 70))
rf <- pull.rf(mapthis2)
lod <- pull.rf(mapthis2, what="lod")
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")

plot(countXO(mapthis2), ylab="Number of crossovers")
hist(countXO(mapthis2), ylab="Number of crossovers")
mean(countXO(mapthis2), ylab="Number of crossovers")

#########  6. Create final map with clean dataset ###############
## Build map
newmap2 <- est.map(mapthis2, error.prob=0.05)
summaryMap(newmap2)
dev.off()
plotMap(newmap2)

#########  7. Output final map in text and graphics ##############
mapthis <- replace.map(mapthis, newmap2)
summaryMap(mapthis)
dev.off()
plotMap(mapthis)

######### 8. Save the R/qtl object with linkage map ##############
save(mapthis, file = '~/Dropbox/Projects/Ongoing/HL_F2_SeqMrk/F2_re_seq/data/HL_F2_mapthis.Rdata')

#### QTL mapping ####

# Phenotype checking
plotPheno(mapthis, pheno.col=2)
plotMissing(mapthis, reorder=TRUE)

# Calculate QTL genotype probabilities
mapthis <- calc.genoprob(mapthis, step=1)
mapthis <- calc.genoprob(mapthis, step=0, off.end = 0)

# ML QTL scan
out.em <- scanone(mapthis, pheno.col = 2)
summary(out.em)
plot(out.em)

## covariates
ac <- pull.pheno(mapthis, c("SEX"))

## scan with additive but not the interactive covariate
out.acovar <- scanone(mapthis, pheno.col = 2, chr = c(1,3,4,7,20), addcovar=ac)
plot(out.acovar)
summary(out.acovar)
mname_c1 <- find.marker(mapthis,chr=1, pos=413)
mname_c4 <- find.marker(mapthis,chr=4, pos=15)
mname_c7 <- find.marker(mapthis,chr=7, pos=72)
effectplot(mapthis,pheno.col = 2, mname1 = mname_c1)
effectplot(mapthis,pheno.col = 2, mname1 = mname_c4)
effectplot(mapthis,pheno.col = 2, mname1 = mname_c7)

# take out several QTLs and make QTL object
qc <- c(1, 4, 7)
qp <- c(413, 15, 72)
fake.f2 <- subset(mapthis, chr=qc)

fake.f2 <- calc.genoprob(fake.f2, step=2, err=0.001)
qtl <- makeqtl(fake.f2, qc, qp, what="prob")

Sex <- data.frame(Sex=pull.pheno(fake.f2, "SEX"))
lod.add <- fitqtl(fake.f2, pheno.col=2, qtl, formula=y~Q1+Q2+Q3+Sex, cov=Sex, method="hk")
summary(lod.add)

# HK QTL scan
out.hk <- scanone(mapthis, pheno.col = 2, chr = c(1,3,4,7,20), method="hk", addcovar = ac)
plot(out.hk)

# Compare
plot(out.em, out.hk, out.acovar, col=c("blue", "red","green"))

#### Additional stuff ########

####### A.1. Test infering linkage groups from the data ########
##   Here this basically splits chromosomes with large distances between markers
lg <- formLinkageGroups(mapthis, max.rf=0.3, min.lod=6) 
table(lg[,2])

###### A.2. Try reordering markers where linkages extend over long regions - possible rearrangements? #######

#### For example chromosome 4 
mapthis <- orderMarkers(mapthis, chr=4)
pull.map(mapthis, chr=4)

###### A.3. Estimate maximum-likelihood estimate of the overall genotyping error from data ############
## Estimate likelihoods
loglik <- err <- c(0.001, 0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02) 
for(i in seq(along=err)) { 
  cat(i, "of", length(err), "\n")
  tempmap <- est.map(mapthis, error.prob=err[i])
  loglik[i] <- do.call(sum,sapply(tempmap, attr, "loglik")) 
}
lod <- (loglik - max(loglik))/log(10)

## Plot likelihoods to identify the most likely overall GBSA genotyping error rate
plot(err, lod, xlab="Genotyping error rate", xlim=c(0,0.02), ylab=expression(paste(log[10], " likelihood")))

###### A.4. Evaluate skewed genotype frequencies   ###############

## We note that some individuals have more genome-wide "AA/AB/BB" genotypes than expected
#### Which are these individuals?
id_gtfreqAA <- mapthis$pheno[which(gfreq[1,]>0.6),1]
id_gtfreqAB <- c(mapthis$pheno[which(gfreq[2,]>0.75),1],mapthis$pheno[which(gfreq[2,]<0.35),1])
id_gtfreqBB <- mapthis$pheno[which(gfreq[3,]>0.6),1]
id_gtfreqALL <- sort(unique(c(id_gtfreqAA,id_gtfreqAB,id_gtfreqBB)))

#### Is it in any way related to coverage? 
coverage_F2[row.names(coverage_F2) %in% id_gtfreqAA,1]
coverage_F2[row.names(coverage_F2) %in% id_gtfreqAB,1]
coverage_F2[row.names(coverage_F2) %in% id_gtfreqBB,1]
#### No, coverage was not low, so check them futher for other problems!

#### Genotype skew is for every chromosome, suggests sample mixup with parent
table(mapthis$geno$`1`$data[which(gfreq[3,]>0.6),])
table(mapthis$geno$`2`$data[which(gfreq[3,]>0.6),])
table(mapthis$geno$`3`$data[which(gfreq[3,]>0.6),])

###### A.5. Look for genotyping errors in individual data ############
## Useful for final quality fix of data
mapthis <- calc.errorlod(mapthis, error.prob=0.05)
print(toperr <- top.errorlod(mapthis, cutoff=5))

## Manually go through chromosomes
plotGeno(mapthis, chr=1, ind=toperr$id[toperr$chr==1], cutoff=5, include.xo=TRUE,cex = .5)
plotGeno(mapthis, chr=14, ind=toperr$id[toperr$chr==14], cutoff=4, include.xo=TRUE)

## Revisit segregation distortion across genome to view these for markers
gt <- geno.table(mapthis, scanone.output=TRUE)
par(mfrow=c(2,1))
plot(gt, ylab=expression(paste(-log[10], " P-value"))) > plot(gt, lod=3:5, ylab="Genotype frequency")
abline(h=c(0.25, 0.5), lty=2, col="gray")

## Identify marker that seem problematic at Bf-threshold = 2 x 10-4
gt[which(gt$neglog10P > -log10(0.05/269)),]

# SD appears due to i) large number of missing genotypes and ii) missing heterozygotes. 
# Possibly due to difficulty in being able to call heterozytotes?? None overlaps with ind of error ind genotypes

####### A.6. Checking of crossovers in individual data ###########

## Which individuals have more crossovers than expected?
id_ex_co <- mapthis$pheno[which(countXO(mapthis) > 70),1]
xo_pop <- cbind(mapthis$pheno[which(countXO(mapthis) > 0),1],countXO(mapthis)[which(countXO(mapthis) > 0)])

## Evaluate their sequence overage
coverage_F2[row.names(coverage_F2) %in% id_ex_co,1]
#### Have good overall coverage, hence likely not coverage related
