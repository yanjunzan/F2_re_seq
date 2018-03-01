#########################################################################################################
# Script for mapping in the HWS x LWS F2 using the GBSA scored genotypes
#
# Author Örjan Carlborg 180227
#########################################################################################################

library(qtl)
setwd(dir = '~/Dropbox/Projects/Ongoing/HL_F2_SeqMrk/F2_re_seq/')

############# 1. Load the R/qtl object with genotype, phenotype & linkage-map data ########
## This object is created using the script 180131_linkagemap.R
save(mapthis, file = '~/Dropbox/Projects/Ongoing/HL_F2_SeqMrk/F2_re_seq/data/HL_F2_mapthis.Rdata')

#### QTL mapping ####

# Phenotype checking
plotPheno(mapthis, pheno.col=2)
plotMissing(mapthis, reorder=TRUE)

# Calculate QTL genotype probabilities
# For all positions according to step-size
mapthis <- calc.genoprob(mapthis, step=1, off.end = 0)
# At marker-locations
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
