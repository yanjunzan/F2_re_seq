
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> AIL IMPUTE<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# From Yanjun 171205
# I documented what I have done in this README.R file
setwd("~/Dropbox/Projects/Ongoing/HL_F2_SeqMrk/From.Yanjun/")
#Yanjun path
setwd("~/Documents/impute/From.Yanjun/")

source("./bin/R/Functions.AIL.impute.R")


###>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>STEP I: SNP calling<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Process: 1.BAM files were generated before using BWA(Group information and duplication have been marked using picard) 
#            those files/file-path are in /tmp/nas/yanjun/AIL/founder.hap/ref.chr1.171030/171030.all.f2.bam.list.txt
#          2.Run R script below in diprotodon will get two vcf for one individual, one with XX_171101_all.vcf.recode.vcf 
#           have all the mrk called and another one with XX_171101_informative.vcf.recode.vcf have all fixed mrk called.
Rscript ./bin/R/171204.call.snp.R ## this script  take in the bam file and paralle them, the computing is done by caling
#                                 bash script "./bin/bash/171204_SNP.call.sh"


###<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<STEP II: Merge indiviual VCF<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# After inviduals are called vcfs need to be merged using bash scripts "./bin/bash/171204_merge.vcf.sh".
# The bash script first merge all individual vcf and then merge it with founder vcf, then scp back the data to local drive

###<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<Analysis in R <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
### A few R packages are needed in this section. data.table, zoo. but you might need more
### Python 3.4 is also required, python library re, sys, gzip are needed(In theory these should come with default installation)
# Prepare input: Two random F0(one from HWS and one from LWS) are selected to represent the founder haplotype at those fixed loci
F_id <- ".*1653.*"
M_id <- ".*2109.*"
#INPUT vcf file with all the founder and F2, need to be ziped and require a space before
vcf.file.f2 <- " ./data/171115_output.vcf.gz"

# The python script below reformat the massive vcf to only select individuals with two founder, need a space after
py.script <- "python3 ./bin/python/Read.vcf.py " #/usr/local/Cellar/python3/3.6.3/bin/

### Analysis start here
###
header.f2 <- vcf.header(vcf = vcf.file.f2) # check the vcf to getout the sample ID. there might be a broken pipe error, but ignore it
f2 <- header.f2[grep(pattern = "F2.*",x = header.f2)][1:30] ## will be all F2 in the end, for now a subset since all not available
#f2 <- f2[1] #If only one individual is desired

#Get chromosomes; currently only the macro + sex chromosomes; TO BE EXTENDED TO INCLUDE ALL CONTIGS
all.chr  <- read.table("~/Dropbox/Projects/Ongoing/HL_F2_SeqMrk/From.yanjun/data/chr_id.match.txt",stringsAsFactors = F,header = T,sep="\t")
chroms <- all.chr$INSDC

## First analysis will be based on marker genotypes estimated in bins of defined physical size
# Physically bin the genome into bins and calculate
bin.size <- 1e6 # Define bin size (bp)
cut.number <- 5 # Define smallest number of markers in bin to calculate average; otherwise get NA
png(filename = paste("./results/Phy.f2_10.sti.png",sep = ""),width = 1024,height = 1024 )

# Run this if you want to plot a certain chromosome and individual
chromosome <- 1 #Define chromosome to plot 
f2.plot <- c(10,11,12,13,14) #Define which f2's to plot

dev.off()
par(mfrow=c(length(f2.plot)+2,1),mar=c(1,1,1,1),xpd=T) #3 panels; updated from mar=c(2,3,1,3) due to plot error on MBP
  
#First plot a panel with locations of markers that are fixed between HWS and LWS on the chromosomes
file.in <- fread("./data/171113.Fixed.sites.founder.all.chr.txt") # read HWS vs LWS fixed mrks from file
pos.chr1 <- file.in$V2[file.in$V1 == chroms[chromosome]] 
plot(pos.chr1,y=rep(0.5,length(pos.chr1)),cex=0.8,pch=19)

#Then plot a panel with a histogram showing the number of markers in each of the bins
hist(pos.chr1,ceiling(max(pos.chr1)/bin.size))

#Then plot the average scores for the selected f2 on the chromosome
for( i in 1:length(f2.plot)){
  pathout <- paste(" ","./result/",f2[i],".vcf",sep="")
  data1 <- read_trio(of_id = f2[f2.plot[i]],F_id = F_id,M_id = M_id,vcf.file = vcf.file.f2,pathout = pathout,generate.input = F)
  setkey(data1,chr)
  data1.chr1 <- data1[chroms[chromosome]] # which chromsome to plot 
  trac <- tracing_physical(input = data.frame(data1.chr1),bin = bin.size,chr.len = max(pos.chr1),cut = cut.number) # trace inheritance in bin
  plot(x = trac$loca,trac$f,col="black",pch=19,ylim=c(0,1),frame.plot = F)
  points(x=trac$loca,y = trac$m,col="red",pch=19)
  #  abline(h=c(0.8,0.9),lty="dashed",col="red")
}

# Create data-structure (file) with scored genotypes across all individuals
k<-chroms[34]
for (k in chroms ){
  # get info on analyzed chromosome
  file.in <- fread("./data/171113.Fixed.sites.founder.all.chr.txt") # this is a file with all the founder fixed mrk
  pos.chr1 <- file.in$V2[file.in$V1 == k] 
  num.bin <- 1+max(pos.chr1)/bin.size
  bin.loca <- seq(0,max(pos.chr1),by=bin.size)+1
#now test with one matrix / chromosome
  genotypes <- matrix(nrow = length(f2),ncol = length(bin.loca))
  nmarkers <- matrix(nrow = length(f2),ncol = length(bin.loca))
  rownames(genotypes)<- unlist(strsplit(f2[],"_"))[c(seq(2,nrow(genotypes)*4,by=4))]
  colnames(genotypes)<-bin.loca[1:length(bin.loca)]+0.5*bin.size-1
  rownames(nmarkers)<- unlist(strsplit(f2[],"_"))[c(seq(2,nrow(genotypes)*4,by=4))]
  colnames(nmarkers)<-bin.loca[1:length(bin.loca)]+0.5*bin.size-1
  l <- 0
  for( i in 1:length(f2)){
    l <- l+1
    f2.id <- unlist(strsplit(f2[i],"_"))[2]
    pathout <- paste(" ","./result/",f2[i],".vcf",sep="")
    data1 <- read_trio(of_id = f2[i],F_id = F_id,M_id = M_id,vcf.file = vcf.file.f2,pathout = pathout,generate.input = F)
    setkey(data1,chr)
    data1.chr1 <- data1[k] # extract mrk from this chromsome 
    trac <- tracing_physical(input = data.frame(data1.chr1),bin = bin.size,chr.len = max(pos.chr1),cut = cut.number)
    genotypes[l,findInterval(trac$loca,bin.loca)] <- trac$f
    nmarkers[l,findInterval(trac$loca,bin.loca)] <- trac$fn
  }
}

# Heatmap informativity
heatmap.2(genotypes,Rowv = FALSE, Colv=FALSE)
heatmap.2(log(nmarkers),Rowv = FALSE, Colv=FALSE)

# Checking density between informative bins
dev.off()
F2_informative_density <- diff(trac$loca[trac$fn > 10])
hist(3*F2_informative_density/bin.size)

#Informativity aross chromosome
mrk_thresh <- 5
dev.off()
par(mfrow=c(3,1),mar=c(1,1,1,1),xpd=T)
plot(x = trac$loca,trac$fn,col="black",pch=19,frame.plot = F)
points(x=trac$loca[trac$fn > mrk_thresh],y = trac$fn[trac$fn > mrk_thresh],col="red",pch=19)

# Get linkage maps!


##########################################20171211-Yanjun ###############################
##This code generate a  data.frame with info on all the chromsomes
# define a few parameter
all.chr  <- read.table("~/Dropbox/Projects/Ongoing/HL_F2_SeqMrk/From.yanjun/data/chr_id.match.txt",stringsAsFactors = F,header = T,sep="\t")
chroms <- all.chr$INSDC
chroms.len <- all.chr$Size.Mb.*1e6
cut.number <- 5 # cutoff on each bin
bin.size <- 1e6 # size of each bin

# get.info is removed to function.R
out.info <- get.info(chroms = chroms,chroms.len = chroms.len,bin.size = bin.size)

## create the big data.frame, store each row for individual and each col for genotype in a bin
f2.id  <- gsub(pattern = "(F2_.*)_S.*",replacement = "\\1",x = f2)
genotype.hap1 <- data.frame(array(NA,dim = c(length(f2),max(out.info$index$end))))
rownames(genotype.hap1) <- f2.id
colnames(genotype.hap1) <- out.info$chr.loca
# here only hap1 data is stored because hap1+hap2 =1

# loop individual first so we read in each individual only once
for( i in 1:length(f2)){
  pathout <- paste(" ","./result/",f2[i],".vcf",sep="")
  data1 <- read_trio(of_id = f2[i],F_id = F_id,M_id = M_id,vcf.file = vcf.file.f2,pathout = pathout,generate.input = F)
  
  for (k in 1:length(chroms) ){
  chr  <- chroms[k]
  setkey(data1,chr)
  data.now <- data1[chr] # extract mrk from this chromsome 
  trac.now <- tracing_physical(input = data.frame(data.now),bin = bin.size,chr.len = chroms.len[k],cut = cut.number)
  # location at current loop
  loca.now <- out.info$loca[out.info$loca.chr==chr]
  idx.now <-c(out.info$index$start[k]:out.info$index$end[k])
  index.rep <- idx.now[findInterval(trac.now$loca/1e6,loca.now)]
  genotype.hap1[i,index.rep] <- trac.now$f
  }
  cat(f2[i],"done","\n")
}

require(gplots)
heatmap.2(data.matrix(genotype.hap1),Rowv = FALSE, Colv=FALSE,trace="none")
heatmap.2(log(nmarkers),Rowv = FALSE, Colv=FALSE)


############## compare the mrk increase from fixed between h/l and with grandparents
i=1
pathout <- paste(" ","./result/",f2[i],".vcf",sep="")
data1 <- read_trio(of_id = f2[f2.plot[i]],F_id = F_id,M_id = M_id,vcf.file = vcf.file.f2,pathout = pathout,generate.input = F)
pathout2 <- paste(" ","./result/",f2[i],"full.vcf",sep="")
py.script.gp <- "python3 ./bin/python/171215_read.vcf.grand.p.py " #/usr/local/Cellar/python3/3.6.3/bin/
##for the first f2 the ped are
#247 1655 1940 m f
#209 1833 2064 f m
# so ped file will be
M_id_h <- ".*1655.*"
F_id_l <- ".*1940.*"

F_id_h <- ".*1833.*"
M_id_l <- ".*2064.*"
# need to change the vcf input, this one is already filtered will not gain anythings
data2 <- read_grand.p(of_id = f2[i],F_id_h = F_id_h,M_id_h = M_id_h,F_id_l = F_id_l,M_id_l = M_id_l,vcf.file = vcf.file.f2,py = py.script.gp,pathout = pathout2)











#file.in <- fread("./data/171113.Fixed.sites.founder.all.chr.txt") # this is a file with all the founder fixed mrk

# Yanjun will do the coverage 

# Currently not used for analyses

## sliding window based estimation of genotypes
wins <- seq(from=10,to=200,by = 20)
dev.off()
for( j in 1:length(f2)){
  pathout <- paste(" ","./result/",f2[j],".vcf",sep="")
  data <- read_trio(of_id = f2[j],F_id = F_id,M_id = M_id,vcf.file = vcf.file.f2,pathout = pathout) # read in data
  # the data we read in only contain sites that are fixed between M-id F-id and have call in of-id
  setkey(data,chr) # using chromsome as key
  data.chr1 <- data["CM000093.4"] # extract mrk from this chromsome 
  
  # Visulise the result
  png(filename = paste("./result/",f2[j],".f2.sliding.png",sep = ""),width = 1024,height = 1024 )
  par(mfrow=c(11,1),mar=c(2,3,1,3),xpd=T)
  for( i in 1:length(wins)){
    win <- wins[i] 
    trc<- tracing(input = data.frame(data.chr1),step = win,stepsize = 1) # step size has to be 1 now, otherwise the physical location of this win is not implemented yet
    trc_step <- trc$counts
    pos.mean <- apply(trc$pos, 1, FUN=mean)/1e6 # using the mid point of start and end mrk in a win to represent the location of win
    hapf <- trc_step[,2]
    hapm <- trc_step[,3]
    
    
    # #based on physical postion
    # plot(pos.mean,hapf,pch=19,cex=0.8,main = win[i],frame.plot = F,ylim = c(0,1))
    # points(pos.mean,hapm,pch=19,cex=0.8,col="red")
    # abline(h=c(0.8,0.9),lty="dashed",col="red")
    # text(x = length(hapf[1:l])/2,y = 0.5,labels = wins[i],cex=2,col="blue")
    ## based on win index
    plot(hapf,pch=19,cex=0.8,main = win[i],frame.plot = F)
    points(hapm,pch=19,cex=0.8,col="red")
    abline(h=(0.8+0.9)/2,lty="dashed",col="red")
    text(x = length(hapf)/2,y = 0.5,labels = wins[i],cex=2,col="blue")
  }
  dev.off()
}

