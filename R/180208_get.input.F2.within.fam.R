# generate a ped input

founder <- read.table("~/Documents/impute/git/data/180222_founder.txt",stringsAsFactors = F,header = T,sep = "\t")
F1 <- read.table("~/Documents/impute/git/data/F1.ped.180208.txt",stringsAsFactors = F,header = T,sep = "\t")
F2 <- read.table("~/Documents/impute/git/data/F2.ped.180208.txt",stringsAsFactors = F,header = T,sep = "\t")

out <- data.frame(array(data = NA,dim = c(nrow(F2),6)))
hid <- founder[founder$Line=="HWS",1]
lid <- founder[founder$Line=="LWS",1]
sink("~/Documents/impute/git/data/Ped.f2.f2.f0.txt")
cat("id.f2","\t","id.f2.ma","\t","id.f2.fa","\t","fa.h","\t","ma.h","\t","ma.l","\t","fa.l","\n")

for( i in 1:nrow(F2)){
    id.f2 <- F2$Identity[i]
    id.f2.fa <- F2$Father[i]
    id.f2.ma <- F2$Mother[i]
    
    if( any(F1$Identity==id.f2.fa)){
      id.f2.fa.fa <- F1[F1$Identity==id.f2.fa,2]
      id.f2.fa.ma <- F1[F1$Identity==id.f2.fa,3]
      if(id.f2.fa.fa %in% hid){
        fa.h <- id.f2.fa.fa
        fa.l <- id.f2.fa.ma
      }else{
        fa.l <- id.f2.fa.fa
        fa.h <- id.f2.fa.ma
      }
      
      if(any(F1$Identity==id.f2.ma)){
        id.f2.ma.fa <- F1[F1$Identity==id.f2.ma,2]
        id.f2.ma.ma <- F1[F1$Identity==id.f2.ma,3]
        if(id.f2.ma.fa %in% hid){
          ma.h <- id.f2.ma.fa
          ma.l <- id.f2.ma.ma
        }else{
          ma.l <- id.f2.ma.fa
          ma.h <- id.f2.ma.ma
        }
      }
    cat(id.f2,"\t",id.f2.ma,"\t",id.f2.fa,"\t",fa.h,"\t",ma.h,"\t",ma.l,"\t",fa.l,"\n")
  
    }
    
}
sink()

# double check 
ped <-  read.table("~/Documents/impute/git/data/Ped.f2.f2.f0.txt",stringsAsFactors = F,header = T,sep = "\t")
sum(ped$fa.h %in% hid)
sum(ped$fa.l %in% lid)
sum(ped$ma.h %in% hid)
sum(ped$ma.l %in% lid)


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< on dip

read_grand.p <- function(of_id,F_id_h,M_id_h,F_id_l,M_id_l,vcf.file,py,pathout=NULL,generate.input=T){
  if(!require(data.table))
    require(data.table)
  cat("1.Please make sure the input F_id_h,M_id_h are from HWS and F_id_l,M_id_l are from LWS","\n")
  #1. check if the ped is right
  # not implemented 
  #2. read in data
  cat("2.Format data using python","\n")
  if(is.null(pathout)){
    if(generate.input){
      bash  <- paste(py,of_id," ",F_id_h," ",M_id_h," ",F_id_l," ",M_id_l,vcf.file ," ~/Documents/impute/data/testplate.F1/test.vcf",sep="")
      system(bash)
    }
    cat("3. Read in data using fread")
   # input <- fread("~/Documents/impute/data/testplate.F1/test.vcf")
  }else{
    if(generate.input){
      bash  <- paste(py,of_id," ",F_id_h," ",M_id_h," ",F_id_l," ",M_id_l,vcf.file ,pathout,sep="")
      system(bash)
    }
    cat("3. Read in data using fread")
    input <- fread(gsub(pattern = "\\s+(.*)",replacement="\\1",pathout))
  }
  colnames(input) <- c("chr","pos","ref","alt","f1","f2","m1","m2","of1","of2")
  return(input)
}




vcf.header <- function(vcf){
  vcf <- gsub(pattern = "\\s+(.*)",replacement="\\1",vcf)
  bash <- paste("gunzip -c ",vcf,"| awk '/#CHR/{print;exit}'",sep = "")
  b <- paste(system(bash,intern = T),collapse="")
  b.s <- strsplit(b,"\t")
  if(class(b.s)=="list"){
    return(unlist(b.s))
  }else{
    return(b.s)
  }
}

#F_id <- ".*1653.*"
#M_id <- ".*2109.*"
#vcf.file.f2 <- " /Users/yanjunzan/Documents/impute/git/data/180208.all.223+700.f2.P60.vcf.gz"
#py.script <- "python3 ~/Documents/impute/git/F2_re_seq/python/171215_read.vcf.grand.p.py " #/usr/local/Cellar/python3/3.6.3/bin/
vcf.file.f2 <- " /home/yanjun/projects/F2_seq/data/180208.all.223+700.f2.P60.vcf.gz"
py.script <- "python3 /home/yanjun/projects/F2_seq/data/171215_read.vcf.grand.p.py " #/usr/local/Cellar/python3/3.6.3/bin/

header.f2 <- vcf.header(vcf = vcf.file.f2) # check the vcf to getout the sample ID. there might be a broken pipe error, but ignore it
f2 <- header.f2[grep(pattern = "F2.*",x = header.f2)] ## will be all F2 in the end, for now a subset since all not available
id2 <- gsub(replacement = "\\1",pattern = "F2_(.*)_(S|m).*",x = f2)
id3<- as.numeric(gsub(replacement = "\\1",pattern = "(.*)_F2_(S|m).*",x = id2))

f2.ped <-  read.table("/home/yanjun/projects/F2_seq/data/Ped.f2.f2.f0.txt",stringsAsFactors = F,header = T,sep = "\t")
#f2.ped <-  read.table("~/Documents/impute/git/data/founder.180208.txt",stringsAsFactors = F,header = T,sep = "\t")

f2.sub <- f2[id3 %in% f2.ped$id.f2]
id2.sub <- gsub(replacement = "\\1",pattern = "F2_(.*)_(S|m).*",x = f2.sub)
id3.sub<- as.numeric(gsub(replacement = "\\1",pattern = "(.*)_F2_(S|m).*",x = id2.sub))
f2.ped.sub <-f2.ped[f2.ped$id.f2 %in% id3.sub,]
f2.ped.sub.2<-  f2.ped.sub[match(table= f2.ped.sub$id.f2,x= id3.sub),]
setwd("~/projects/F2_seq")
library(doSNOW)
library(foreach)
cl<-makeCluster(20) #change 7 to the number of CPU cores you want 
registerDoSNOW(cl)

foreach(i=1:nrow(f2.ped.sub.2)) %dopar% {
  pathout <- paste(" ","/home/yanjun/projects/F2_seq/data/with.fam.f2/",f2.ped.sub.2$id.f2[i],".vcf",sep="")
  #pathout <- paste(" ","~/Documents/impute/data/180215_",f2.ped.sub.2$id.f2[i],".vcf",sep="")
  of_id.1 <- paste("(^F2_",f2.ped.sub.2$id.f2[i],".*)|(^",f2.ped.sub.2$id.f2[i],"_F2.*)",sep="")
  of_id.2 <- header.f2[grep(pattern = of_id.1,x = header.f2)]
  F_id_h <- paste(".*",f2.ped.sub.2$fa.h[i],".*",sep="")
  M_id_h <- paste(".*",f2.ped.sub.2$ma.h[i],".*",sep="")
  F_id_l <- paste(".*",f2.ped.sub.2$ma.l[i],".*",sep="")
  M_id_l <- paste(".*",f2.ped.sub.2$fa.l[i],".*",sep="")
  read_grand.p(of_id = of_id.2,F_id_h =F_id_h,M_id_h = M_id_h,F_id_l=F_id_l,M_id_l=M_id_l,vcf.file = vcf.file.f2,pathout = pathout,py = py.script,generate.input = T)
}
stopCluster(cl)




source("~/projects/F2_seq/F2_re_seq/R/Functions.AIL.impute.R")
all.chr  <- read.table("/home/yanjun/projects/F2_seq/F2_re_seq/data/chr_id.match.txt",stringsAsFactors = F,header = T,sep="\t")
chroms <- all.chr$INSDC
chroms.len <- all.chr$Size.Mb.*1e6
cut.number <- 50 # cutoff on each bin
bin.size <- 1e6 # size of each bin

# get.info is removed to function.R
out.info <- get.info(chroms = chroms,chroms.len = chroms.len,bin.size = bin.size)

## create the big data.frame, store each row for individual and each col for genotype in a bin
f2.id  <- f2.ped.sub.2$id.f2 #gsub(pattern = "(.*)_(S|m).*",replacement = "\\1",x = f2)
genotype.hap1 <- data.frame(array(NA,dim = c(length(f2.id),max(out.info$index$end))))
rownames(genotype.hap1) <- f2.id
colnames(genotype.hap1) <- out.info$chr.loca
# here only hap1 data is stored because hap1+hap2 =1
# 208 is empty

# loop individual first so we read in each individual only once
for( i in 1:length(f2.ped.sub.2$id.f2)){
  pathout <- paste(" ","/home/yanjun/projects/F2_seq/data/with.fam.f2/",f2.ped.sub.2$id.f2[i],".vcf",sep="")
  of_id <- paste(".*",f2.ped.sub.2$id.f2[i],".*",sep="")
  F_id_h <- paste(".*",f2.ped.sub.2$fa.h[i],".*",sep="")
  M_id_h <- paste(".*",f2.ped.sub.2$ma.h[i],".*",sep="")
  F_id_l <- paste(".*",f2.ped.sub.2$ma.l[i],".*",sep="")
  M_id_l <- paste(".*",f2.ped.sub.2$fa.l[i],".*",sep="")
  
  data1 <-  try(silent = T,read_grand.p(of_id = of_id,F_id_h =F_id_h,M_id_h = M_id_h,F_id_l=F_id_l,M_id_l=M_id_l,vcf.file = vcf.file.f2,pathout = pathout,py = py.script,generate.input = F))

  if(!inherits(data1,"try-error")){
    
   
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
  }
  cat("\n",i,"-",f2.ped.sub.2$id.f2[i],"done","\n")
}

dim(genotype.hap1)
write.table(genotype.hap1,file = "~/projects/F2_seq/data/F2_827.within.fam.bin50_bugcorrected.180223.txt",quote = F,sep = "\t")

genotype.hap1 <- read.table("~/Documents/impute/data/F2_827.within.fam.bin50_bugcorrected.180223.txt",header=T,stringsAsFactors = F)
geno <- genotype.hap1
upper.cut <- 0.8
lower.cut <- 0.2
arbitary.cut <- function(geno,upper.cut,lower.cut){
  for( i in 1:nrow(geno)){
    id.h <- which(geno[i,] >= upper.cut)
    geno[i,id.h] <- "A"
    id.l <- which(geno[i,] <= lower.cut)
    geno[i,id.l] <- "B"
    id.hl <- which(geno[i,] > lower.cut & geno[i,] < upper.cut)
    geno[i,id.hl] <- "H"
    cat(i,"\n")
  }
  return(geno)
}

arbitary.cut.missing <- function(geno,upper.cut,lower.cut,missing.upper,missing.lower){
  for( i in 1:nrow(geno)){
    #set homo
    id.h <- which(geno[i,] >= upper.cut)
    geno[i,id.h] <- 1
    id.l <- which(geno[i,] <= lower.cut)
    geno[i,id.l] <- -1
    #set het
    id.hl <- which(geno[i,] >= missing.lower & geno[i,] <= missing.upper)
    geno[i,id.hl] <- 0
    #set missing
    id.missing1 <- which(geno[i,] > missing.upper & geno[i,] < upper.cut)
    geno[i,id.missing1] <- NA
    id.missing2 <- which(geno[i,] > lower.cut & geno[i,] < missing.lower)
    geno[i,id.missing2] <- NA
    
    cat(i,"\n")
  }
  return(geno)
}

test <- arbitary.cut(geno = genotype.hap1,upper.cut = 0.8,lower.cut = 0.2) # this is cut above 0.8 to homozygous high below 0.2 to homozygous low and in the middle to me heterzogous

# add first and second line
line1 <- c("body_weight",colnames(test))
chr <- gsub(pattern = "(.*\\.{1}\\d+)\\.\\d{1,3}\\.\\d",replacement="\\1",x=colnames(test))
all.chr  <- read.table("~/Documents/AIL/From.yanjun/data/chr_id.match.txt",stringsAsFactors = F,header = T,sep="\t")
chr <- all.chr$Name[ match(chr,all.chr$INSDC)]
chr[grep(pattern = "W",chr)] <- "34"
chr[grep(pattern = "Z",chr)] <- "35"
chr[grep(pattern = "LGE64",chr)] <- "36"
chr[grep(pattern = "MT",chr)] <- "37"
line2 <- c("",chr)


test2 <- data.frame(array(NA,dim = c(nrow(test)+2,length(line1))))
test2[1,] <-line1
test2[2,] <- line2
test2[3:nrow(test2),2:length(line1)] <- test
test2[3:nrow(test2),1] <- rownames(test)
#rbind.data.frame(line1,line2,test1,)
write.table(test2,file = "~/Documents/impute/data/F2_827.within.fam.mat.bin50.cut.bugcorrected.180223.csv",quote = F,sep = ";",row.names = F,col.names = F)

# below is 1-0.85 homozygous hig 0.85-0.7 NA,0.7-0.3 heterzogous, 0.3-0.15 NA,below 0.15 hom low
test2 <-arbitary.cut.missing(geno = genotype.hap1,upper.cut = 0.85,lower.cut = 0.15,missing.upper = 0.7,missing.lower = 0.3)

write.table(test,file = "~/projects/F2_seq/data/F2_827.within.fam.mat.bin50.cut.txt",quote = F,sep = "\t")
require(gplots)
heatmap.2(data.matrix(test),Rowv = FALSE, Colv=FALSE,trace="none")
############# check the raw data
id <- 380
folder <- "~/Documents/impute/data/"
chr <- "CM000093.4"
start <- 45
end <- 51
step <- 2
extract.raw <- function(id,folder,chr,start,end,step){
  if(!require(data.table))
    require(data.table)
  file <- list.files(path=folder,pattern = paste(id,".vcf",sep=""))
  data.in <- fread(paste(folder,file,sep=""),header = T)
  colnames(data.in) <- c("chr","pos","ref","alt","f1","f2","m1","m2","of1","of2")
  
  setkey(data.in,chr)
  data.chr <- data.in[chr]
  out <- data.chr[which(data.chr$pos > (start-step)*1e6 & data.chr$pos < (end+step)*1e6),]
  return(data.frame(out))
}

test.out <- extract.raw(id = 380,folder = "~/Documents/impute/data/",chr = "CM000093.4",start = 47.5,end = 49.5,step = 2)
View(test.out)








test <- read.table("data/F2_827.within.fam.mat.bin15.txt",sep = "\t",header = T,stringsAsFactors = F)
require(gplots)
heatmap.2(data.matrix(test),Rowv = FALSE, Colv=FALSE,trace="none")
count.na <- function(mat,y){
  return(apply(mat, y, FUN = function(x) sum(!is.na(x))))
}
rna <- count.na(mat = test,y = 1)
cna <- count.na(mat = test,y=2)



# remove F1 first
f1 <- read.table("~/Documents/impute/From.yanjun/data/Pedigree.F1.txt",stringsAsFactors = F,header = T,sep = "\t")
r.id <- rownames(test)
r.id.num  <- as.numeric(sub(pattern = "(F2_|_F2)",replacement = "",x =r.id ))
f2.all <- r.id[!(r.id.num %in% f1$F1)]
test.f2 <- test[f2.all,]
write.table(test.f2,file = "~/Documents/impute/From.yanjun/data/F2_n707.test.geno.mat.txt",quote = F,sep = "\t")

# read in phynotype
phe <- read.table(file = "~/Documents/impute/From.yanjun/data/061106_hl_f2_phenotypes.txt",stringsAsFactors = F,header = T,sep = "\t",fill=T)
phe <- phe[!is.na(phe[,2]),]
id.match <- intersect(r.id.num,phe$id)

f2.phe <- r.id[(r.id.num %in% id.match)]
f2.phe.id <- as.numeric(sub(pattern = "(F2_|_F2)",replacement = "",x =f2.phe ))
pheno <- phe[match(f2.phe.id,phe$id),]
geno <- test.f2[f2.phe,cna>50]

# get sex
ped.f2 <- read.table("~/Documents/impute/From.yanjun/data/f2.ped.txt",stringsAsFactors =F,header = T,sep = "\t")
sex <- ped.f2$sex[match(f2.phe.id,ped.f2$id)]
View(cbind(ped.f2[match(f2.phe.id,ped.f2$id),]$id,f2.phe.id,rownames(geno)))
get.chr.loca <- function(geno){
  chr.t<-sub(pattern = "\\.",replacement="-",x=colnames(geno))
  chr <- gsub(pattern = "(.*-\\d+)(.*)",replacement="\\1",x=chr.t)
  chr <- sub(pattern = "-",replacement=".",chr)
  loca <- as.numeric(gsub(pattern = "(.*-\\d+)\\.(.*)",replacement="\\2",x=chr.t))
  return(list("chr"=chr,"loca"=loca))
}

p <- numeric(ncol(geno))
names(p) <- colnames(geno)
for( i in 1:ncol(geno)){
  lm <- summary(lm(formula = pheno$BW8~sex+geno[,i]) )
  p[i] <- lm$coefficients[3,4]
  cat(i,"\n")
}

# plot function
mat.full <- test.f2
mat.scan <- geno # test.f2[f2.phe.id,cna>50]
index <- out.info$num.bin
all.chr  <- read.table("~/Documents/impute/From.yanjun/data/chr_id.match.txt",stringsAsFactors = F,header = T,sep="\t")

plot.scan.f2 <- function(mat.scan,p,out.info){
  xlim <- c(0,sum(out.info$num.bin))
  ylim <- c(0,max(-log10(p)))
  p.full <- rep(NA,sum(out.info$num.bin))
  names(p.full) <- sub(pattern = "-",replacement = "\\.",out.info$chr.loca)
  p.full[colnames(mat.scan)] <- p
  
  plot(x = 1:sum(out.info$num.bin),y = -log10(p.full),xlim=xlim,ylim=ylim,frame.plot=F,ylab="-Log10 (p-value)",xlab="Chromosome",xaxt="n")
  loca.tick <- c(1,out.info$index[,2])
  axis(side = 1,at = loca.tick,labels = T)
}