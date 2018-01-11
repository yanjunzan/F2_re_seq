setwd("~/projects/F2_seq")
library(doSNOW)
library(foreach)
cl<-makeCluster(30) #change 7 to the number of CPU cores you want 
registerDoSNOW(cl)
read_trio <- function(of_id,F_id,M_id,vcf.file,pathout=NULL){
  #if(!require(data.table))
  #  require(data.table)
  bash  <- paste(py.script,of_id," ",F_id," ",M_id,vcf.file ," ",pathout,sep="")
  system(bash)
  cat(paste(of_id,"done","\n",sep="\t"))
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
F_id <- ".*1653.*"
M_id <- ".*2109.*"
vcf.file.f2 <- " ~/projects/F2_seq/data/171215_all.780.F0.output.recode.vcf.gz"
py.script <- "python3 ~/projects/F2_seq/F2_re_seq/python/Read.vcf.py " #/usr/local/Cellar/python3/3.6.3/bin/
header.f2 <- vcf.header(vcf = vcf.file.f2) # check the vcf to getout the sample ID. there might be a broken pipe error, but ignore it
f2 <- header.f2[grep(pattern = "F2.*",x = header.f2)] ## will be all F2 in the end, for now a subset since all not available


foreach(i=1:length(f2),.export = "read_trio") %dopar% {
  pathout <- paste(" ","~/projects/F2_seq/results/",f2[i],".vcf",sep="")
  read_trio(of_id = f2[i],F_id = F_id,M_id = M_id,vcf.file = vcf.file.f2,pathout = pathout)
}
stopCluster(cl)




source("~/projects/F2_seq/F2_re_seq/R/Functions.AIL.impute.R")
all.chr  <- read.table("/home/yanjun/projects/F2_seq/F2_re_seq/data/chr_id.match.txt",stringsAsFactors = F,header = T,sep="\t")
chroms <- all.chr$INSDC
chroms.len <- all.chr$Size.Mb.*1e6
cut.number <- 5 # cutoff on each bin
bin.size <- 1e6 # size of each bin

# get.info is removed to function.R
out.info <- get.info(chroms = chroms,chroms.len = chroms.len,bin.size = bin.size)

## create the big data.frame, store each row for individual and each col for genotype in a bin
f2.id  <- gsub(pattern = "(.*)_(S|m).*",replacement = "\\1",x = f2)
genotype.hap1 <- data.frame(array(NA,dim = c(length(f2),max(out.info$index$end))))
rownames(genotype.hap1) <- f2.id
colnames(genotype.hap1) <- out.info$chr.loca
# here only hap1 data is stored because hap1+hap2 =1
# 208 is empty

# loop individual first so we read in each individual only once
for( i in 209:length(f2)){
  pathout <- paste(" ","~/projects/F2_seq/results/",f2[i],".vcf",sep="")
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

dim(genotype.hap1)
write.table(genotype.hap1,file = "~/projects/F2_seq/results/F2_n780.test.count.mat.txt",quote = F,sep = "\t")
geno <- genotype.hap1
upper.cut <- 0.8
lower.cut <- 0.2
arbitary.cut <- function(geno,upper.cut,lower.cut){
  for( i in 1:nrow(geno)){
    id.h <- which(geno[i,] >= upper.cut)
    geno[i,id.h] <- 1
    id.l <- which(geno[i,] <= lower.cut)
    geno[i,id.l] <- -1
    id.hl <- which(geno[i,] > lower.cut & geno[i,] < upper.cut)
    geno[i,id.hl] <- 0
    cat(i,"\n")
  }
  return(geno)
}

test <- arbitary.cut(geno = genotype.hap1,upper.cut = 0.8,lower.cut = 0.2)

write.table(test,file = "~/projects/F2_seq/results/F2_n780.test.geno.mat.txt",quote = F,sep = "\t")

test <- read.table("data/F2_n780.test.geno.mat.txt",sep = "\t",header = T,stringsAsFactors = F)
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
