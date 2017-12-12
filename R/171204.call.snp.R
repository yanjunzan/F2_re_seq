#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< SNP calling using GATK UG<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# This reads in a bamfile list and parallel the bash script to call snp for each bam file
library(doSNOW)
library(foreach)
cl<-makeCluster(10) #change 7 to the number of CPU cores you want 
registerDoSNOW(cl)
file  <- read.table(file = "/tmp/nas/yanjun/AIL/founder.hap/ref.chr1.171030/171030.all.f2.bam.list.txt",sep="\t",header = F,stringsAsFactors = F)
name.all <- gsub(pattern = ".*bamfile/(.*)_(S.*|mer.*).*",replacement = "\\1",x = file[,1]) 
file.done <- list.files("/tmp/nas/yanjun/AIL/1711recall.vcf/vcffile",pattern = "*_all.vcf.log")
name.done <- gsub(pattern = "(.*)_S\\d+_L.*",replacement = "\\1",x = file.done) 
file.todo <- file[!(name.all %in% name.done),1]

foreach(i=1:length(file.todo)) %dopar% {
  line <- file.todo[i]
  bash <- paste("bash /tmp/nas/yanjun/AIL/1711recall.vcf/171103.recall.snp.sh ",line,sep="") # this pass the bam file to the bash scripts, and parallel it
  system(bash)
}
stopCluster(cl)

