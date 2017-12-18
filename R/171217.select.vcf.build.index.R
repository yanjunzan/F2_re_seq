path1 <- "/tmp/nas/yanjun/AIL/1711recall.vcf/vcffile"
vcf <- list.files(path = path1,pattern = "all.vcf.recode.vcf$")
vcf.idx <- list.files(path = path1,pattern = "all.vcf.recode.vcf.idx$")
id.vcf <- gsub(pattern = "(.*)_S.*",replacement = "\\1",x = vcf)
id.vcf.idx <- gsub(pattern = "(.*)_S.*",replacement = "\\1",x = vcf.idx)
out <- vcf[!(id.vcf %in% id.vcf.idx)]
#write.table(out,file="/tmp/nas/yanjun/AIL/1711recall.vcf/vcffile/buildindex.txt",sep="\t",col.names=FALSE,row.names = F,quote = F)

library(doSNOW)
library(foreach)
cl<-makeCluster(10) #change 7 to the number of CPU cores you want 
registerDoSNOW(cl)

foreach(i=1:length(file.todo)) %dopar% {
  line <- out[i]
  bash <- paste("bash /tmp/nas/yanjun/AIL/1711recall.vcf/171218_ge.index.sh ",line,sep="") # this pass the bam file to the bash scripts, and parallel it
  system(bash)
}
stopCluster(cl)
