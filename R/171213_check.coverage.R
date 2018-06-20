require(data.table)

cov <- file("~/Documents/impute/F2.all.coverage1212.txt")
open(cov)
num.line <-  unlist(strsplit(system("wc -l ~/Documents/impute/F2.all.coverage1212.txt",intern = T),split = " "))[5]
num.line <- as.numeric(num.line)/2
id <- numeric(num.line)
coverage1 <- numeric(num.line)
coverage2 <- numeric(num.line)
for ( i in 1:num.line){
  line <- readLines(cov, n = 2, warn = FALSE)
  id[i] <- gsub(replacement = "\\1",pattern = ".*bamfile/(.*)_.*",x = line[1])
  coverage1[i] <- as.numeric(gsub(replacement = "\\1",pattern = ".*=\\s{1}(.*\\d+).*c.*",x = line[2]))
  coverage2[i] <- as.numeric(gsub(replacement = "\\1",pattern = ".*c.*=\\s{1}(.*\\d+)",x = line[2]))
}
close(cov)
#id2 <- 
id2 <- gsub(replacement = "\\1",pattern = "F2_(.*)_(S|m).*",x = id)
id3<- as.numeric(gsub(replacement = "\\1",pattern = "(.*)_F2_(S|m).*",x = id2))


par(mfrow=c(2,2))
plot(c(1:num.line),y = coverage1,col="red",ylab = " coverage",xlab = "indivdual",pch=19,cex=0.5,frame.plot = F)
abline(h=mean(coverage1,na.rm=T),col="red")
abline(h=0.4,lty="dashed")
abline(h=0.6,lty="dashed")
abline(h=0.8,lty="dashed")

plot(c(1:num.line),y = coverage2,col="red",ylab = " site coverage",xlab = "indivdual",pch=19,cex=0.5,frame.plot = F)
abline(h=mean(coverage2,na.rm=T),col="red")
abline(h=0.4,lty="dashed")
abline(h=0.6,lty="dashed")
abline(h=0.8,lty="dashed")
#legend("topright",fill = c("red","green"),legend = c("mean.coverage","site covered"))
plot(coverage2,coverage1,frame.plot = F,xlab = "site coverage",ylab = "coverage")
close(cov)