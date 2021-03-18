data<- read.csv('NRAP_16S_SAE3_Raw_Data_195.csv', header = TRUE, stringsAsFactors = FALSE)
data=data[which(data$sum>50),]
head(data)
n.gen=length(data$Genus)
gens <-data$Genus
gens

data=data.frame(data)
dat0=log(data[,8:18]+.01)
dat0=data.matrix(dat0)
#data0=dat0/mean(dat0)
dat5=log(data[,19:29]+.01)
dat5=data.matrix(dat5)
#data5=dat5/mean(dat5)
dat23=log(data[,30:40]+.01)
dat23=data.matrix(dat23)
#data23=dat23/mean(dat23)
dat37=log(data[,41:51]+.01)
dat37=data.matrix(dat37)
#data37=dat37/mean(dat37)
dat68=log(data[,52:56]+.01)
dat68=data.matrix(dat68)
#data68=dat68/mean(dat68)
data1=cbind(dat0,dat5,dat23,dat37,dat68)
data1=data.frame(data1)
allmeans =cbind(apply(dat0,1,mean), apply(dat5,1,mean), apply(dat23,1,mean), apply(dat37,1,mean),apply(dat68,1,mean))
normdat=data1/allmeans
zero=rep(0,11)
five=rep(5,11)
tweentythree=rep(23,11)
thirtyseven=rep(37,11)
sixtyeigth=rep(68,5)
timegroup=c(zero,five,tweentythree,thirtyseven, sixtyeigth) 
##Transpose
# The first column is actually the row names
rownames(data1) <- gens
#colnames(data1)<-timegroup
head(data1)
#my_data<-data1  # this is for SAE only
#head(my_data)
library(data.table)
mydata_transpose <- transpose(data1)
head(mydata_transpose)
#However we have to get back the row and col names
#rownames(mydata_transpose) <- colnames(data1)
colnames(mydata_transpose) <- rownames(data1)
mydata=mydata_transpose
head(mydata)

pval=rep(1,n.gen)
for (i in 1:n.gen){
y= mydata[,i]

dat=cbind(timegroup,y)
dat=data.frame(dat)
# Compute the analysis of variance
anv <- aov(y~factor(timegroup), data=dat)
# Summary of the analysis
summary(anv)
pval[i]=summary(anv)[[1]][["Pr(>F)"]]



}
print(pval)

w1 = which(pval<0.05/n.gen)    
gens[w1]
w2=which(pval>0.05/n.gen)
gens[w2]
#Significant Genus 123
mydata$"Acinetobacter"=NULL
mydata$"Polynucleobacter"=NULL
length(mydata)
install.packages("lawstat")
library(lawstat)
pvalue=rep(1,length(w1))
#####PLots
pdf(file=paste("plotsofgenus.pdf"))

     for (i in 1:length(w1)){
          thegenus = gens[w1][i]
          y= mydata[,i]
          data2=cbind(timegroup,y)
          data2=data.frame(data2)
##Anova 
anv <- aov(y~factor(timegroup), data=data2)

# Summary of the analysis
               summary(anv)
##Test of homogenity of variance
#Levene test
lvtest=levene.test(y, factor(timegroup), location="median", correction.method="zero.correction")
pvalue[i]=lvtest[[2]][1]
###pvalue should be less than 0.05.
               mse=summary(anv)[[1]][["Mean Sq"]]
               mse= mse[2]
               b= length(timegroup)-5 # Df=N-5
		   month=c(0,5,23,37,68)
Mean_log_reads=c(mean(data2[,2][timegroup==0]),mean(data2[,2][timegroup==5]),mean(data2[,2][timegroup==23]),mean(data2[,2][timegroup==37]),mean(data2[,2][timegroup==68]))
plot(month,Mean_log_reads, main=paste("Genus",thegenus), type="o" ,xlim = c(0, 70), ylim = c(-10, 10))
abline(h = 1, col = "blue")
lines(c(0,0), c(Mean_log_reads[1]- qt(0.95,b)*sqrt(mse/11) ,Mean_log_reads[1]+qt(0.95,b)*sqrt(mse/11)))
lines(c(5,5), c(Mean_log_reads[2]- qt(0.95,b)*sqrt(mse/11) ,Mean_log_reads[2]+qt(0.95,b)*sqrt(mse/11)))
lines(c(23,23), c(Mean_log_reads[3]- qt(0.95,b)*sqrt(mse/11) ,Mean_log_reads[3]+qt(0.95,b)*sqrt(mse/11)))
lines(c(37,37), c(Mean_log_reads[4]- qt(0.95,b)*sqrt(mse/11) ,Mean_log_reads[4]+qt(0.95,b)*sqrt(mse/11)))
lines(c(68,68), c(Mean_log_reads[5]- qt(0.95,b)*sqrt(mse/5) ,Mean_log_reads[5]+qt(0.95,b)*sqrt(mse/5)))
}

dev.off()


print(pvalue)
gens[which(pvalue>0.05)]
sink("confidenceInterval.txt")
####CI
     for (i in 1:length(w1)){
          thegenus = gens[w1][i]
          y= mydata[,i]
          data2=cbind(timegroup,y)
          data2=data.frame(data2)
##Anova 
anv <- aov(y~factor(timegroup), data=data2)
##CI
CI=confint(anv, method="Wald")

print(CI)
}
sink()


