**Biostatistics and data visualization using R**

**Contributors**

Qingxiang Guo

**About**

The following contents contains source code for the statistical analyses using R. Some codes of the statistical analyses and plotting using R are listed in specific section like evolutionary and comparative genomics. So only general analyses are listed here.

**Notes**

Written by Qingxiang Guo, qingxiang.guo@outlook.com, distributed without any guarantees or restrictions.

**Codes**

**1. T-test**

\# A t-test is necessary for small samples because their distributions are not normal (n < 30).

\# Test whether the data fit with normal distribution prior to T-test

**1.1 One Sample t-test**

\# To determine whether an unknown population mean is different from a specific value

x=c(187,183,175,178,180,175,174,178)

t.test(x,mu=169,alternative="greater")   

**1.2 Paired Samples t Test**

\# The Paired Samples t Test compares the means of two measurements taken from the same individual, object, or related units

x=c(9.7 ,6.2 ,7.0 ,5.3 ,8.1 ,9.9 ,4.7 ,5.8 ,7.8 ,8.6 ,6.1,9.9)

y=c(6.7 ,5.4 ,5.7 ,5.0 ,7.5 ,8.3 ,4.6 ,4.2 ,7.5 ,7.0 ,5.3,10.3)

t.test(x,y,paired=TRUE)

**1.3 Two-sample independent t-test**

\# To test whether the unknown population means of two groups are equal or not

x=c(78,82,76,74,77,78,76,77,81,83,79,85)

y=c(75,76,72,75,73,71,70,73,80,75,72)

t.test(x,y,alternative="greater")

**1.4 Use T-test to analyze the differential expression between group C and group M**

a=read.table("all\_fpkm",header=T,sep="\t")

Pvalue<-c(rep(0,nrow(a)))

log2\_FC<-c(rep(0,nrow(a)))

for(i in 1:nrow(a)){

y=t.test(a[i,2:4],a[i,5:7])

Pvalue[i]<-y$p.value # Output the P\_value

log2\_FC[i]<-log2(mean(as.numeric(a[i,2:4]))/mean(as.numeric(a[i,5:7]))) # output the fold\_change

}

out<-cbind(a,log2\_FC,Pvalue)

write.table(out,file="ttest.out.xls",quote=FALSE,sep="\t",row.names=FALSE)

keep=out$Pvalue<0.05 & !is.na(out$Pvalue) & abs(out$log2\_FC)>1 

out1=out[keep,]

write.table(out1,file="ttest.out.diff.xls",quote=FALSE,sep="\t",row.names=FALSE)

**2. ANOVA**

\# Do homogeneity test for variance before ANOVA, pass if P > 0.05

bartlett.test(exp~type)

\# General ANOVA

x=c(2,4,3,2,4,7,7,2,5,4)

y=c(5,6,8,5,10,7,12,6,6)

z=c(7,11,6,6,7,9,5,10,6,3,10)

exp=c(x,y,z)

type=factor(c(rep(“a”,length(x)),rep(“b”,length(y)),rep(“c”,length(z))))

ba.an=aov(exp~type)

TukeyHSD(ba.an) 

summary(ba.an)

\# Use ANOVA to analyze the differential expression between group C and group M

a=read.table("all\_fpkm",header=T,sep="\t")

Pvalue<-c(rep(0,nrow(a)))

log2\_FC<-c(rep(0,nrow(a)))

type<-factor(c(rep("c",3),rep("m",3)))

for(i in 1:nrow(a)){

y=aov(as.numeric(a[i,2:7])~type)

Pvalue[i]<-summary(y)[[1]][,5][1]

log2\_FC[i]<-log2(sum(a[i,2:4])/sum(a[i,5:7]))

}

out<-cbind(a,log2\_FC,Pvalue)

write.table(out,file="aov.out.xls",quote=FALSE,sep="\t",row.names=FALSE)

keep=out$Pvalue<0.05 & !is.na(out$Pvalue) &abs(out$log2\_FC)>1

out1=out[keep,]

write.table(out1,file="aov.out.diff.xls",quote=FALSE,sep="\t",row.names=FALSE)

**3. Chi-square test**

\# For independence compares two variables in a contingency table to see if they are related.

x=c(60,3,32,11)

dim(x)=c(2,2)

chisq.test(x, correct=F)

**4. Mann–Whitney U test**

\# The Mann-Whitney U test is used to compare differences between two independent groups when the dependent variable is either ordinal or continuous, but not normally distributed

x=c(24,26,29,34,43,58,63,72,87,101)

y=c(82,87,97,121,164,208,213)

wilcox.test(x, y, exact=FALSE, correct=FALSE)

\# Use Mann–Whitney U test to analyze the differential expression between group C and group M

a=read.table("all\_fpkm",header=T,sep="\t")

Pvalue<-c(rep(0,nrow(a)))

log2\_FC<-c(rep(0,nrow(a)))

for(i in 1:nrow(a)){

y=wilcox.test (as.numeric(a[i,2:4]), as.numeric(a[i,5:7]),exact=FALSE,correct=FALSE)

Pvalue[i]<-y$p.value

log2\_FC[i]<-log2(sum(a[i,2:4])/sum(a[i,5:7]))

}

out<-cbind(a,log2\_FC,Pvalue)

write.table(out,file="wilcox.out.xls",quote=FALSE,sep="\t",row.names=FALSE)

keep=out$Pvalue<0.05 & !is.na(out$Pvalue) & abs(out$log2\_FC)>1

out1=out[keep,]

write.table(out1,file="wilcox.out.diff.xls",quote=FALSE,sep="\t",row.names=FALSE)

**5. Pearson and Spearman correlation coefficients**

a=read.table("all\_fpkm",header=T,sep="\t",row.names=1)

cor(a$C2\_FPKM,a$C3\_FPKM,method="pearson")

cor(a$C2\_FPKM,a$C3\_FPKM,method="spearman")

cor(a$C2\_FPKM,a$M2\_FPKM,method="spearman")

cor(a$C2\_FPKM,a$M2\_FPKM,method="pearson")

**6. Cosine analysis in R**

library(ggthemes)

library(ggplot2)

library(cosinor)

data<-read.csv("ccc.csv",header=T) 

data <- data.frame(A=c(data[1]),B=c(data[2]), C=c(data[3]))

\# Make dot plot

qplot(time, Y,data=data,xlab ="",ylab= "", ylim=c(0,2),width=0.5,colour="red")+geom\_line(size=1)+theme(panel.grid.major = element\_blank(),panel.grid.minor = element\_blank())+geom\_errorbar(aes(ymin = Y- sd, ymax = Y +sd),width=0.2)+ theme(legend.position='none')

fit <- cosinor.lm(Y ~ time(time), period=7, data = data)

summary(fit)

ggplot.cosinor.lm(fit)+ theme\_set(theme\_bw())+theme(panel.grid.major = element\_blank(),panel.grid.minor = element\_blank())

**7. PCA analysis for RPKM and FPKM gene expression data**

**7.1 RPKM**

library(gmodels)

`	`inname = "rpkm.txt"

`	`# out PCA figure name \*\* 

`	`outname = "all\_DEGs.PCA.pdf"

`	`# define the color for points  \*\*

`	`mycolors <- c(rep("red",2),rep("green",2),rep("yellow",2),rep("blue",2))

`	`# read the data

`	`expr1 <- read.table(inname, header=T, row.names=1)

`	`expr = t(scale(t(expr1)))

`	`# transpose the data

`	`data <- t(expr)

`	`# do PCA 

`	`data.pca <- fast.prcomp(data)

`	`#data.pca <- fast.prcomp(data,retx=T,scale=F,center=T)

`	`# fetch the proportion of PC1 and PC2

`	`a <- summary(data.pca)

`	`tmp <- a[4]$importance

`	`pro1 <- as.numeric(sprintf("%.3f",tmp[2,1]))\*100

`	`pro2 <- as.numeric(sprintf("%.3f",tmp[2,2]))\*100

`	`# fetch the x axis min and max value (PC1)

`	`xmax <- max(data.pca$x[,1])

`	`xmin <- min(data.pca$x[,1])

`	`# fetch the y axis min and max value (PC2)

`	`ymax <- max(data.pca$x[,2])

`	`ymin <- min(data.pca$x[,2])

`	`# fetch sample names

`	`samples =rownames(data.pca$x)

`	`# draw PCA plot figure

`	`pdf(outname)

`	`plot(

`	`data.pca$x[,1],

`	`data.pca$x[,2],

`	`xlab=paste("PC1","(",pro1,"%)",sep=""),

`	`ylab=paste("PC2","(",pro2,"%)",sep=""),

`	`main="PCA",

`	`xlim=c(1.1\*xmin,1.1\*xmax),

`	`ylim=c(1.1\*ymin,1.1\*ymax),

`	`pch=16,col=mycolors)

`	`abline(h=0,col="gray")

`	`abline(v=0,col="gray")

`	`text(data.pca$x[,1],data.pca$x[,2],labels=samples)

`	`dev.off()

**7.2 FPKM**

library(gmodels)

\# input file name \*\*

inname = "degs.fpkm"

\# out PCA figure name \*\* 

outname = "TAIR\_DEGs.PCA.pdf"

\# define the color for points  \*\*

mycolors <- c(rep("red",3), rep("blue",3))

\# read the expr data

expr <- read.table(inname, header=T, row.names=1)

\# transpose the data

data <- t(expr)

\# do PCA 

data.pca <- fast.prcomp(data)

#data.pca <- fast.prcomp(data,retx=T,scale=F,center=T)

\# fetch the proportion of PC1 and PC2

a <- summary(data.pca)

tmp <- a[4]$importance

pro1 <- as.numeric(sprintf("%.3f",tmp[2,1]))\*100

pro2 <- as.numeric(sprintf("%.3f",tmp[2,2]))\*100

\# fetch the x axis min and max value (PC1)

xmax <- max(data.pca$x[,1])

xmin <- min(data.pca$x[,1])

\# fetch the y axis min and max value (PC2)

ymax <- max(data.pca$x[,2])

ymin <- min(data.pca$x[,2])

\# fetch sample names

samples =rownames(data.pca$x)

\# draw PCA plot figure

pdf(outname)

plot(

data.pca$x[,1],

data.pca$x[,2],

xlab=paste("PC1","(",pro1,"%)",sep=""),

ylab=paste("PC2","(",pro2,"%)",sep=""),

main="PCA",

xlim=c(1.1\*xmin,1.1\*xmax),

ylim=c(1.1\*ymin,1.1\*ymax),

pch=16,col=mycolors)

abline(h=0,col="gray")

abline(v=0,col="gray")

text(data.pca$x[,1],data.pca$x[,2],labels=samples)

dev.off()

**8. Barplot**

**8.1 One variable Barplot**

data = read.table("test.barplot.dat",header=T,row.names=1)

matrix\_data = as.matrix(data)

t\_matrix\_data = t(matrix\_data)

mycolors = c("#fb6a4a","#74c476")

barplot(t\_matrix\_data,col=mycolors,legend=c("Up","Down"))

barplot(t\_matrix\_data,col=mycolors,legend=c("Up","Down"),horiz=T,beside=T)

par(mar=c(5,7,4,2)+0.1)

barplot(t\_matrix\_data,col=mycolors,legend=c("Up","Down"),horiz=T,beside=T,las=1)

**8.2 Two variable Barplot**

library(reshape2)

dt <- data.frame(Species=c('M.honghuensis','T.kitauei','K.iwatai','E. leei','S. zaharoni','Hydra','Nematostella'),GC=c(27.4,37.5,28,33.5,28,29,41),Genome\_size=c(205.56,150.7,22.5,67.98,173.59,1005,450))

dt.long <- melt(dt)

ggplot(dt.long,aes(Species,value,fill=variable))+geom\_bar(stat="identity",position="dodge")+theme(axis.text.x=element\_text(angle=90))

**9. Histogram**

**9.1 Histogram for gene expression**

diff = read.table("gene\_exp.diff",header=T,sep="\t",comment.char="")

values = diff$log2FC

len = length(values)

for (i in 1 : len)

{

`	`if (values[i] > 10 || values[i] == "Inf")

`	`{

`		`values[i] = 10;

`	`}

`	`else if (values[i] < -10 || values[i] == "-Inf")

`	`{

`		`values[i] = -10;

`	`}

}

hist(values,col="lightgreen",breaks=seq(-10,10,0.5),freq=F)

lines(density(values),lwd=2,col="red")

**9.2 Histogram for sequence species distribution**

library(ggplot2)

library(ggthemes)

data <- read.csv("Species\_distribution.csv",header=T)

data <- data[seq(30),]

data <- data.frame(A=c(data[1]),B=c(data[2]),C=c(data[3]))

data$Species <- factor(data$Species, levels = data$Species[order(data$Count)])

p <- ggplot(data,aes(y=Count, x=Species))

p+geom\_bar(stat="identity", width=0.8, fill= "dodgerblue3",colour="black", cex=0.4)+coord\_flip()+ theme\_set(theme\_bw())+labs(title='The top hit species distribution')+theme(line=element\_line(size=0.5))+ ylab("BLAST Hit")+geom\_text(data=data, aes(x=Species, y=Count, label=paste0(round((Count/829)\*100, digits=2), "%")), hjust=-0.2 ,size=3)+ scale\_y\_continuous(limits=c(0,250))

ggsave("plot2.pdf", width=8, height=8)

**10. Heatmap**

library(gplots)

data=read.table("gene\_exp.xls", header=T, row.names=1)

dim(data)

tmp<-data[,1:11]

dim(tmp)

input = as.matrix(tmp)

\# Start drawing

heatmap.2(input)

\# col

heatmap.2(input,col=greenred(255))

\# No trace and density info

heatmap.2(input,col=greenred(255),trace="none",density.info="none")

heatmap.2(input,col=greenred(255),trace="none",density.info="none",scale="row")

heatmap.2(input,col=greenred(255),trace="none",density.info="none",scale="col")

\# No hierarchical clustering

heatmap.2(input,col=greenred(255),trace="none",density.info="none",scale="row",Rowv=FALSE,Colv=FALSE,dendrogram="none")

\# Only do hierarchical clustering for column

heatmap.2(input,col=greenred(255),trace="none",density.info="none",scale="row",Rowv=FALSE,dendrogram="column")

\# Only do hierarchical clustering for row

heatmap.2(input,col=greenred(255),trace="none",density.info="none",scale="row",Colv=FALSE,dendrogram="row")

\# Do hierarchical clustering for both row and column

heatmap.2(input,col=greenred(255),trace="none",density.info="none",scale="row",dendrogram="none")

\# margins

heatmap.2(input,col=greenred(255),trace="none",density.info="none",scale="row",margins=c(5,14))

\# RowSideColors, ColSideColors

colColors = c(rep("#fb6a4a",5),rep("#74c476",6))

heatmap.2(input,col=greenred(255),trace="none",density.info="none",scale="row",margins=c(5,14),ColSideColors=colColors)

\# Normalization on row

tmp1 = t(scale(t(tmp)))

heatmap.2(as.matrix(tmp1),col=greenred(255),trace="none",density.info="none",margins=c(5,14))

\# Normalization on column

tmp2 = scale(tmp)

heatmap.2(as.matrix(tmp2),col=greenred(255),trace="none",density.info="none",margin=c(5,14))

\## Normalization on both

\# Change the matrix into vector

tlist=c()

for(i in 1:nrow(tmp))

{

`	`tlist=append(tlist,as.numeric(tmp[i,]))

}

\# Normalization

cl=scale(tlist)

\# Change the vector into matrix

tmp3 = matrix(cl,nrow=nrow(tmp),byrow=T)

rownames(tmp3)=rownames(tmp)

colnames(tmp3)=colnames(tmp)

\# Drawing

heatmap.2(as.matrix(tmp3),col=greenred(255),trace="none",density.info="none",margin=c(5,14))

**11. Venn diagram**

aov = read.table("aov.out.diff.xls",header=T,row.names=1)

#aovpass = rownames(aov)[aov$Pvalue<=0.05]

aovpass = rownames(aov)

\# ttest

ttest = read.table("ttest.out.diff.xls",header=T,row.names=1)

#ttestpass = rownames(ttest)[ttest$Pvalue<=0.05]

ttestpass = rownames(ttest)

\# wlicox 

wlicox = read.table("wilcox.out.diff.xls",header=T,row.names=1)

#wlicoxpass = rownames(wlicox)[wlicox$Pvalue<=0.05]

wlicoxpass = rownames(wlicox)

\# cuffdiff

cuffdiff = read.table("gene\_exp.diff",header=T,row.names=1,sep="\t",comment.char="")

cuffdiffpass = rownames(cuffdiff)[cuffdiff$significant=="yes" & abs(cuffdiff$log2FC)>1]

library(VennDiagram)

venn.diagram(list(aov=aovpass,ttest=ttestpass,wlicox=wlicoxpass,cuffdiff=cuffdiffpass),

fill=c("red","green","blue","yellow"),cex=2, filename="DEGs.VennDiagram.tiff")

**12. Scatter graph**

expr = read.table("all.fpkm", header=T, row.names=1)

\# draw plot

png("C1\_C2.plot.png", width=480, height=480)

plot(log(expr$C1\_FPKM), log(expr$C2\_FPKM), pch=16,cex=0.5, col="red", xlim=c(-5,10),ylim=c(-5,10), xlab="log10(C1\_FPKM)", ylab="log10(C2\_FPKM)")

dev.off()

**License**

All source code, i.e. scripts/\*.pl, scripts/\*.sh or scripts/\*.py are under the MIT license.
