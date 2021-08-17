if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install(c("ggplot2", "tidyverse", "BBmisc", "limma"))

library(ggplot2)
library(tidyverse)
library('BBmisc')
library(limma)

setwd("PROJECT FOLDER")

tablet<-read.csv('covidData.txt') #read table from txt file
normalized<-normalize(tablet) # normalize data to 0-1

value<-c() 
name<-c()
for (t in 1:10) {
  value<-append(value,as.numeric(normalized[t,])) #get line values by 1-10 sample id
  name<-append(name,rep(as.character(t),3001))
  
}

data <- data.frame(name,value)
data %>%                        # boxplot using ggplot2
  ggplot( aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9)

group <- factor(c('treat','treat','control','control', 'treat','control','treat','treat','control','control'), levels = c('control', 'treat'))
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(group)
exprSet<-read.csv('data.csv')
rownames(exprSet)<-exprSet[,1] # dataframe modify
exprSet<-exprSet[,-1] 
fit <- lmFit(exprSet, design, method = 'ls')
contrast <- makeContrasts('treat-control', levels = design)


#Bayes model
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

qqt(fit2$t, df = fit2$df.prior+fit2$df.residual, pch = 16, cex = 0.2)
abline(0,1)

#pvalue
diff_gene <- topTable(fit2, number = Inf, adjust.method = 'fdr')
head(diff_gene, 100)


diff_gene[which(diff_gene$logFC >= 1 & diff_gene$adj.P.Val < 0.01),'sig'] <- 'red'
diff_gene[which(diff_gene$logFC <= -1 & diff_gene$adj.P.Val < 0.01),'sig'] <- 'blue'
diff_gene[which(abs(diff_gene$logFC) < 1 | diff_gene$adj.P.Val >= 0.01),'sig'] <- 'gray'

log2FoldChange <- diff_gene$logFC
FDR <- diff_gene$adj.P.Val
plot(log2FoldChange, -log10(FDR), pch = 20, col = diff_gene$sig) #volcano plot
abline(v = 1, lty = 2)
abline(v = -1, lty = 2)
abline(h = -log(0.01, 10), lty = 2)

