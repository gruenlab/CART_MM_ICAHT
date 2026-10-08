sc<-readRDS("PATH/TO/SC_OBJ")
meta.data<-readRDS("PATH/TO/METADATA")

exp<-sc@expdata

umi<-colSums(exp)
head(umi)
head(sc@expdata)

exgene <- exp != 0
exgene[1:5,1:5]
ngene<-colSums(exgene)


mtratio<- colSums(exp[grep("^MT",rownames(exp)),])/umi

meta.data$umi<-umi
meta.data$ngene<-ngene
meta.data$mt<-mtratio

df<-data.frame(cells=meta.data$cells,umi=meta.data$umi,ngene=meta.data$ngene,mt=meta.data$mt)
library(reshape2)
df_long <- melt(df, id.vars = "cells", variable.name = "metric", value.name = "value")

###Suppl_Fig1a
p_mt <- ggplot(subset(df_long, metric == "mt"), aes(x = metric, y = value)) +
  geom_violin(fill = "pink", color = "darkred") +
  theme_minimal() +geom_boxplot(width=0.1, fill="white")+
  labs(title = "MT Violin Plot", x = "MT", y = "Value")


###Suppl_Fig1b
p_ngenes <- ggplot(subset(df_long, metric == "ngene"), aes(x = metric, y = value)) +
  geom_violin(fill = "lightgreen", color = "darkgreen") +
  theme_minimal() +geom_boxplot(width=0.1, fill="white")+
  labs(title = "nGenes Violin Plot", x = "nGenes", y = "Value")

###Suppl_Fig1c
p_umi <- ggplot(subset(df_long, metric == "umi"), aes(x = metric, y = value)) +
  geom_violin(fill = "skyblue", color = "darkblue") +
  theme_minimal() +geom_boxplot(width=0.1, fill="white")+
  labs(title = "UMI Violin Plot", x = "UMI", y = "Value")

###patient level
df<-data.frame(cells=meta.data$cells,samples=meta.data$sample,umi=meta.data$umi,ngene=meta.data$ngene,mt=meta.data$mt)
library(ggplot2)
library(reshape2)

df_long <- melt(df, id.vars = c("cells", "samples"), 
                variable.name = "metric", 
                value.name = "value")

###Suppl_Fig1d
p_mt <- ggplot(subset(df_long, metric == "mt"), aes(x = metric, y = value)) +
  geom_violin(fill = "pink", color = "darkred") +
  theme_minimal() +geom_boxplot(width=0.1, fill="white")+
  labs(title = "MT Violin Plot", x = "MT", y = "Value")

###Suppl_Fig1e
p_ngenes <- ggplot(subset(df_long, metric == "ngene"), aes(x = metric, y = value)) +
  geom_violin(fill = "lightgreen", color = "darkgreen") +
  theme_minimal() +geom_boxplot(width=0.1, fill="white")+
  labs(title = "nGenes Violin Plot", x = "nGenes", y = "Value")

###Suppl_Fig1f
p_umi <- ggplot(subset(df_long, metric == "umi"), aes(x = metric, y = value)) +
  geom_violin(fill = "skyblue", color = "darkblue") +
  theme_minimal() +geom_boxplot(width=0.1, fill="white")+
  labs(title = "UMI Violin Plot", x = "UMI", y = "Value")

