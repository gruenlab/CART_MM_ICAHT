reticulate::use_python("/opt/python/bin/python",required=T)
modules <- reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")
reticulate::py_available()
require(RaceID)
require(Matrix)
require(readr)
library(stringr)

x <- readMM("./C1_counts/result_dev.mtx")
f <- read.csv ("./C1_counts/result_dev.genes.txt", sep = "\t", header = FALSE)
b <- read.csv ("./C1_counts/result_dev.barcodes.txt", sep = "\t", header = FALSE)
xM <- as (x, "dgCMatrix")

dimnames(xM) <- list(as.character(b$V1), as.character(f$V1))



xM <- t(xM)
cs <- colSums(xM)
C1 <- xM [, cs>1000]

x <- readMM("./C2_counts/result_dev.mtx")
f <- read.csv ("./C2_counts/result_dev.genes.txt", sep = "\t", header = FALSE)
b <- read.csv ("./C2_counts/result_dev.barcodes.txt", sep = "\t", header = FALSE)
xM <- as (x, "dgCMatrix")

dimnames(xM) <- list(as.character(b$V1), as.character(f$V1))



xM <- t(xM)
cs <- colSums(xM)
C2 <- xM [, cs>1000]

x <- readMM("./C3_counts/result_dev.mtx")
f <- read.csv ("./C3_counts/result_dev.genes.txt", sep = "\t", header = FALSE)
b <- read.csv ("./C3_counts/result_dev.barcodes.txt", sep = "\t", header = FALSE)
xM <- as (x, "dgCMatrix")

dimnames(xM) <- list(as.character(b$V1), as.character(f$V1))



xM <- t(xM)
cs <- colSums(xM)
C3 <- xM [, cs>1000]

x <- readMM("./C4_counts/result_dev.mtx")
f <- read.csv ("./C4_counts/result_dev.genes.txt", sep = "\t", header = FALSE)
b <- read.csv ("./C4_counts/result_dev.barcodes.txt", sep = "\t", header = FALSE)
xM <- as (x, "dgCMatrix")

dimnames(xM) <- list(as.character(b$V1), as.character(f$V1))



xM <- t(xM)
cs <- colSums(xM)
C4 <- xM [, cs>1000]

x <- readMM("./NC1_counts/result_dev.mtx")
f <- read.csv ("./NC1_counts/result_dev.genes.txt", sep = "\t", header = FALSE)
b <- read.csv ("./NC1_counts/result_dev.barcodes.txt", sep = "\t", header = FALSE)
xM <- as (x, "dgCMatrix")

dimnames(xM) <- list(as.character(b$V1), as.character(f$V1))



xM <- t(xM)
cs <- colSums(xM)
NC1 <- xM [, cs>1000]

x <- readMM("./NC2_counts/result_dev.mtx")
f <- read.csv ("./NC2_counts/result_dev.genes.txt", sep = "\t", header = FALSE)
b <- read.csv ("./NC2_counts/result_dev.barcodes.txt", sep = "\t", header = FALSE)
xM <- as (x, "dgCMatrix")

dimnames(xM) <- list(as.character(b$V1), as.character(f$V1))



xM <- t(xM)
cs <- colSums(xM)
NC2 <- xM [, cs>1000]


x <- readMM("./NC3_counts/result_dev.mtx")
f <- read.csv ("./NC3_counts/result_dev.genes.txt", sep = "\t", header = FALSE)
b <- read.csv ("./NC3_counts/result_dev.barcodes.txt", sep = "\t", header = FALSE)
xM <- as (x, "dgCMatrix")

dimnames(xM) <- list(as.character(b$V1), as.character(f$V1))



xM <- t(xM)
cs <- colSums(xM)
NC3 <- xM [, cs>1000]

x <- readMM("./NC4_counts/result_dev.mtx")
f <- read.csv ("./NC4_counts/result_dev.genes.txt", sep = "\t", header = FALSE)
b <- read.csv ("./NC4_counts/result_dev.barcodes.txt", sep = "\t", header = FALSE)
xM <- as (x, "dgCMatrix")

dimnames(xM) <- list(as.character(b$V1), as.character(f$V1))



xM <- t(xM)
cs <- colSums(xM)
NC4 <- xM [, cs>1000]

x <- readMM("./NC5_counts/result_dev.mtx")
f <- read.csv ("./NC5_counts/result_dev.genes.txt", sep = "\t", header = FALSE)
b <- read.csv ("./NC5_counts/result_dev.barcodes.txt", sep = "\t", header = FALSE)
xM <- as (x, "dgCMatrix")

dimnames(xM) <- list(as.character(b$V1), as.character(f$V1))



xM <- t(xM)
cs <- colSums(xM)
NC5 <- xM [, cs>1000]

x <- readMM("./NC6_counts/result_dev.mtx")
f <- read.csv ("./NC6_counts/result_dev.genes.txt", sep = "\t", header = FALSE)
b <- read.csv ("./NC6_counts/result_dev.barcodes.txt", sep = "\t", header = FALSE)
xM <- as (x, "dgCMatrix")

dimnames(xM) <- list(as.character(b$V1), as.character(f$V1))



xM <- t(xM)
cs <- colSums(xM)
NC6 <- xM [, cs>1000]

x <- readMM("./NC7_counts/result_dev.mtx")
f <- read.csv ("./NC7_counts/result_dev.genes.txt", sep = "\t", header = FALSE)
b <- read.csv ("./NC7_counts/result_dev.barcodes.txt", sep = "\t", header = FALSE)
xM <- as (x, "dgCMatrix")

dimnames(xM) <- list(as.character(b$V1), as.character(f$V1))



xM <- t(xM)
cs <- colSums(xM)
NC7 <- xM [, cs>1000]

x <- readMM("./NC8_counts/result_dev.mtx")
f <- read.csv ("./NC8_counts/result_dev.genes.txt", sep = "\t", header = FALSE)
b <- read.csv ("./NC8_counts/result_dev.barcodes.txt", sep = "\t", header = FALSE)
xM <- as (x, "dgCMatrix")

dimnames(xM) <- list(as.character(b$V1), as.character(f$V1))



xM <- t(xM)
cs <- colSums(xM)
NC8 <- xM [, cs>1000]

x <- readMM("./NC9_counts/result_dev.mtx")
f <- read.csv ("./NC9_counts/result_dev.genes.txt", sep = "\t", header = FALSE)
b <- read.csv ("./NC9_counts/result_dev.barcodes.txt", sep = "\t", header = FALSE)
xM <- as (x, "dgCMatrix")

dimnames(xM) <- list(as.character(b$V1), as.character(f$V1))



xM <- t(xM)
cs <- colSums(xM)
NC9 <- xM [, cs>1000]

x <- readMM("./NC10_counts/result_dev.mtx")
f <- read.csv ("./NC10_counts/result_dev.genes.txt", sep = "\t", header = FALSE)
b <- read.csv ("./NC10_counts/result_dev.barcodes.txt", sep = "\t", header = FALSE)
xM <- as (x, "dgCMatrix")

dimnames(xM) <- list(as.character(b$V1), as.character(f$V1))



xM <- t(xM)
cs <- colSums(xM)
NC10 <- xM [, cs>1000]

x <- readMM("./NC11_counts/result_dev.mtx")
f <- read.csv ("./NC11_counts/result_dev.genes.txt", sep = "\t", header = FALSE)
b <- read.csv ("./NC11_counts/result_dev.barcodes.txt", sep = "\t", header = FALSE)
xM <- as (x, "dgCMatrix")

dimnames(xM) <- list(as.character(b$V1), as.character(f$V1))



xM <- t(xM)
cs <- colSums(xM)
NC11 <- xM [, cs>1000]

x <- readMM("./NC12_counts/result_dev.mtx")
f <- read.csv ("./NC12_counts/result_dev.genes.txt", sep = "\t", header = FALSE)
b <- read.csv ("./NC12_counts/result_dev.barcodes.txt", sep = "\t", header = FALSE)
xM <- as (x, "dgCMatrix")

dimnames(xM) <- list(as.character(b$V1), as.character(f$V1))



xM <- t(xM)
cs <- colSums(xM)
NC12 <- xM [, cs>1000]

x <- readMM("./P1_counts/result_dev.mtx")
f <- read.csv ("./P1_counts/result_dev.genes.txt", sep = "\t", header = FALSE)
b <- read.csv ("./P1_counts/result_dev.barcodes.txt", sep = "\t", header = FALSE)
xM <- as (x, "dgCMatrix")

dimnames(xM) <- list(as.character(b$V1), as.character(f$V1))



xM <- t(xM)
cs <- colSums(xM)
P1 <- xM [, cs>1000]

colnames(C1)<-paste(colnames(C1), "C1", sep = "/")
colnames(C2)<-paste(colnames(C2), "C2", sep = "/")
colnames(C3)<-paste(colnames(C3), "C3", sep = "/")
colnames(C4)<-paste(colnames(C4), "C4", sep = "/")
colnames(NC1)<-paste(colnames(NC1), "NC1", sep = "/")
colnames(NC2)<-paste(colnames(NC2), "NC2", sep = "/")
colnames(NC3)<-paste(colnames(NC3), "NC3", sep = "/")
colnames(NC4)<-paste(colnames(NC4), "NC4", sep = "/")
colnames(NC5)<-paste(colnames(NC5), "NC5", sep = "/")
colnames(NC6)<-paste(colnames(NC6), "NC6", sep = "/")
colnames(NC7)<-paste(colnames(NC7), "NC7", sep = "/")
colnames(NC8)<-paste(colnames(NC8), "NC8", sep = "/")
colnames(NC9)<-paste(colnames(NC9), "NC9", sep = "/")
colnames(NC10)<-paste(colnames(NC10), "NC10", sep = "/")
colnames(NC11)<-paste(colnames(NC11), "NC11", sep = "/")
colnames(NC12)<-paste(colnames(NC12), "NC12", sep = "/")
colnames(P1)<-paste(colnames(P1), "P1", sep = "/")

prdata <- cbind(C1,C2,C3,C4,NC1,NC2,NC3,NC4,NC5,NC6,NC7,NC8,NC9,NC10,NC11,NC12,P1)

sc <-SCseq(prdata)
sc<-filterdata(sc,mintotal = 1000,CGenes=rownames(sc@expdata)[grep("^(MT|RP(L|S)|GM\\D|GYPA)",rownames(sc@expdata))])
expData  <- getExpData(sc)
res   <- pruneKnn(expData,large=TRUE,regNB=TRUE,knn=25,seed=12345, no_cores=32,do.prune = FALSE,pcaComp = 40)
plotPC(res)

cl    <- graphCluster(res,pvalue=0.01, use.leiden = T, leiden.resolution=1)

sc <- updateSC(sc,res=res,cl=cl,flo=.1)
sc <- compumap(sc)
plotmap(sc, um=T, cex=0.5)

C <- names(sc@cpart)[sc@cpart%in%c(9,13)]
TF<-!colnames(sc@ndata)%in%C
prdata_clean<-prdata[,TF]

sc <-SCseq(prdata_clean)
sc<-filterdata(sc,mintotal = 1000,CGenes=rownames(sc@expdata)[grep("^(MT|RP(L|S)|GM\\D|GYPA)",rownames(sc@expdata))])
expData  <- getExpData(sc)
res   <- pruneKnn(expData,large=TRUE,regNB=TRUE,knn=25,seed=12345, no_cores=32,do.prune = FALSE,pcaComp = 40)
plotPC(res)

cl    <- graphCluster(res,pvalue=0.01, use.leiden = T, leiden.resolution=1)

sc <- updateSC(sc,res=res,cl=cl,flo=.1)
sc <- compumap(sc)
plotmap(sc, um=T, cex=0.5)

cluster_no <- sc@cpart
length(unique(cluster_no))

cluster_names <- c('B.prog', 'MDP/cDC', 'Mesenchymal', 'T.cell', 'B.prog', 'GMP',
                   'Erythroid', 'Monocyte', 'HSC', 'MEP/Megakaryocyte',
                   'Neutrophil', 'B.prog', 'NP', 'Erythroid', 'pDC', 
                   'B.cell', 'Erythroid', 'Neutrophil', 
                   'B.prog', 'Neutrophil', 'Neutrophil', 'Mast.cell'
                   , 'Inflammatory.mesenchymal', 'Myeloma.mesenchymal',
                   'Monocyte', 'NK', 'Endothelial', 'T.cell',
                   'Mast.cell', 'Macrophage', 'Osteolineage', 'SMC', "Erythroid"
                   
                   
)
length(cluster_names)

cluster_Ann <- cluster_no

value = 0
repeat{
  value = value + 1
  cluster_Ann <- replace(cluster_Ann,cluster_no==value, cluster_names[value])
  if (value == 33){
    print("repeat loop ends");
    break
  }
}
meta.data<-data.frame(cells=names(sc@cpart))
meta.data$sample<-sub(".+/", "", meta.data$cells)
meta.data$cell_type<-cluster_Ann


meta.data_imm <- meta.data[!meta.data$cell_type%in%c("Mesenchymal","SMC","Erythroid",
                                                     "Endothelial" ,"Inflammatory.mesenchymal",
                                                     "Myeloma.mesenchymal","Osteolineage" ),]


######

CAR<-sc@expdata["Cilta_Braun_CAR_construct",]
CAR<-names(CAR)[CAR>0]

meta.data_imm[meta.data_imm$cells%in%CAR,]$cell_type<-"CAR.T"
ifng <- sc@ndata["IFNG",meta.data_imm$cells]
df <- data.frame(cells = meta.data_imm$cells, cell_type = meta.data_imm$cell_type)
df$IFNG_expression <- ifng[match(df$cells, names(ifng))]

library(dplyr)

dot.df <- df %>%
  group_by(cell_type) %>%
  summarise(
    pct_expressing = mean(IFNG_expression > 0) * 100,
    mean_expression = mean(IFNG_expression),
    mean_expression_nonzero = mean(IFNG_expression[IFNG_expression > 0],
                                   na.rm = TRUE)
  )

dot.df


library(ggplot2)

###Fig7f
ggplot(dot.df,
       aes(x = "IFNG",
           y = cell_type,
           size = pct_expressing,
           color = mean_expression_nonzero)) +
  geom_point() +
  scale_size(range = c(4, 15)) +
  theme_classic() +
  labs(
    x = NULL,
    y = NULL,
    size = "% expressing",
    color = "Mean IFNG\n(non-zero cells)"
  )
