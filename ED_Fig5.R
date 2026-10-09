reticulate::use_python("/opt/python/bin/python",required=T)
modules <- reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")
reticulate::py_available()
require(RaceID)
require(Matrix)
require(readr)
library(stringr)

sc<-readRDS("PATH/TO/SC_OBJ")
meta.data<-readRDS("PATH/TO/supple_table8")
Tcell <-meta.data$cells[grep("T.cell",meta.data$cell_type)]
prdata_t<-sc@expdata[,Tcell]

sc <-SCseq(prdata_t)
sc<-filterdata(sc,mintotal = 1000,CGenes=rownames(sc@expdata)[grep("^(MT|RP(L|S)|GM\\D|GYPA)",rownames(sc@expdata))])
expData  <- getExpData(sc)
res   <- pruneKnn(expData,large=TRUE,regNB=TRUE,knn=25,seed=12345, no_cores=32,do.prune = FALSE)
cl    <- graphCluster(res,pvalue=0.01, use.leiden = T, leiden.resolution=1)
sc <- updateSC(sc,res=res,cl=cl,flo=.1)
sc <- compumap(sc)
plotmap(sc, um=T, cex=0.5)

C <- names(sc@cpart)[sc@cpart%in%c(15,16,17,18)] 
TF<-!colnames(sc@ndata)%in%C
prdata_clean<-prdata_t[,TF]

sc <-SCseq(prdata_clean)
sc<-filterdata(sc,mintotal = 1000,CGenes=rownames(sc@expdata)[grep("^(MT|RP(L|S)|GM\\D|GYPA)",rownames(sc@expdata))])
expData  <- getExpData(sc)
res   <- pruneKnn(expData,large=TRUE,regNB=TRUE,knn=25,seed=12345, no_cores=32,do.prune = FALSE)
cl    <- graphCluster(res,pvalue=0.01, use.leiden = T, leiden.resolution=1)
sc <- updateSC(sc,res=res,cl=cl,flo=.1)
sc <- compumap(sc)
plotmap(sc, um=T, cex=0.5)

C <- names(sc@cpart)[sc@cpart%in%c(15)] 
TF<-!colnames(sc@ndata)%in%C
prdata_clean1<-prdata_clean[,TF]

sc <-SCseq(prdata_clean1)
sc<-filterdata(sc,mintotal = 1000,CGenes=rownames(sc@expdata)[grep("^(MT|RP(L|S)|GM\\D|GYPA)",rownames(sc@expdata))])
expData  <- getExpData(sc)
res   <- pruneKnn(expData,large=TRUE,regNB=TRUE,knn=25,seed=12345, no_cores=32,do.prune = FALSE)
cl    <- graphCluster(res,pvalue=0.01, use.leiden = T, leiden.resolution=1)
sc <- updateSC(sc,res=res,cl=cl,flo=.1)
sc <- compumap(sc)
plotmap(sc, um=T, cex=0.5)

CD3e<-sc@expdata["CD3E",]
CD4<-sc@expdata["CD4",]
CD8A<-sc@expdata["CD8A",]
NK<-sc@expdata["KLRC1",]
FOXP3<-sc@expdata["FOXP3",]

CD3e<-names(CD3e)[CD3e>mean(CD3e)]
CD4<-names(CD4)[CD4>mean(CD4)]
CD8A<-names(CD8A)[CD8A>mean(CD8A)]
NK<-names(NK)[NK>mean(NK)]
reg<-names(FOXP3)[FOXP3>mean(FOXP3)]

CD4_T<-intersect(CD3e,CD4)
CD8A_T<-intersect(CD3e,CD8A)
NK_T<-intersect(CD3e,NK)
CD8_NKT<-intersect(NK_T,CD8A)
CD4_NKT<-intersect(NK_T,CD4)

meta.data_T[meta.data_T$cells%in%CD8A,]$broad_ct<-"CD8_T"
meta.data_T[meta.data_T$cells%in%CD4,]$broad_ct<-"CD4_T"
meta.data_T[meta.data_T$cells%in%reg,]$broad_ct<-"Treg"
meta.data_T[meta.data_T$cells%in%NK,]$broad_ct<-"NK"
meta.data_T[meta.data_T$cells%in%CD8_NKT,]$broad_ct<-"CD8_NKT"
m <- meta.data_T[!meta.data_T$broad_ct%in%"Immune_T.cell",]

####ED5a
ggplot(m,
       aes(x=sample, fill=broad_ct))+ theme_classic()+ geom_bar(position = "fill") +   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


####CAR T 
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

####ED5b
plotsymbolsmap(sc,types = meta.data$cell_type,um = T,map = TRUE, leg = F)
plotsymbolsmap(sc,types = meta.data$cell_type,um = T,map = F, leg = T)

####ED5c, example for patient NC12
plotexpmap(sc,"Cilta_Braun_CAR_construct",logsc=T,fr=F, um=T, cex=1,
           cells = colnames(sc@ndata)[grep("NC12",colnames(sc@ndata))])

####ED5d D90_clinical_score vs CAR T frequency in blood 8 weeks post CAR
a <- read.csv("PATH/TO/SUPPL_TABLE7")
b <- a[,c("Hb", "PLT",  "ANC")]
bn <- t(t(b)/colSums(b))

gmean <- function(x, na.rm = FALSE) {
  n <- if (na.rm) sum(!is.na(x)) else length(x)
  prod(x, na.rm = na.rm)^(1/n)
}

bng <- apply(bn,1,gmean)
v <- log(bng)
u <- log(a$X8W_postCAR_CD3)

fit <- lm(v ~ u)
model_summary <- tidy(fit)
p_value <- model_summary$p.value[2]
uv <- data_frame(clinical_score = v, X8W_postCAR_CD31 = u)

library(viridis)
ggplot(uv, 
       aes(x = X8W_postCAR_CD31, y = clinical_score)) +
  geom_point(color = "grey30") +  
  theme_classic() +
  geom_smooth(method = lm, color = "#1f78b4", fill = "#a6cee3") +  
  scale_fill_viridis(discrete = TRUE) +
  sm_statCorr(color = "#1f78b4", fill = "#a6cee3", 
              text_size = 5) +
  sm_classic()

####ED5d HSC senescence vs CAR T frequency in blood 8 weeks post CAR

v <- log(a$HSC_sen + 1)
u <- log(a$X8W_postCAR_CD3)

fit <- lm(v ~ u)
model_summary <- tidy(fit)
p_value <- model_summary$p.value[2]
uv <- data_frame(HSC_sen = v, X8W_postCAR_CD31 = u)

library(viridis)
ggplot(uv, 
       aes(x = X8W_postCAR_CD31, y = HSC_sen)) +
  geom_point(color = "grey30") +  
  theme_classic() +
  geom_smooth(method = lm, color = "#1f78b4", fill = "#a6cee3") +  
  scale_fill_viridis(discrete = TRUE) +
  sm_statCorr(color = "#1f78b4", fill = "#a6cee3", 
              text_size = 5) +
  sm_classic()

####ED5d Mesenchymal inflammation vs CAR T frequency in blood 8 weeks post CAR

v <- log(a$Mes_infl + 1)
u <- log(a$X8W_postCAR_CD3)

fit <- lm(v ~ u)
model_summary <- tidy(fit)
p_value <- model_summary$p.value[2]
uv <- data_frame(Mes_infl = v, X8W_postCAR_CD31 = u)

library(viridis)
ggplot(uv, 
       aes(x = X8W_postCAR_CD31, y = Mes_infl)) +
  geom_point(color = "grey30") +  
  theme_classic() +
  geom_smooth(method = lm, color = "#1f78b4", fill = "#a6cee3") +  
  scale_fill_viridis(discrete = TRUE) +
  sm_statCorr(color = "#1f78b4", fill = "#a6cee3", 
              text_size = 5) +
  sm_classic()

####ED5e D90_clinical_score vs CAR T frequency in blood 12 weeks post CAR
b <- a[,c("Hb", "PLT",  "ANC")]
bn <- t(t(b)/colSums(b))

gmean <- function(x, na.rm = FALSE) {
  n <- if (na.rm) sum(!is.na(x)) else length(x)
  prod(x, na.rm = na.rm)^(1/n)
}

bng <- apply(bn,1,gmean)
v <- log(bng)
u <- log(a$X12W_postCAR_CD3)

fit <- lm(v ~ u)
model_summary <- tidy(fit)
p_value <- model_summary$p.value[2]
uv <- data_frame(clinical_score = v, X12W_postCAR_CD3 = u)

library(viridis)
ggplot(uv, 
       aes(x = X12W_postCAR_CD3, y = clinical_score)) +
  geom_point(color = "grey30") +  
  theme_classic() +
  geom_smooth(method = lm, color = "#1f78b4", fill = "#a6cee3") +  
  scale_fill_viridis(discrete = TRUE) +
  sm_statCorr(color = "#1f78b4", fill = "#a6cee3", 
              text_size = 5) +
  sm_classic()

####ED5d HSC senescence vs CAR T frequency in blood 8 weeks post CAR
v <- log(a$HSC_sen + 1)
u <- log(a$X12W_postCAR_CD3)

fit <- lm(v ~ u)
model_summary <- tidy(fit)
p_value <- model_summary$p.value[2]
uv <- data_frame(HSC_sen = v, X12W_postCAR_CD3 = u)

library(viridis)
ggplot(uv, 
       aes(x = X12W_postCAR_CD3, y = HSC_sen)) +
  geom_point(color = "grey30") +  
  theme_classic() +
  geom_smooth(method = lm, color = "#1f78b4", fill = "#a6cee3") +  
  scale_fill_viridis(discrete = TRUE) +
  sm_statCorr(color = "#1f78b4", fill = "#a6cee3", 
              text_size = 5) +
  sm_classic()


####ED5e D90_clinical_score vs CAR T frequency in BM 12 weeks post CAR
b <- a[,c("Hb", "PLT",  "ANC")]
bn <- t(t(b)/colSums(b))

gmean <- function(x, na.rm = FALSE) {
  n <- if (na.rm) sum(!is.na(x)) else length(x)
  prod(x, na.rm = na.rm)^(1/n)
}

bng <- apply(bn,1,gmean)
v <- log(bng)
u <- log(a$CART_IF_stain)

fit <- lm(v ~ u)
model_summary <- tidy(fit)
p_value <- model_summary$p.value[2]
uv <- data_frame(clinical_score = v, CART_IF_stain = u)

library(viridis)
ggplot(uv, 
       aes(x = CART_IF_stain, y = clinical_score)) +
  geom_point(color = "grey30") +  
  theme_classic() +
  geom_smooth(method = lm, color = "#1f78b4", fill = "#a6cee3") +  
  scale_fill_viridis(discrete = TRUE) +
  sm_statCorr(color = "#1f78b4", fill = "#a6cee3", 
              text_size = 5) +
  sm_classic()


