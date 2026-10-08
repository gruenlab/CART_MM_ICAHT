reticulate::use_python("/opt/python/bin/python",required=T)
modules <- reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")
reticulate::py_available()
require(RaceID)
require(Matrix)


sc<-readRDS("PATH/TO/SC_OBJ")
meta.data<-readRDS("PATH/TO/METADATA")
meta.data$broad_ct

library(stringr)
meta.data$broad_ct <-sub("Immune_","",meta.data$broad_ct)
meta.data$broad_ct <-sub("Niche_","",meta.data$broad_ct)

unique(meta.data$broad_ct)

cDB<-readRDS("PATH/TO/CELLPHONEDB")

#CellChatDB
library(CellChat)

DB <- CellChatDB.human

f <- DB$interaction$interaction_name_2
h <- strsplit(f," - ")
h_L <- c()
h_R <- c()
H_L <- c()
H_R <- c()
for (i in 1:length(h)){
  h_L[i] <- h[[i]][1]
  h_R[i] <- h[[i]][2]
  
  n <- strsplit(h_L[i],"")[[1]] == " "
  if ( sum(n) > 0 ){  for ( j in 1:sum( n ) )  h_L[i] <- sub("\\s+", "", h_L[i]) }
  h_L[i] <- sub("\\(", "", h_L[i])
  h_L[i] <- sub("\\)", "", h_L[i])
  
  n <- strsplit(h_R[i],"")[[1]] == " "
  if ( sum(n) > 0 ){  for ( j in 1:sum( n ) )  h_R[i] <- sub("\\s+", "", h_R[i]) }
  h_R[i] <- sub("\\(", "", h_R[i])
  h_R[i] <- sub("\\)", "", h_R[i])
  
  lv <- strsplit(h_L[i],"\\+")[[1]]
  rv <- strsplit(h_R[i],"\\+")[[1]]
  
  for ( il in 1:length(lv) ){
    for ( ir in 1:length(rv) ){
      H_L <- c(H_L,lv[il])
      H_R <- c(H_R,rv[ir])
    }
  }
}

x <- unique(data.frame(Ligand=H_L,receptor=H_R))
colnames(cDB)<-colnames(x)
x <- unique( rbind(cDB,x) )

## subsetting of cells to be included (if needed), e.g., when analyzing contributions of different samples to the same clustering
###subset cytopenic cells here and re-run the script below to generate another dfa and dff
f<-meta.data$cells[meta.data$cyt%in%"No_Cytopenia"]
y<-meta.data$broad_ct[meta.data$cyt%in%"No_Cytopenia"]
names(y)<-meta.data$cells[meta.data$cyt%in%"No_Cytopenia"]
head(y)
## expression matrix
z <- as.matrix(sc@ndata[,f])

## discard ligand and receptor pairs with ligand or receptor not contained in expression matrix
f <- x[,1] %in% rownames(z) & x[,2] %in% rownames(z)
xf <- x[f,]

## derive enriched ligand-receptor pairs
p2  <- as.vector(aggregate(rep(1,length(y)),by=list(y),sum)[,-1])
cl  <- sort(unique(y))
cl1 <- as.vector(t(sapply(cl,function(x){rep(x,length(cl))})))
cl2 <- as.vector(sapply(cl,function(x){rep(x,length(cl))}))    
N   <- sum(p2)
## p2: number of cells for each cluster
## N: total number of cells
thr <- -Inf # keep pairs with enrichment score > .9 (max. enrichment score is 1); keep all: -Inf
flag <- TRUE
for ( i in 1:nrow(xf) ){
  #n1 <- "Ccl4"
  #n2 <- "Ackr2"
  
  n1 <- xf[i,1]
  n2 <- xf[i,2]
  a1 <- aggregate(t(z[c(n1,n2),])>0,by=list(y),sum)
  rn <- a1[,1]
  ## p1: number of cells positive for ligand (column 1) or receptor (column 2) for each cluster
  p1 <- a1[,-1]
  
  pv1 <- apply(cbind(p1[,1],p2),1,function(x){ fisher.test(matrix(c(x[1],x[2] - x[1],sum(p1[,1]),N - sum(p1[,1])),ncol=2),alternative = "g")$p.value })
  pv2 <- apply(cbind(p1[,2],p2),1,function(x){ fisher.test(matrix(c(x[1],x[2] - x[1],sum(p1[,2]),N - sum(p1[,2])),ncol=2),alternative = "g")$p.value })
  
  pv1 <- p.adjust(pv1,method="bonferroni")
  pv2 <- p.adjust(pv2,method="bonferroni")
  fr <- p1/p2
  p <- (1 - pv1)%*%t(1-pv2)
  ind <- which(p > thr)
  
  
  ##  rows: cl1 columns: cl2
  ## n1g: fraction of ligand-expressing cells in each cluster
  n1g <- p1[,1]/p2
  ## p1g: fold change of the fraction of ligand-expressing cells in each cluster compared to the average across all clusters (baseline)
  p1g <- n1g/mean(n1g)
  p1gs <- n1g/sum(n1g)
  
  v <-  p1gs*log(p1gs)/log(length(p1gs)) 
  E1g <- -sum(v[ p1gs >0 ])
  ## n2g: fraction of receptor-expressing cells in each cluster
  n2g <- p1[,2]/p2
  ## p2g: fold change of the fraction of receptor-expressing cells in each cluster compared to the average across all clusters (baseline)
  p2g <- n2g/mean(n2g)
  p2gs <- n2g/sum(n2g)
  v <-  p2gs*log(p2gs)/log(length(p2gs)) 
  E2g <- -sum(v[ p2gs >0 ])
  
  
  rownames(fr) <- names(p1g) <- names(p2g) <- rn 
  
  if (length(ind) > 0){
    resN <- cbind(rep(n1,length(ind)),rep(n2,length(ind)))
    c1i <- as.character(cl1[ind])
    c2i <- as.character(cl2[ind])
    res <- cbind(cl1[ind],cl2[ind],as.vector(p)[ind],fr[c1i,1],fr[c2i,2],p1g[c1i],p2g[c2i],.5*(p1g[c1i]+p2g[c2i]),rep(E1g,length(ind)),rep(E2g,length(ind)),.5*rep(E1g+E2g,length(ind)))
    colnames(resN) <- c("n1","n2")
    colnames(res) <- c("cl1","cl2","score","fr1","fr2","p1g","p2g","p12g","E1g","E2g","E12g")
    
    if ( flag ){
      int <- res
      intN <- resN
      flag <- FALSE
    }else{
      int <- rbind(int,res)
      intN <- rbind(intN,resN)
    }
  }
}

## data.frame with enriched ligand-receptor pairs
## columns:
## n1 ligand gene symbol
## n2 receptor gene symbol
## cl1 cluster enriched in ligand-expressing cells
## cl2 cluster enriched in receptor-expressing cells
## score interaction score  = (1 - pv(ligand))*(1 - pv(receptor))
##       where pv(ligand) is the hypergeometric p-value for overrepresentation of ligand-expressing cells in cl1 and
##       and pv(receptor) is the hypergeometric p-value for overrepresentation of receptor-expressing cells in cl2
##     The maximum score is 1 and the minimum score is 0
## fr1 fraction of ligand-expressing cells in cl1
## fr2 fraction of receptor-expressing cells in cl2
## p1g enrichment of ligand-expressing cells in cl1 
## p2g enrichment of receptor-expressing cells in cl2 
## p12g enrichment score for ligand-receptor pair in cl1 and cl2 ( = 1/2*(p1g + p2g) )
## E1g entropy of ligand across clusters (min:0 max:1)
## E2g entropy of receptor across clusters (min:0 max:1)
## E12g entropy of ligand-receptor pair across clusters ( = 1/2*(E1g + E2g) ); small values correspond to cluster-specific expression

## Remark: to obtain predictions with cluster-specific expression, filter by p12g

d <- unique(data.frame(intN,int))
rownames(d) <- 1:nrow(d)
f <- apply(is.na(d),1,sum) == 0
d <- d[f,]


## The ouput can be computed for many sample; initialize a list:
df <- list()
df[["WT"]] <- d
## add more samples if needed


### define cell types (here, in this example: cell type equals clutser number)
CT <- list()
for ( i in sort(unique(y)) ){ CT[[i]] <- i }
#names(CT) <- paste("cl",sort(unique(y)),sep=".")
colCT <- sc@fcol
## here: cell types defined by names (example with RaceID v0.2.3)



CT2 <- CT
a1 <- a2 <- c()
for ( i in names(CT2) ){
  a1 <- c(a1, rep(i,length(CT2[[i]])) )
  a2 <- c(a2,CT2[[i]])
}
ordered_cell_pops <- a1[order(a2)]
names(ordered_cell_pops) <- paste("cl",a2[order(a2)],sep=".")

## make new list with cell type info added
dfa <- df
for ( i in names(df) ){
  dfa[[i]][,"ct1"] <- ordered_cell_pops[paste("cl",dfa[[i]]$cl1,sep=".")]
  dfa[[i]][,"ct2"] <- ordered_cell_pops[paste("cl",dfa[[i]]$cl2,sep=".")]
}


## make filtered list
## CHANGE FILTERING PARAMETERS FOR DESIRED OUTPUT
dff <- list()
for ( i in 1:length(dfa) ){
  samp <- names(dfa)[i]
  tmp <- dfa[[samp]]
  ## here...
  f <- tmp$score > 0 & !is.na(tmp$fr1) & !is.na(tmp$fr2) & tmp$score > 0 & tmp$fr1 > 0 & tmp$fr2 >0& tmp$p12g >0
  dff[[samp]] <- tmp[f,]
}


## write output
out <- dff
for ( i in names(out) ){
  x <- out[[i]][,c("n1","n2","ct1","ct2","score","fr1","fr2","cl1","cl2")]
  colnames(x) <- c("Ligand","Receptor","Cell Type Ligand","Cell Type Receptor","Score (max. equals 1)","Fraction Ligand+","Fraction Receptor+","cl1","cl2")
  out[[i]] <- x
}

out <- out[["WT"]]


# you can show the statistics for all cell type pairs if you run the pipeline without filtering (thr > -Inf), i.e., you keep all pairs in the initial object df


LRDotPlot <- function (l, cap = Inf, flo = -Inf){
  l <- l[order(l$ct2),]
  l <- l[order(l$ct1),]
  
  ct <-  paste(l$ct1,l$ct2,sep="_")
  lr <- paste(l$n1,l$n2,sep="_")
  
  
  data <- data.frame(LR = factor(lr), 
                     CellType = factor(ct), Enrichment = l$p12g, 
                     Score = l$score)
  data[which(data$Enrichment > cap), "Enrichment"] <- cap
  data[which(data$Enrichment < flo), "Enrichment"] <- flo
  colorPalette = c("darkblue", "blue", "grey", "red", "darkred")
  ColorRamp <- colorRampPalette(colorPalette)(100)
  
  print(ggplot(data, aes_string(x = "CellType", y = "LR")) + 
          geom_point(aes_string(size = "Enrichment", color = "Score")) + 
          scale_colour_gradientn(colours = ColorRamp) + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(axis.text.x = element_text(angle = 90,hjust = 1)))
  
}


dff_cyt <- list()
for ( i in 1:length(dfa_cyt) ){
  samp <- names(dfa_cyt)[i]
  tmp <- dfa_cyt[[samp]]
  f <- tmp$score > .9 & !is.na(tmp$fr1) & !is.na(tmp$fr2) & tmp$score > 0.9 & tmp$fr1 > .2 & tmp$fr2 > .2& tmp$p12g > 4
  dff_cyt[[samp]] <- tmp[f,]
}

dff_nocyt <- list()
for ( i in 1:length(dfa_nocyt) ){
  samp <- names(dfa_nocyt)[i]
  tmp <- dfa_nocyt[[samp]]
  f <- tmp$score > .9 & !is.na(tmp$fr1) & !is.na(tmp$fr2) & tmp$score > 0.9 & tmp$fr1 > .2 & tmp$fr2 > .2& tmp$p12g > 4
  dff_nocyt[[samp]] <- tmp[f,]
}


LRDotPlot <- function (l, cap = Inf, flo = -Inf){
  l <- l[order(l$ct2),]
  l <- l[order(l$ct1),]
  
  ct <-  paste(l$ct1,l$ct2,sep="_")
  lr <- paste(l$n1,l$n2,sep="_")
  
  
  data <- data.frame(LR = factor(lr), 
                     CellType = factor(ct), Enrichment = l$p12g, 
                     Score = l$score,cyt=l$cyt)
  data[which(data$Enrichment > cap), "Enrichment"] <- cap
  data[which(data$Enrichment < flo), "Enrichment"] <- flo
  colorPalette = c("darkblue", "blue", "grey", "red", "darkred")
  ColorRamp <- colorRampPalette(colorPalette)(100)
  
  print(ggplot(data, aes_string(x = "CellType", y = "LR")) + 
          geom_point(aes_string(size = "Enrichment", color = "Score")) + facet_wrap(~ cyt)+
          scale_colour_gradientn(colours = ColorRamp) + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(axis.text.x = element_text(angle = 90,hjust = 1)))
  
}

dff_nocyt[[1]]$ct1<-sub("cl.","",dff_nocyt[[1]]$ct1)
dff_nocyt[[1]]$ct2<-sub("cl.","",dff_nocyt[[1]]$ct2)
dff_cyt[[1]]$ct1<-sub("cl.","",dff_cyt[[1]]$ct1)
dff_cyt[[1]]$ct2<-sub("cl.","",dff_cyt[[1]]$ct2)


l_nocyt <- dff_nocyt[[1]][dff_nocyt[[1]]$ct1%in%c("T.cell","pDC","Monocyte","Macrophage","Neutrophil" )&dff_nocyt[[1]]$ct2%in%c("Infl.Mesenchymal","Myeloma.mesenchymal","Mesenchymal" ) & !dff_nocyt[[1]]$ct1%in%c("Infl.Mesenchymal","Myeloma.mesenchymal","Mesenchymal" ),]
l_nocyt$cyt <- "No_Cyt"
l_cyt <- dff_cyt[[1]][dff_cyt[[1]]$ct1%in%c("T.cell","pDC","Mocyte","Macrophage","Neutrophil" )&dff_cyt[[1]]$ct2%in%c("Infl.Mesenchymal","Myeloma.mesenchymal","Mesenchymal" ) & !dff_cyt[[1]]$ct1%in%c("Infl.Mesenchymal","Myeloma.mesenchymal","Mesenchymal" ),]
l_cyt$cyt <- "Yes_Cyt"

ll_nocyt <- dfa_nocyt[[1]][dfa_nocyt[[1]]$ct1%in%c("T.cell","pDC","Monocyte","Macrophage","Neutrophil"  )&dfa_nocyt[[1]]$ct2%in%c("Infl.Mesenchymal","Myeloma.mesenchymal","Mesenchymal" )& !dfa_nocyt[[1]]$ct1%in%c("Infl.Mesenchymal","Myeloma.mesenchymal","Mesenchymal" ),]
ll_nocyt$cyt <- "No_Cyt"

ll_cyt <- dfa_cyt[[1]][dfa_cyt[[1]]$ct1%in%c("T.cell","pDC","Monocyte","Macrophage","Neutrophil"  )&dfa_cyt[[1]]$ct2%in%c("Infl.Mesenchymal","Myeloma.mesenchymal","Mesenchymal" )& !dfa_cyt[[1]]$ct1%in%c("Infl.Mesenchymal","Myeloma.mesenchymal","Mesenchymal" ),]
ll_cyt$cyt <- "Yes_Cyt"
ll<-rbind(ll_nocyt,ll_cyt)
l<-rbind(l_nocyt,l_cyt)

ll$n <- paste(ll$n1,ll$n2,sep="_")
l$n <- paste(l$n1,l$n2,sep="_")
f <- ll$n %in% unique(l$n)
ll <- ll[f,]
ll$e<-paste(ll$ct1,ll$ct2,ll$cyt,ll$n,sep = "_")
ll <- ll[!duplicated(ll$e), ]

LRDotPlot(ll,flo=3,cap=4)


