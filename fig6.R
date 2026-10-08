reticulate::use_python("/opt/python/bin/python",required=T)
modules <- reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")
reticulate::py_available()
require(RaceID)
require(Matrix)
sc<-readRDS("PATH/TO/SC_OBJ")
meta.data<-readRDS("PATH/TO/METADATA")

plotexpmap(sc,"IFNG",logsc=T,fr=F, um=T, cex=40,leg = FALSE,cells = meta.data[grep("No",meta.data$cyt),]$cells)
plotexpmap(sc,"IFNG",logsc=T,fr=F, um=T, cex=40,leg = TRUE,map = FALSE,cells = meta.data[grep("No",meta.data$cyt),]$cells)


plotexpmap(sc,"IFNG",logsc=T,fr=F, um=T, cex=40,leg = FALSE,cells = meta.data[grep("Yes",meta.data$cyt),]$cells)
plotexpmap(sc,"IFNG",logsc=T,fr=F, um=T, cex=40,leg = TRUE,map = FALSE,cells = meta.data[grep("Yes",meta.data$cyt),]$cells)


plotexpmap(sc,"AREG",logsc=T,fr=F, um=T, cex=40,leg = FALSE,cells = meta.data[grep("No",meta.data$cyt),]$cells)
plotexpmap(sc,"AREG",logsc=T,fr=F, um=T, cex=40,leg = TRUE,map = FALSE,cells = meta.data[grep("No",meta.data$cyt),]$cells)

plotexpmap(sc,"AREG",logsc=T,fr=F, um=T, cex=40,leg = FALSE,cells = meta.data[grep("Yes",meta.data$cyt),]$cells)
plotexpmap(sc,"AREG",logsc=T,fr=F, um=T, cex=40,leg = TRUE,map = FALSE,cells = meta.data[grep("Yes",meta.data$cyt),]$cells)

plotexpmap(sc,"EREG",logsc=T,fr=F, um=T, cex=40,leg = FALSE,cells = meta.data[grep("No",meta.data$cyt),]$cells)
plotexpmap(sc,"EREG",logsc=T,fr=F, um=T, cex=40,leg = TRUE,map = FALSE,cells = meta.data[grep("No",meta.data$cyt),]$cells)

plotexpmap(sc,"EREG",logsc=T,fr=F, um=T, cex=40,leg = FALSE,cells = meta.data[grep("Yes",meta.data$cyt),]$cells)
plotexpmap(sc,"EREG",logsc=T,fr=F, um=T, cex=40,leg = TRUE,map = FALSE,cells = meta.data[grep("Yes",meta.data$cyt),]$cells)

plotexpmap(sc,"TGFB1",logsc=T,fr=F, um=T, cex=40,leg = FALSE,cells = meta.data[grep("No",meta.data$cyt),]$cells)
plotexpmap(sc,"TGFB1",logsc=T,fr=F, um=T, cex=40,leg = TRUE,map = FALSE,cells = meta.data[grep("No",meta.data$cyt),]$cells)

plotexpmap(sc,"TGFB1",logsc=T,fr=F, um=T, cex=40,leg = FALSE,cells = meta.data[grep("Yes",meta.data$cyt),]$cells)
plotexpmap(sc,"TGFB1",logsc=T,fr=F, um=T, cex=40,leg = TRUE,map = FALSE,cells = meta.data[grep("Yes",meta.data$cyt),]$cells)

library(stringr)
meta.data$broad_ct <-sub("Immune_","",meta.data$broad_ct)
meta.data$broad_ct <-sub("Niche_","",meta.data$broad_ct)
meta.data$cells <-sub("/.*", "", meta.data$cells)


CAR<-sc@expdata["Cilta_Braun_CAR_construct",]
CAR<-names(CAR)[CAR>0]
meta.data[meta.data$cells%in%CAR,]$broad_ct<-"CAR.T"

cDB<-readRDS("PATH/to/CellphoneDB")

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

y<-meta.data$broad_ct
names(y)<-meta.data$cells
head(y)
## expression matrix
z <- as.matrix(sc@ndata)

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

dff <- list()
for ( i in 1:length(dfa) ){
  samp <- names(dfa)[i]
  tmp <- dfa[[samp]]
  f <- tmp$score > .9 & !is.na(tmp$fr1) & !is.na(tmp$fr2) & tmp$score > 0.9 & tmp$fr1 > .15 & tmp$fr2 > .15& tmp$p12g > 2
  dff[[samp]] <- tmp[f,]
}

unique(dff$WT$ct1)
l <- dff$WT[dff$WT$ct1%in%c("CAR.T","T.cell","pDC","Monocyte","Macrophage","Neutrophil" )&dff[[1]]$ct2%in%c("Infl.mesenchymal")& !dff[[1]]$ct1%in%c("Infl.mesenchymal"),]

library(circlize)
library(dplyr)
library(scales)

plotInteractionsSimple_Arrow <- function(x, pop_cols = NULL, 
                                         show_types = TRUE, show_legend = TRUE, 
                                         show_clust = TRUE, lab.cex = 0.5, 
                                         track.height = 0.1, legend.cex = .5) {
  library(circlize)
  
  # Ensure numeric
  x[,6] <- as.numeric(as.character(x[,6]))
  x[,7] <- as.numeric(as.character(x[,7]))
  
  # Ligand and receptor names
  lig_names <- paste(x[,3], x[,1], sep = "_")
  rec_names <- paste(x[,4], x[,2], sep = "_")
  
  all_lig <- unique(lig_names)
  all_rec <- unique(rec_names)
  
  # Interaction matrix
  m <- matrix(0, nrow = length(all_lig), ncol = length(all_rec))
  rownames(m) <- all_lig
  colnames(m) <- all_rec
  for (i in seq_len(nrow(x))) {
    m[lig_names[i], rec_names[i]] <- x[i,6] * x[i,7]
  }
  
  # Remove empty rows/cols
  m <- m[rowSums(m) > 0, colSums(m) > 0, drop = FALSE]
  
  # Colors
  all_cells <- sort(unique(c(x[,3], x[,4])))
  if (is.null(pop_cols)) {
    set.seed(123)
    pop_cols <- setNames(rainbow(length(all_cells)), all_cells)
  }
  
  lig_cell <- sapply(strsplit(rownames(m), "_"), `[`, 1)
  rec_cell <- sapply(strsplit(colnames(m), "_"), `[`, 1)
  sector_col <- c(pop_cols[lig_cell], pop_cols[rec_cell])
  names(sector_col) <- c(rownames(m), colnames(m))
  
  link_col <- pop_cols[lig_cell]
  link_col_mat <- matrix(rep(link_col, ncol(m)), nrow = nrow(m))
  
  circos.clear()
  
  # Correct gap.after based on actual number of sectors
  n_lig <- nrow(m)
  n_rec <- ncol(m)
  gap.after <- c(rep(1, n_lig - 1), 10, rep(1, n_rec - 1), 10)
  gap.after <- gap.after[1:(n_lig + n_rec)]  
  
  circos.par(
    start.degree = -5, 
    gap.after = gap.after, 
    cell.padding = c(0, 0, 0, 0), 
    points.overflow.warning = FALSE
  )
  
  chordDiagram(
    m, 
    col = link_col_mat, 
    grid.col = sector_col, 
    directional = 1,
    direction.type = c("arrows", "diffHeight"),
    annotationTrack = "grid", 
    preAllocateTracks = list(track.height = track.height)
  )
  
  # Labels
  if (show_clust) {
    circos.trackPlotRegion(
      track.index = 2, 
      panel.fun = function(x, y) {
        sn <- CELL_META$sector.index
        circos.text(
          CELL_META$xcenter, CELL_META$ycenter, 
          gsub("^[^_]+_", "", sn), 
          facing = "clockwise", niceFacing = TRUE, 
          adj = c(0, 0.5), cex = lab.cex
        )
      }, 
      bg.border = NA
    )
  }
  
  if (show_legend) {
    legend("topleft", legend = names(pop_cols), col = pop_cols, pch = 15, 
           cex = legend.cex, bty = "n")
  }
}


plotInteractionsSimple_Arrow(l)







