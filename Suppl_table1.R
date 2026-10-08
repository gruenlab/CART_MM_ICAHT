reticulate::use_python("/opt/python/bin/python",required=T)
modules <- reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")
reticulate::py_available()
require(RaceID)
require(Matrix)

sc<-readRDS("PATH/TO/SC_OBJ")
meta.data<-readRDS("PATH/TO/METADATA")
plotmap(sc, um=T, cex=0.5)

a <- setdiff(rownames(sc@expdata),sc@genes)

unique(meta.data$cell_type)
B <- meta.data[meta.data$cell_type%in%"Infl.mesenchymal",]$cells
A <- meta.data[!meta.data$cell_type%in%A,]$cells

x <- diffexpnb(getfdata(sc,n=c(A,B),g=c(sc@genes,"CXCL12","LEPR")), A=A, B=B)
infl_vs_all <- x$res


B <- meta.data[meta.data$cell_type%in%c("Mesenchymal","Infl.mesenchymal","Myeloma.mesenchymal"),]$cells
A <- meta.data[!meta.data$cell_type%in%A,]$cells

x <- diffexpnb(getfdata(sc,n=c(A,B),g=c(sc@genes,"CXCL12","LEPR")), A=A, B=B )
allmes_vs_all <- x$res

B <- meta.data[meta.data$cell_type%in%c("Mesenchymal"),]$cells
A <- meta.data[!meta.data$cell_type%in%A,]$cells

x <- diffexpnb(getfdata(sc,n=c(A,B),g=c(sc@genes,"CXCL12","LEPR")), A=A, B=B )
mes_vs_all <- x$res

B <- meta.data[meta.data$cell_type%in%"Myeloma.mesenchymal",]$cells
A <- meta.data[!meta.data$cell_type%in%A,]$cells

x <- diffexpnb(getfdata(sc,n=c(A,B),g=c(sc@genes,"CXCL12","LEPR")), A=A, B=B )
mye_vs_all <- x$res

