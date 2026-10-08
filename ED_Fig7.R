sc<-readRDS("PATH/TO/SC_OBJ")
meta.data<-readRDS("PATH/TO/METADATA")
reticulate::use_python("/opt/python/bin/python",required=T)
modules <- reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")
reticulate::py_available()
require(RaceID)
require(Matrix)

A <- intersect(meta.data$cells[grep("onocyte",meta.data$cell_type)],meta.data$cells[grep("No",meta.data$cyt)])
B <- intersect(meta.data$cells[grep("onocyte",meta.data$cell_type)],meta.data$cells[grep("Yes",meta.data$cyt)])

x <- diffexpnb(getfdata(sc,n=c(A,B)), A=A, B=B )
plotdiffgenesnb(x,pthr=.05,lthr=.5,mthr=-1,Aname="No_cytopenia",Bname="Yes_cytopenia",show_names=TRUE,padj=TRUE)

mono<-x$res[order(-x$res$log2FoldChange),]
mono_FC <- mono["log2FoldChange"]
mono_FC<-mono_FC$log2FoldChange
names(mono_FC) <- rownames(mono)

library(msigdbr)
all_gene_sets = msigdbr(species = "human", category = "H")
all_gene_setsC2 = msigdbr(species = "human", category = "C2")
unique(all_gene_setsC2$gs_subcat)
all_gene_setsC2<-all_gene_setsC2[grep("CP",all_gene_setsC2$gs_subcat),]
msigdbr_list = split(x = all_gene_sets$gene_symbol, f = all_gene_sets$gs_name)
msigdbr_listC2 = split(x = all_gene_setsC2$gene_symbol, f = all_gene_setsC2$gs_name)
msigdbr_list<-c(msigdbr_list,msigdbr_listC2)

require(readr)
require(fgsea)
senMayo<-read_tsv("~/human/reference/senMayo/SAUL_SEN_MAYO.v2023.2.Hs.tsv")
senMayo_gs<- senMayo[17,2]
senMayo_gs<-str_split_fixed(senMayo_gs,",",125)
abc <-senMayo_gs
names(abc)<-NULL
msigdbr_list_sen<-msigdbr_list
msigdbr_list_sen[["SAUL_SEN_MAYO"]]<-abc

fgseaRes_mes <- fgsea( pathways = msigdbr_list_sen,
                       stats    = mono_FC,
                       
                       minSize  = 15,
                       maxSize  = 500)
topPathwaysUp <- fgseaRes_mes[ES > 0][head(order(pval), n=30), pathway]

topPathwaysDown <- fgseaRes_mes[ES < 0][head(order(pval), n=20), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(msigdbr_list_sen[topPathways], mono_FC, fgseaRes_mes, 
              gseaParam=0.5)


topPathwaysUp_df<-fgseaRes_mes[fgseaRes_mes$pathway%in%topPathwaysUp,]

ggplot(topPathwaysUp_df, aes(x = NES, y = reorder(pathway, NES))) +
  geom_point(aes(size = size, color = -log10(padj))) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway",
    size = "Gene Set Size",
    color = "-log10(FDR)"
  ) +
  theme_minimal()
