sc<-readRDS("PATH/TO/SC_OBJ")
meta.data<-readRDS("PATH/TO/METADATA")

reticulate::use_python("/opt/python/bin/python",required=T)
modules <- reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")
reticulate::py_available()
require(RaceID)
require(Matrix)

A <- intersect(meta.data$cells[grep("HSC",meta.data$cell_type)],meta.data$cells[grep("No",meta.data$cyt)])
B <- intersect(meta.data$cells[grep("HSC",meta.data$cell_type)],meta.data$cells[grep("Yes",meta.data$cyt)])

x <- diffexpnb(getfdata(sc,n=c(A,B)), A=A, B=B )

###Fig8a
plotdiffgenesnb(x,pthr=.05,lthr=.5,mthr=-1,Aname="No_cytopenia",Bname="Yes_cytopenia",show_names=TRUE,padj=TRUE)

HSC<-x$res[order(-x$res$log2FoldChange),]
HSC_FC <- HSC["log2FoldChange"]
HSC_FC<-HSC_FC$log2FoldChange
names(HSC_FC) <- rownames(HSC)

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
senMayo<-read_tsv("./SAUL_SEN_MAYO.v2023.2.Hs.tsv") # Suppl_Table7 of Saul et al., 2022 Pubmed 35974106 
senMayo_gs<- senMayo[17,2]
senMayo_gs<-str_split_fixed(senMayo_gs,",",125)
abc <-senMayo_gs
names(abc)<-NULL
msigdbr_list_sen<-msigdbr_list
msigdbr_list_sen[["SAUL_SEN_MAYO"]]<-abc

fgseaRes_HSC <- fgsea( pathways = msigdbr_list_sen,
                       stats    = HSC_FC,
                       
                       minSize  = 15,
                       maxSize  = 500)
topPathwaysUp <- fgseaRes_HSC[ES > 0][head(order(pval), n=30), pathway]

topPathwaysUp_df<-fgseaRes_HSC[fgseaRes_HSC$pathway%in%topPathwaysUp,]

###Fig8b
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

###Fig8d HSC senescence ~ mes infl

a <- read.csv("PATH/to/suppl_table_7")
u <- log(a$Mes_infl+1)
v <- log(a$HSC_sen+1)
fit <- lm(v ~ u)
summary(fit)
plot(u,v)
abline(coef(fit)[1],coef(fit)[2])
R2 <- cor(u,v,use="pairwise.complete.obs")
legend("topright",paste("R=",round(R2,2),sep=""))
uv <- data_frame(HSC_sen = v, Mes_infl = u)


fit <- lm(HSC_sen ~ Mes_infl,
          data = uv)

pval <- summary(fit)$coefficients[2, 4]
p_label <- paste0("p = ", formatC(pval, format = "f", digits = 5))


ggplot(uv, 
            aes(x = Mes_infl, y = HSC_sen)) +
  geom_point(color = "grey30") +   
  theme_classic() +
  geom_smooth(method = lm, color = "#1f78b4", fill = "#a6cee3") +  
  scale_fill_viridis(discrete = TRUE) +
  sm_statCorr(color = "#1f78b4", fill = "#a6cee3", 
              label_x = -0.2, label_y = 0.1, text_size = 5) +
  sm_classic()

###Fig8d HSC senescence ~ BM in situ CAR freq

u <- log(a$CART_IF_stain)
v <- log(a$HSC_sen+1)
fit <- lm(v ~ u)
summary(fit)
plot(u,v)
abline(coef(fit)[1],coef(fit)[2])
R2 <- cor(u,v,use="pairwise.complete.obs")
legend("topright",paste("R=",round(R2,2),sep=""))
uv <- data_frame(HSC_sen = v, BM_CAR = u)


fit <- lm(HSC_sen ~ BM_CAR,
          data = uv)

pval <- summary(fit)$coefficients[2, 4]
p_label <- paste0("p = ", formatC(pval, format = "f", digits = 5))


ggplot(uv, 
       aes(x = BM_CAR, y = HSC_sen)) +
  geom_point(color = "grey30") +   
  theme_classic() +
  geom_smooth(method = lm, color = "#1f78b4", fill = "#a6cee3") +  
  scale_fill_viridis(discrete = TRUE) +
  sm_statCorr(color = "#1f78b4", fill = "#a6cee3", 
              label_x = -0.2, label_y = 0.1, text_size = 5) +
  sm_classic()


