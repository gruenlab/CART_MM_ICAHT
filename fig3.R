reticulate::use_python("/opt/python/bin/python",required=T)
modules <- reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")
reticulate::py_available()
require(RaceID)
require(Matrix)
sc<-readRDS("PATH/TO/SC_OBJ")
meta.data<-readRDS("PATH/TO/suppl_table8")

plotexpmap(sc,"IL6",logsc=T,fr=F, um=T, cex=40,leg = FALSE,cells = meta.data[grep("Yes",meta.data$cyt),]$cells)
plotexpmap(sc,"IL6",logsc=T,fr=F, um=T, cex=40,leg = TRUE,map = FALSE,cells = meta.data[grep("Yes",meta.data$cyt),]$cells)

plotexpmap(sc,"IL6",logsc=T,fr=F, um=T, cex=40,leg = FALSE,cells = meta.data[grep("No",meta.data$cyt),]$cells)
plotexpmap(sc,"IL6",logsc=T,fr=F, um=T, cex=40,leg = TRUE,map = FALSE,cells = meta.data[grep("No",meta.data$cyt),]$cells)

plotexpmap(sc,"CCL2",logsc=T,fr=F, um=T, cex=40,leg = FALSE,cells = meta.data[grep("Yes",meta.data$cyt),]$cells)
plotexpmap(sc,"CCL2",logsc=T,fr=F, um=T, cex=40,leg = TRUE,map = FALSE,cells = meta.data[grep("Yes",meta.data$cyt),]$cells)

plotexpmap(sc,"CCL2",logsc=T,fr=F, um=T, cex=40,leg = FALSE,cells = meta.data[grep("No",meta.data$cyt),]$cells)
plotexpmap(sc,"CCL2",logsc=T,fr=F, um=T, cex=40,leg = TRUE,map = FALSE,cells = meta.data[grep("No",meta.data$cyt),]$cells)

plotexpmap(sc,"CXCL8",logsc=T,fr=F, um=T, cex=40,leg = FALSE,cells = meta.data[grep("Yes",meta.data$cyt),]$cells)
plotexpmap(sc,"CXCL8",logsc=T,fr=F, um=T, cex=40,leg = TRUE,map = FALSE,cells = meta.data[grep("Yes",meta.data$cyt),]$cells)

plotexpmap(sc,"CXCL8",logsc=T,fr=F, um=T, cex=40,leg = FALSE,cells = meta.data[grep("No",meta.data$cyt),]$cells)
plotexpmap(sc,"CXCL8",logsc=T,fr=F, um=T, cex=40,leg = TRUE,map = FALSE,cells = meta.data[grep("No",meta.data$cyt),]$cells)

plotexpmap(sc,"CXCL3",logsc=T,fr=F, um=T, cex=40,leg = FALSE,cells = meta.data[grep("Yes",meta.data$cyt),]$cells)
plotexpmap(sc,"CXCL3",logsc=T,fr=F, um=T, cex=40,leg = TRUE,map = FALSE,cells = meta.data[grep("Yes",meta.data$cyt),]$cells)

plotexpmap(sc,"CXCL3",logsc=T,fr=F, um=T, cex=40,leg = FALSE,cells = meta.data[grep("No",meta.data$cyt),]$cells)
plotexpmap(sc,"CXCL3",logsc=T,fr=F, um=T, cex=40,leg = TRUE,map = FALSE,cells = meta.data[grep("No",meta.data$cyt),]$cells)

plotexpmap(sc,"CXCL5",logsc=T,fr=F, um=T, cex=40,leg = FALSE,cells = meta.data[grep("Yes",meta.data$cyt),]$cells)
plotexpmap(sc,"CXCL5",logsc=T,fr=F, um=T, cex=40,leg = TRUE,map = FALSE,cells = meta.data[grep("Yes",meta.data$cyt),]$cells)

plotexpmap(sc,"CXCL5",logsc=T,fr=F, um=T, cex=40,leg = FALSE,cells = meta.data[grep("No",meta.data$cyt),]$cells)
plotexpmap(sc,"CXCL5",logsc=T,fr=F, um=T, cex=40,leg = TRUE,map = FALSE,cells = meta.data[grep("No",meta.data$cyt),]$cells)

plotexpmap(sc,"CCL26",logsc=T,fr=F, um=T, cex=40,leg = FALSE,cells = meta.data[grep("Yes",meta.data$cyt),]$cells)
plotexpmap(sc,"CCL26",logsc=T,fr=F, um=T, cex=40,leg = TRUE,map = FALSE,cells = meta.data[grep("Yes",meta.data$cyt),]$cells)

plotexpmap(sc,"CCL26",logsc=T,fr=F, um=T, cex=40,leg = FALSE,cells = meta.data[grep("No",meta.data$cyt),]$cells)
plotexpmap(sc,"CCL26",logsc=T,fr=F, um=T, cex=40,leg = TRUE,map = FALSE,cells = meta.data[grep("No",meta.data$cyt),]$cells)

plotexpmap(sc,"PTGS2",logsc=T,fr=F, um=T, cex=40,leg = FALSE,cells = meta.data[grep("Yes",meta.data$cyt),]$cells)
plotexpmap(sc,"PTGS2",logsc=T,fr=F, um=T, cex=40,leg = TRUE,map = FALSE,cells = meta.data[grep("Yes",meta.data$cyt),]$cells)

plotexpmap(sc,"PTGS2",logsc=T,fr=F, um=T, cex=40,leg = FALSE,cells = meta.data[grep("No",meta.data$cyt),]$cells)
plotexpmap(sc,"PTGS2",logsc=T,fr=F, um=T, cex=40,leg = TRUE,map = FALSE,cells = meta.data[grep("No",meta.data$cyt),]$cells)

###Fig3b DEG + GSEA
A <- intersect(meta.data$cells[grep("esenchy",meta.data$cell_type)],meta.data$cells[grep("No",meta.data$cyt)])
B <- intersect(meta.data$cells[grep("esenchy",meta.data$cell_type)],meta.data$cells[grep("Yes",meta.data$cyt)])

x <- diffexpnb(getfdata(sc,n=c(A,B)), A=A, B=B )
plotdiffgenesnb(x,pthr=.05,lthr=.5,mthr=-1,Aname="No_cytopenia",Bname="Yes_cytopenia",show_names=TRUE,padj=TRUE)

Mesen<-x$res[order(-x$res$log2FoldChange),]
Mesen_FC <- Mesen["log2FoldChange"]
Mesen_FC<-Mesen_FC$log2FoldChange
names(Mesen_FC) <- rownames(Mesen)

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
                       stats    = Mesen_FC,
                       
                       minSize  = 15,
                       maxSize  = 500)
topPathwaysUp <- fgseaRes_mes[ES > 0][head(order(pval), n=20), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(msigdbr_list_sen[topPathways], Mesen_FC, fgseaRes_mes, 
              gseaParam=0.5)

topPathways <- fgseaRes_mes %>%
  arrange(-NES) %>%
  head(20)
ggplot(topPathways, aes(x = NES, y = reorder(pathway, NES))) +
  geom_point(aes(size = size, color = -log10(padj))) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway",
    size = "Gene Set Size",
    color = "-log10(FDR)"
  ) +
  theme_minimal()

###Fig3c correlation of clinical scores and mes infl scores
z <- read.csv("PATH/to/suppl_table_6")
y <- z[,c("Hb", "PLT",  "ANC")]
yn <- t(t(y)/colSums(y))

gmean <- function(x, na.rm = FALSE) {
  n <- if (na.rm) sum(!is.na(x)) else length(x)
  prod(x, na.rm = na.rm)^(1/n)
}

yng <- apply(yn,1,gmean)

i <- 6
x <- z[,i] + .1
f <- x > -Inf
u <- log(x[f])
v <- log(yng[f])
fit <- lm(v ~ u)
summary(fit)
plot(u,v)
abline(coef(fit)[1],coef(fit)[2])
R2 <- cor(u,v,use="pairwise.complete.obs")
legend("topright",paste("R=",round(R2,2),sep=""))

library(smplot2)
library(ggpmisc)
library(broom)
library(tibble)
model_summary <- tidy(fit)
p_value <- model_summary$p.value[2]
uv <- data_frame(clinical_score = v, mes_infl = u)

ggplot(uv, aes(x = mes_infl, y = clinical_score)) +
  geom_point() +
  theme_classic()+
  geom_smooth(method=lm , color="red", fill="#69b3a2") +
  sm_statCorr(color="red", fill="#69b3a2",label_x = -4,
              label_y = -2.2,
              text_size = 5) +sm_classic()


###Fig 3d correlation of blood IL6 peak levels and mes infl
a <- read.csv("PATH/to/suppl_table_7")

u <- log(a$IL6_peak)
v <- a$Mes_infl

fit <- lm(v ~ u)
model_summary <- tidy(fit)
p_value <- model_summary$p.value[2]
uv <- data_frame(Mes_infl = v, IL6_peak = u)

library(viridis)
ggplot(uv, 
       aes(x = IL6_peak, y = Mes_infl)) +
  geom_point(color = "grey30") +   
  theme_classic() +
  geom_smooth(method = lm, color = "#1f78b4", fill = "#a6cee3") +  
  scale_fill_viridis(discrete = TRUE) +
  sm_statCorr(color = "#1f78b4", fill = "#a6cee3", 
              text_size = 5) +
  sm_classic()

###Fig3d correlation of blood IL6 peak levels and clinical score

b <- a[,c("Hb", "PLT",  "ANC")]
bn <- t(t(b)/colSums(b))

gmean <- function(x, na.rm = FALSE) {
  n <- if (na.rm) sum(!is.na(x)) else length(x)
  prod(x, na.rm = na.rm)^(1/n)
}

bng <- apply(bn,1,gmean)
v <- log(bng)
u <- log(a$IL6_peak)

fit <- lm(v ~ u)
model_summary <- tidy(fit)
p_value <- model_summary$p.value[2]
uv <- data_frame(Mes_infl = v, IL6_peak = u)

library(viridis)
ggplot(uv, 
       aes(x = IL6_peak, y = Mes_infl)) +
  geom_point(color = "grey30") +   
  theme_classic() +
  geom_smooth(method = lm, color = "#1f78b4", fill = "#a6cee3") +  
  scale_fill_viridis(discrete = TRUE) +
  sm_statCorr(color = "#1f78b4", fill = "#a6cee3", 
              text_size = 5) +
  sm_classic()

### Fig3e mes infl and CRS grade
df <- data_frame(clinical_score = a$Mes_infl, treatment_intensity = a$CRS_grade)
library(ggplot2)
library(dplyr)

df <- df %>%
  mutate(treatment_intensity = as.factor(treatment_intensity))

kw <- kruskal.test(clinical_score ~ treatment_intensity, data = df)

p_label <- paste0("p = ", signif(kw$p.value, 3))

ggplot(df, aes(x = treatment_intensity, y = clinical_score)) +
  
  
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = "white", color = "black") +
  
  
  geom_jitter(width = 0.12, size = 2.5, alpha = 0.8) +
  
  
  stat_summary(fun = median, geom = "point", size = 3, color = "red") +
  

  labs(
    x = "treatment_intensity",
    y = "IL6_peak",
    title = NULL,
    subtitle = p_label
  ) +
  
 
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

###Fig3f correlation of mes infl and baseline clinical score
b <- z[,c("Hb_baseline", "PLT_baseline",  "ANC_baseline")]

bn <- t(t(b)/colSums(b))

gmean <- function(x, na.rm = FALSE) {
  n <- if (na.rm) sum(!is.na(x)) else length(x)
  prod(x, na.rm = na.rm)^(1/n)
}

bng <- apply(bn,1,gmean)
v <- log(bng)
u <- log(z$mes_infl_score + .1)

fit <- lm(v ~ u)
model_summary <- tidy(fit)
p_value <- model_summary$p.value[2]
uv <- data_frame(baseline_clinical_score = v, Mes_infl = u)

ggplot(uv, aes(x = Mes_infl, y = baseline_clinical_score)) +
  geom_point() +
  theme_classic()+
  geom_smooth(method=lm , color="red", fill="#69b3a2") +
  sm_statCorr(color="red", fill="#69b3a2", text_size = 5) +sm_classic()

###Fig3g correlation of mes infl and baseline blood IL6level

u <- log(a$IL6_pre)
v <- a$Mes_infl

fit <- lm(v ~ u)
model_summary <- tidy(fit)
p_value <- model_summary$p.value[2]
uv <- data_frame(Mes_infl = v, IL6_peak = u)

library(viridis)
ggplot(uv, 
       aes(x = IL6_peak, y = Mes_infl)) +
  geom_point(color = "grey30") +   
  theme_classic() +
  geom_smooth(method = lm, color = "#1f78b4", fill = "#a6cee3") +  
  scale_fill_viridis(discrete = TRUE) +
  sm_statCorr(color = "#1f78b4", fill = "#a6cee3", 
              text_size = 5) +
  sm_classic()




