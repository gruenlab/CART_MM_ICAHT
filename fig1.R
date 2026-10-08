reticulate::use_python("/opt/python/bin/python",required=T)
modules <- reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")
reticulate::py_available()
require(RaceID)
require(Matrix)
sc<-readRDS("PATH/TO/SC_OBJ")
meta.data<-readRDS("PATH/TO/suppl_table8")

plotsymbolsmap2 <- function (object, types, subset = NULL, samples_col = NULL, cex = 0.5, 
                             fr = FALSE, um = FALSE, leg = TRUE, map = TRUE, cex.legend = 0.75, 
                             leg.pos = "topleft", rand = FALSE, seed=123) 
{
  if (length(object@tsne) == 0 & length(object@fr) == 0 & length(object@umap) == 
      0) 
    stop("run comptsne/compfr/compumap before plotlabelsmap")
  if (!is.logical(fr)) 
    stop("fr has to be TRUE or FALSE")
  if (!is.logical(um)) 
    stop("um has to be TRUE or FALSE")
  if (fr == FALSE & um == FALSE & dim(object@tsne)[1] == 0) {
    if (dim(object@fr)[1] != 0) {
      fr <- TRUE
    }
    else if (dim(object@umap)[1] != 0) {
      um <- TRUE
    }
  }
  if (is.null(subset)) 
    subset <- unique(types)
  h <- sort(unique(types)) %in% subset
  if (!is.null(subset)) {
    fp <- rep(FALSE, length(types))
    fp[types %in% subset] <- TRUE
  }
  if (is.null(samples_col)) {
    samples_col <- rainbow(length(unique(types[fp])))
  }
  else {
    samples_col <- samples_col[h]
  }
  if (fr) {
    d <- object@fr
  }
  else if (um) {
    d <- object@umap
  }
  else {
    d <- object@tsne
  }
  if (map) {
    if (rand){
      plot(d, xlab = "", ylab = "", axes = FALSE, cex = cex, 
           pch = 20, col = "grey")
      set.seed(seed)
      indx <- sample(1:length(types))
      tp   <- sort(unique(types[fp]))
      for ( i in indx ){
        j <- which( tp == types[i] )
        points(d[i, 1], d[i, 2], col = samples_col[j], pch = 20, cex = cex)
      }
    }else{
      plot(d, xlab = "", ylab = "", axes = FALSE, cex = cex, 
           pch = 20, col = "grey")
      for (i in 1:length(unique(types[fp]))) {
        f <- types == sort(unique(types[fp]))[i]
        points(d[f, 1], d[f, 2], col = samples_col[i], pch = 20, 
               cex = cex)
      }
    }   
  }
  else {
    plot(d, xlab = "", ylab = "", axes = FALSE, cex = 0, 
         pch = 20, col = "grey", xlim = c(min(d[, 1]), max(d[, 
                                                             1])), ylim = c(min(d[, 2]), max(d[, 2])))
  }
  if (leg) 
    legend(leg.pos, legend = sort(unique(types[fp])), col = samples_col, 
           pch = 20, cex = cex.legend, bty = "n")
}

###Fig1c
plotsymbolsmap(sc,types = meta.data$cell_type,um = T,leg = FALSE, cex=5)
plotsymbolsmap(sc,types = meta.data$cell_type,um = T,map = FALSE)
###Fig1e
plotsymbolsmap(sc,types = meta.data$sample,um = T,leg = FALSE, cex=5)
plotsymbolsmap(sc,types = meta.data$sample,um = T,map = FALSE)
###Fig1f
col2<-c("red","blue","yellow")
plotsymbolsmap2(sc,types = meta.data$cyt,um = T,leg = FALSE, cex=5, samples_col=col2, rand = TRUE, seed = 12345)
plotsymbolsmap(sc,types = meta.data$cyt,um = T,leg = TRUE, map = FALSE, cex=0.2, samples_col=col2)

###Fig1d
cluster_order <- c("HSC","B.prog",  "B.cell","T.cell","GMP",              
                   "NP", "FCGR3B+.neutrophil","Defensin+.neutrophil", "LTF+.Neutrophil" ,
                   "MMP9+.neutrophil","MDP" ,"Monocyte","Macrophage","pDC","cDC",
                   "MEP","Megakaryocyte","Erythrocyte","Mast.cell",
                   "Mesenchymal" , "Infl.mesenchymal", "Myeloma.mesenchymal" ,
                   "SMC",
                   "Osteolineage", 
                   "Endothelial")

cluster_order<-rev(cluster_order)
genes=c("HLF","SPINK2","AVP","CD34","CD38","VPREB1","CD19","CD3E","MPO","ELANE","ITGAM","FCGR3B","DEFA4",
        "LTF","MMP9","CD14","C1QB","CLEC4C","LILRA4","CD1C","CLEC10A","PPBP","PF4","GYPA","HDC",
        "CPA3","PDGFRA","IL6","CCL2","PTX3","FST","ACTA2","IBSP","IFITM5","CDH5")


fractDotPlot(sc, genes, cluster=cluster_order,logscale=TRUE)
###Fig1g
cluster_order <- c("Mesenchymal" , "Infl.mesenchymal", "Myeloma.mesenchymal"
)
cluster_order<-rev(cluster_order)
genes=c("PDGFRA","CXCL12","LEPR","IL6","CCL2","CXCL3","CXCL8","PTX3","FST")
fractDotPlot(sc, genes, cluster=cluster_order,logscale=T)

###Fig1h
library(tidyverse)
library(dplyr)

####prog
meta.data_prog<-meta.data[grep("P$|HSC|B.prog",meta.data$broad_ct),]
counts <- meta.data_prog %>%
  group_by(sample, broad_ct) %>%
  summarise(Count = n(), .groups = 'drop')
total_counts <- counts %>%
  group_by(sample) %>%
  summarise(Total = sum(Count), .groups = 'drop')
proportions <- counts %>%
  left_join(total_counts, by = "sample") %>%
  mutate(Proportion = Count / Total) %>%
  dplyr::select(sample, broad_ct, Proportion)

proportions$group<-proportions$sample
proportions$group[grep("^NC",proportions$group)]<-"No_Cytopenia"
proportions$group[grep("^C",proportions$group)]<-"Cytopenia"
proportions<-proportions %>%
  filter(group != "P1")
#proportions<-proportions %>%
# filter(sample != "C1")

# Perform Mann-Whitney U test for each cell type
results <- proportions %>%
  group_by(broad_ct) %>%
  summarise(
    p_value_great = t.test(Proportion ~ group, alternative = "greater")$p.value,
    p_value_less = t.test(Proportion ~ group, alternative = "less")$p.value,
    p_value_two = t.test(Proportion ~ group)$p.value
    
  )
print(results)


ggplot(proportions, aes(x = broad_ct, y = Proportion, fill = group)) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.7, outlier.shape = NA) + 
  geom_jitter(aes(color = group), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 2) +  
  theme_minimal() +
  labs(title = "Cell Type Proportions by Sample", x = "Cell Type", y = "Proportion") +
  theme(legend.position = "top")

###C45+ CD34-
meta.data_mye<-meta.data[grep("DC|Macroph|Mast|Mega|Mono|Neu|T.cell|B.cell",meta.data$broad_ct),]
counts <- meta.data_mye %>%
  group_by(sample, broad_ct) %>%
  summarise(Count = n(), .groups = 'drop')
total_counts <- counts %>%
  group_by(sample) %>%
  summarise(Total = sum(Count), .groups = 'drop')
proportions <- counts %>%
  left_join(total_counts, by = "sample") %>%
  mutate(Proportion = Count / Total) %>%
  dplyr::select(sample, broad_ct, Proportion)

proportions$group<-proportions$sample
proportions$group[grep("^NC",proportions$group)]<-"No_Cytopenia"
proportions$group[grep("^C",proportions$group)]<-"Cytopenia"
proportions<-proportions %>%
  filter(group != "P1")

# Perform Mann-Whitney U test for each cell type
results <- proportions %>%
  group_by(broad_ct) %>%
  summarise(
  
    p_value_two = t.test(Proportion ~ group)$p.value
  )


results_t_imm <- proportions %>%
  group_by(broad_ct) %>%
  summarise(
    p_value_two = t.test(Proportion ~ group)$p.value
  )

results_wil_imm <- proportions %>%
  group_by(broad_ct) %>%
  summarise(
    p_value_great = wilcox.test(Proportion ~ group, alternative = "greater")$p.value,
    p_value_less = wilcox.test(Proportion ~ group, alternative = "less")$p.value,
    p_value_two = wilcox.test(Proportion ~ group)$p.value
  )

print(results)
ggplot(proportions, aes(x = broad_ct, y = Proportion, fill = group)) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.7, outlier.shape = NA) +  # Side by side boxplots
  geom_jitter(aes(color = group), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 2) +  # Jittered points
  theme_minimal() +
  labs(title = "Cell Type Proportions by Sample", x = "Cell Type", y = "Proportion") +
  theme(legend.position = "top")

meta.data_prog<-meta.data[grep("Niche",meta.data$broad_ct),]
counts <- meta.data_prog %>%
  group_by(sample, cell_type) %>%
  summarise(Count = n(), .groups = 'drop')
total_counts <- counts %>%
  group_by(sample) %>%
  summarise(Total = sum(Count), .groups = 'drop')
proportions <- counts %>%
  left_join(total_counts, by = "sample") %>%
  mutate(Proportion = Count / Total) %>%
  dplyr::select(sample, cell_type, Proportion)

proportions$group<-proportions$sample
proportions$group[grep("^NC",proportions$group)]<-"No_Cytopenia"
proportions$group[grep("^C",proportions$group)]<-"Cytopenia"
proportions<-proportions %>%
  filter(group != "P1")

# Perform Mann-Whitney U test for each cell type
results <- proportions %>%
  group_by(broad_ct) %>%
  summarise(
    p_value_great = t.test(Proportion ~ group, alternative = "greater")$p.value,
    p_value_less = t.test(Proportion ~ group, alternative = "less")$p.value,
    p_value_two = t.test(Proportion ~ group)$p.value
    
  )
print(results)

results_w <- proportions %>%
  group_by(cell_type) %>%
  summarise(
    p_value_great = wilcox.test(Proportion ~ group, alternative = "greater")$p.value,
    p_value_less = wilcox.test(Proportion ~ group, alternative = "less")$p.value,
    p_value_two = wilcox.test(Proportion ~ group)$p.value
  )

ggplot(proportions, aes(x = cell_type, y = Proportion, fill = group)) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.7, outlier.shape = NA) +  # Side by side boxplots
  geom_jitter(aes(color = group), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), size = 2) +  # Jittered points
  theme_minimal() +
  labs(title = "Cell Type Proportions by Sample", x = "Cell Type", y = "Proportion") +
  theme(legend.position = "top")


