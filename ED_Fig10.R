###ED10a

library(smplot2)

library(viridis)
library(ggplot2)
library(broom)

a <- read.csv("PATH/to/suppl_table_7")

b <- a[,c("Hb", "PLT",  "ANC")]
bn <- t(t(b)/colSums(b))

gmean <- function(x, na.rm = FALSE) {
  n <- if (na.rm) sum(!is.na(x)) else length(x)
  prod(x, na.rm = na.rm)^(1/n)
}

bng <- apply(bn,1,gmean)
v <- log(bng)

u <- a$age
v <- log(a$HSC_sen + 1)

fit <- lm(v ~ u)
model_summary <- tidy(fit)
p_value <- model_summary$p.value[2]
uv <- data.frame(log_HSC_sen = v, age = u)


ggplot(uv, 
       aes(x = age, y = log_HSC_sen)) +
  geom_point(color = "grey30") +   
  theme_classic() +
  geom_smooth(method = lm, color = "red", fill = "red") +  
  scale_fill_viridis(discrete = TRUE) +
  sm_statCorr(color = "red", fill = "red", 
              text_size = 5) +
  sm_classic()

##ED10c
library(Seurat)
seu <- readRDS("PATH/to/seurat_obj")
DimPlot(seu, reduction = "umap")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes



seu <- CellCycleScoring(seu, 
                        s.features = s.genes, 
                        g2m.features = g2m.genes, 
                        set.ident = TRUE)
library(stringr)
meta.data<-readRDS("PATH/to/suppl_table_8")
seu@meta.data$sample<-meta.data$sample

metadata<- seu@meta.data
metadata_HSC <- metadata[metadata$cell_types%in%"HSC",]

library(ggplot2)

ggplot(metadata_HSC,
       aes(x=sample, fill=Phase))+ theme_classic()+geom_bar(position = "fill") +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,  size=14))

###ED10d
seurat_obj <- readRDS("PATH/to/seurat_obj")
Idents(seurat_obj) <- seurat_obj$cell_types
Idents(seurat_obj) <- "celltype"
unique(Idents(seurat_obj))

seurat_obj = SetIdent(seurat_obj, value = seurat_obj$cell_types)


library(nichenetr)
library(RColorBrewer)
library(tidyverse)
library(Seurat) 
library(dplyr)
organism <- "human"

if(organism == "human"){
  lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
  weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
} else if(organism == "mouse"){
  lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
  weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
  
}

lr_network <- lr_network %>% distinct(from, to)
head(lr_network)

ligand_target_matrix[1:5,1:5] 

head(weighted_networks$lr_sig) 
head(weighted_networks$gr) 

#####step by step

receiver = "HSC"
expressed_genes_receiver <- get_expressed_genes(receiver, seurat_obj, pct = 0.05)
all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

sender_celltypes <- c("Monocyte",
                      
                      "Infl.mesenchymal")

# Use lapply to get the expressed genes of every sender cell type separately here
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_obj, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

###geneset of interest
condition_oi <-  "Yes_Cytopenia"
condition_reference <- "No_Cytopenia"

seurat_obj_receiver <- subset(seurat_obj, idents = receiver)

DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "cytopenia",
                                  min.pct = 0.05) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

require(readr)
senMayo<-read_tsv("./SAUL_SEN_MAYO.v2023.2.Hs.tsv") # Suppl_Table7 of Saul et al., 2022 Pubmed 35974106 
senMayo_gs<- senMayo[17,2]
senMayo_gs<-str_split_fixed(senMayo_gs,",",125)

geneset_oi<-intersect(geneset_oi,senMayo_gs)
geneset_oi<-unique(c(senMayo_gs,senescence_genes))
geneset_oi<-senescence_genes

#####from reactome
library(msigdbr)
library(dplyr)
reactome_sets <- msigdbr(
  species = "Homo sapiens",
  category = "C2",
  subcategory = "CP:REACTOME"
)
sen <- reactome_sets %>%
  filter(gs_name == "REACTOME_CELLULAR_SENESCENCE") %>%
  pull(gene_symbol) %>%
  unique()
######

background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>%
  arrange(desc(aupr_corrected)) %>%
  mutate(rank = rank(-aupr_corrected))


best_upstream_ligands <- ligand_activities %>% top_n(40, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)


vis_ligand_aupr <- ligand_activities %>%
  dplyr::filter(test_ligand %in% best_upstream_ligands) %>%
  tibble::column_to_rownames("test_ligand") %>%
  dplyr::select(aupr_corrected) %>%
  arrange(aupr_corrected) %>%
  as.matrix()

(make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands", "Ligand activity", 
                     legend_title = "AUPR", color = "darkorange") + 
    theme(axis.text.x.top = element_blank()))  

active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

# Dotplot of sender-focused approach
p_dotplot <- DotPlot(subset(seurat_obj, cell_types %in% sender_celltypes),
                     features = rev(best_upstream_ligands), cols = "RdYlBu") + 
  coord_flip() +
  scale_y_discrete(position = "right")

p_dotplot

####LFC


celltype_order <- levels(Idents(seurat_obj)) 

DE_table_top_ligands <- lapply(
  celltype_order[celltype_order %in% sender_celltypes],
  get_lfc_celltype, 
  seurat_obj = seurat_obj,
  condition_colname = "cytopenia",
  condition_oi = condition_oi,
  condition_reference = condition_reference,
  celltype_col = "cell_types",
  min.pct = 0, logfc.threshold = 0,
  features = best_upstream_ligands 
) 

DE_table_top_ligands <- DE_table_top_ligands %>%
  purrr::reduce(dplyr::full_join) %>%
  tibble::column_to_rownames("gene")

vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(best_upstream_ligands), , drop = FALSE])

p_lfc <- make_threecolor_heatmap_ggplot(vis_ligand_lfc,
                                        "Prioritized ligands", "LFC in Sender",
                                        low_color = "midnightblue", mid_color = "white",
                                        mid = median(vis_ligand_lfc), high_color = "red",
                                        legend_title = "LFC")

p_lfc

make_line_plot(
  ligand_activities = as.data.frame(ligand_activities_all),
  potential_ligands = potential_ligands_focused
)

ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 
vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                     y_name = "Ligands", x_name = "Receptors",  
                     color = "mediumvioletred", legend_title = "Prior interaction potential")

