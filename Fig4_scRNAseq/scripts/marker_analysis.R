#update.packages(ask=FALSE)
set.seed(0)  # random seed for reproducibility
rm(list = ls())
gc()
load_pkg <- rlang::quos(tidyverse, RColorBrewer, pheatmap, ggplot2, dplyr, stringr, heatmaply, Seurat,
                        svglite, DataCombine, biomaRt, GOSemSim, org.Hs.eg.db, clusterProfiler, DOSE,
                        enrichplot, ggplot2, ggrepel, scales)
invisible(lapply(lapply(load_pkg, rlang::quo_name), library, character.only = TRUE))
`%nin%` = Negate(`%in%`)

# Load Seurat_Integrated Datasets and New Meta with NGN1
seurat_integrated <- readRDS("../data/JH_org_sc_seurat_integrated_CCA_2batch_clustered.rds")
new.meta <- readRDS("../meta/JH_org_sc_seurat_integrated_CCA_2batch_clustered_meta_with_NGN1.rds")

# Load Marker Lists Endo and Neuro
marker_df <- read.csv("../meta/marker_list.csv", fileEncoding="UTF-8-BOM")
neuro_marker_df <- read.csv("../meta/marker_list_neuro.csv", fileEncoding="UTF-8-BOM")

# Known Markers 
marker_list <- marker_df$Gene
neuro_marker_list <- neuro_marker_df$Gene

# Edit new.meta to include Endothelial in 'is_NGN1'
new.meta$is_NGN1 <- ifelse(new.meta$integrated_snn_res.0.1==7, "Endothelium", new.meta$is_NGN1)
new.meta$sample_clust <- paste(new.meta$sample_type, new.meta$integrated_snn_res.0.1, sep="_")
new.meta$sample_NGN1 <- paste(new.meta$sample_type, new.meta$is_NGN1, sep="_")
seurat_integrated@meta.data <- new.meta

# Get All Genes and find out Entrez and Ensembl for GO 
all_genes <- seurat_integrated@assays$RNA@counts@Dimnames[[1]]
all_genes_df <- seurat_integrated@assays$RNA@counts@Dimnames[[1]] %>% as.data.frame()

# Find ensembl and entrez IDs
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene_IDs <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
                  filters = "external_gene_name", values = all_genes, mart = mart)

# Create output directory
output_dir = paste("../results/data_analysis_scRNAseq_", format(Sys.Date(), "%Y%m%d"), sep="")
dir.create(output_dir, recursive = TRUE)

# Create Plot PDF Fn and Volcano
pdfFN <- function(path, object, w, h){
  pdf(path, width = w, height = h)
  print(object)
  dev.off()
}

#################################### Section 1.5: Counting Endothelium: Fig4|S14 ##################################
# Return counts in integrated_snn_res.0.1 for 12 samples which were ran. 
DefaultAssay(seurat_integrated) <- "integrated" # Unsure
orig.ident_list <- unique(new.meta$orig.ident)
numcells_df <- data.frame(matrix(ncol=12, nrow=8))
colnames(numcells_df) <- orig.ident_list
rownames(numcells_df) <- 0:7

for(x in orig.ident_list){
  subset = new.meta[new.meta$orig.ident == x,]
  for (y in 0:7){
    temp <- table(subset$integrated_snn_res.0.1)
    numcells_df[as.character(y),x] <- temp[names(temp)==y]}
}

write.csv(numcells_df, paste(output_dir, "/Fig4_S14_counts_snn_res0.1.csv", sep=""))   #
################################## Section 1.6 is_NGN1 Counts: Fig4|S14 #############################################
#is_NGN1 is thresholded above log2x=4 (NGN1 barcode count > 16)
all_ngn.meta <- new.meta[which(new.meta$NGN1_count>0),]
ngncells_df <- data.frame(matrix(ncol=12, nrow=8))
colnames(ngncells_df) <- orig.ident_list
rownames(ngncells_df) <- 0:7
for(x in orig.ident_list){
  subset = all_ngn.meta[all_ngn.meta$orig.ident == x,]
  for (y in 0:7){
    temp <- table(subset$integrated_snn_res.0.1)
    ngncells_df[as.character(y),x] <- temp[names(temp)==y]}
}

write.csv(ngncells_df, paste(output_dir,"/Fig4_allngn1_counts.csv", sep=""))

########################################### Section II: Heatmaps ###############################################
# Make Diagonal
draw_colnames_45 <- function(coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...))
  return(res)}
assignInNamespace(x = "draw_colnames", value = "draw_colnames_45", ns = asNamespace("pheatmap"))

#################### Heatmaps for Endo using Average Normalized Expression FigS15 #############################

DefaultAssay(seurat_integrated) <- "RNA"
mural_list <- c(7, 0, 1)
Idents(object = seurat_integrated) <- "integrated_snn_res.0.1"
df_ane <- AverageExpression(object = seurat_integrated, assays = "RNA")$RNA
df_EndoMarkers <- df_ane[(rownames(df_ane) %in% marker_list),] %>% as.data.frame()
df_EndoMarkers <- df_EndoMarkers[,(colnames(df_EndoMarkers) %in% mural_list)]
df_EndoMarkers <- df_EndoMarkers[match(marker_list, rownames(df_EndoMarkers)),]
color = colorRampPalette(c("blue", "white", "red"))(499)
df_EndoMarkers <- as.data.frame(t(as.matrix(df_EndoMarkers)))
hmp = pheatmap(df_EndoMarkers, color = color, cluster_rows = F, cluster_cols = F, show_rownames = T, 
  scale = "column", border_color = NA, fontsize = 7, fontsize_row = 7, height = 20)

ggsave(file=paste(output_dir,"/FigS15_Endo_Norm_Counts.svg", sep=""), plot=hmp, width=3.75, height=1.5) #  Endothelium heatmap
write.csv(df_EndoMarkers, paste(output_dir,"/FigS15_Endo_Norm_Counts.csv", sep=""))

##################### Heatmaps for NGN1 using Average Normalized Expression| FigS15 ###########################
neural_list <- c(1, 2, 5, 6)
Idents(object = seurat_integrated) <- "is_NGN1"
df_neuro <- AverageExpression(object = seurat_integrated, assays = "RNA")$RNA
df_NGN1 <- df_neuro[(rownames(df_neuro) %in% neuro_marker_list),] %>% as.data.frame()
df_NGN1 <- df_NGN1[match(neuro_marker_list, rownames(df_NGN1)),]

# Get Average Others
df_NeuroMarkers <- df_ane[(rownames(df_ane) %in% neuro_marker_list),] %>% as.data.frame()
df_NeuroMarkers <- df_NeuroMarkers[,(colnames(df_NeuroMarkers) %in% neural_list)]
df_NeuroMarkers <- df_NeuroMarkers[match(neuro_marker_list, rownames(df_NeuroMarkers)),]
df_NeuroMarkers$NGN1 <- df_NGN1$NGN1
df_NeuroMarkers <- df_NeuroMarkers[, c(1, 5, 2, 6, 'NGN1')] %>%  rename('iNeuron' = 'NGN1')
df_NeuroMarkers <- as.data.frame(t(as.matrix(df_NeuroMarkers)))
color = colorRampPalette(c("blue", "white", "red"))(499)
nmp = pheatmap(df_NeuroMarkers, color = color, cluster_rows = F, cluster_cols = F, show_rownames = T, 
  scale = "column", border_color = NA, fontsize = 7, fontsize_row = 7, height = 20)

ggsave(file=paste(output_dir,"/FigS15_Neuron_Norm_Counts.svg", sep=""), plot=nmp, width=3.75, height=1.5) #  iNeuron heatmap
write.csv(df_NeuroMarkers, paste(output_dir,"/FigS15_Neuron_Norm_Counts.csv", sep=""))

################################ Section III: GO Analysis and Volcano Plots comparing ##########################

### Mixed_Css_Markers DE Analysis 
Idents(object = seurat_integrated) <- "sample_clust"
DefaultAssay(seurat_integrated) <- "RNA"
CSS_MIX = list()
CSS_WT = list()
MIX_WT = list()

for (x in 0:7){
  ml_temp <- FindMarkers(seurat_integrated, ident.1 = paste("CSS", x, sep="_"), ident.2 = paste("Mix", x, sep="_"), 
    logfc.threshold=0.5) %>% arrange(desc(avg_log2FC))
  ml_temp$gene <- rownames(ml_temp)
  ml_temp <- left_join(ml_temp, gene_IDs, by = c("gene"="external_gene_name"))
  CSS_MIX[[as.character(x)]] <- ml_temp[!duplicated(ml_temp$gene),]
    if (x!=7) {
    ml_temp <- FindMarkers(seurat_integrated, ident.1 = paste("CSS", x, sep="_"), ident.2 = paste("WT", x, sep="_"), 
      logfc.threshold=0.5) %>% arrange(desc(avg_log2FC))
    ml_temp$gene <- rownames(ml_temp)
    ml_temp <- left_join(ml_temp, gene_IDs, by = c("gene"="external_gene_name"))
    CSS_WT[[as.character(x)]] <- ml_temp[!duplicated(ml_temp$gene),]
    ml_temp <- FindMarkers(seurat_integrated, ident.1 = paste("Mix", x, sep="_"), ident.2 = paste("WT", x, sep="_"), 
      logfc.threshold=0.5) %>% arrange(desc(avg_log2FC))
    ml_temp$gene <- rownames(ml_temp)
    ml_temp <- left_join(ml_temp, gene_IDs, by = c("gene"="external_gene_name"))
    MIX_WT[[as.character(x)]] <- ml_temp[!duplicated(ml_temp$gene),]
  }
}

# CREATE ENSEMBL GO Analysis Function
ont_list = c("CC", "MF", "BP")
pair_list = c("CSS_MIX", "CSS_WT", "MIX_WT")
range_vector = c(7, 6, 6)

GOFn <- function(input, pair_name, z){
    GO_Ensembl_list = list()
    dir.create(paste(output_dir,"/ENSEMBL_GO/", pair_name, "/CC" , sep=""), recursive=TRUE)
    dir.create(paste(output_dir,"/ENSEMBL_GO/", pair_name, "/MF" , sep=""))
    dir.create(paste(output_dir,"/ENSEMBL_GO/", pair_name, "/BP" , sep=""))
    for (x in 0:z){
    allOE_genes <- gene_IDs$ensembl_gene_id
    sigOE_genes <- input[[as.character(x)]]$ensembl_gene_id
    if (!is_empty(sigOE_genes)) {
      for (y in ont_list){
        ego <- enrichGO(gene = sigOE_genes, universe = allOE_genes, keyType = "ENSEMBL", OrgDb = org.Hs.eg.db, ont = y, 
          pAdjustMethod = "BH", readable = TRUE, qvalueCutoff = 0.05, pvalueCutoff = 0.05)
        d <- godata('org.Hs.eg.db', ont=y)
        if (!is_empty(ego)) {
          if (!nrow(ego)==0){
            ego <- pairwise_termsim(ego, method = "Wang", semData = d)
          }
        }
        print(paste(output_dir, "/ENSEMBL_GO/", pair_name, "/", y, "/", x, "_ego.rda", sep=""))
        save(ego, file=paste(output_dir, "/ENSEMBL_GO/", pair_name, "/", y, "/", x, "_ego.rda", sep=""))
        write.csv(ego, paste(output_dir, "/ENSEMBL_GO/", pair_name, "/", y, "/", x,"_ego.csv", sep=""))
        ID = paste(x, y, sep="_")
        GO_Ensembl_list[[ID]] <- ego}}}
    return(GO_Ensembl_list)
}

# Create Volcano Plotting Function
pdfVol <- function(folder, object, x){
  padj.cutoff <- 0.05
  lfc.cutoff <- 0.5
  dir.create(folder, recursive = TRUE)

  res_tableOE_tb <- object %>% arrange(desc(avg_log2FC))
  res_tableOE_tb <- res_tableOE_tb %>% 
    mutate(threshold_OE = p_val_adj < padj.cutoff & abs(avg_log2FC) >= lfc.cutoff)
  res_tableOE_tb <- res_tableOE_tb %>% arrange(desc(abs(avg_log2FC)))
  res_tableOE_tb <- mutate(res_tableOE_tb, genelabels = ifelse(threshold_OE == TRUE, gene, ""))
  write.csv(res_tableOE_tb,paste(folder,x,"_volcano.csv", sep=""))
  pdfFN(paste(folder,x,"_volcano.pdf", sep=""),
             ggplot(res_tableOE_tb, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
               geom_point(aes(colour = threshold_OE)) +
               geom_text_repel(aes(label = genelabels), max.overlaps = Inf) +
               ggtitle(paste(x, "Volcano Plot", sep=" ")) +
               xlab("log2 fold change") + 
               ylab("-log10 adjusted p-value") +
               geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "blue") +
               geom_vline(xintercept=-0.5, linetype="dashed", color = "blue") +
               geom_vline(xintercept=0.5, linetype="dashed", color = "blue") +
               theme(legend.position = "none",
                     plot.title = element_text(size = rel(1.5), hjust = 0.5),
                     axis.title = element_text(size = rel(1.25))),
                     6, 8 )
}

# Run GO, Save output.
go_Ensembl_CSS_MIX <- GOFn(CSS_MIX, pair_list[1], range_vector[1])
go_Ensembl_CSS_WT <- GOFn(CSS_WT, pair_list[2], range_vector[2])
go_Ensembl_MIX_WT <- GOFn(MIX_WT, pair_list[3], range_vector[3])
listGO <- list(go_Ensembl_CSS_MIX, go_Ensembl_CSS_WT, go_Ensembl_MIX_WT)
listOE <- list (CSS_MIX, CSS_WT, MIX_WT)

## GO Print
for (t in 1:3){
  pair_name <- pair_list[t]
  GO_temp <- listGO[[t]]
  GO_OE<- listOE[[t]]
  for (x in 0:range_vector[t]){
    for (y in ont_list){
      go_Ensembl <- GO_temp[[paste(x,y,sep="_")]]
      pair_OE <- GO_temp[[as.character(x)]]
      OE_foldchanges <- pair_OE$avg_log2FC %>% as.vector()
      names(OE_foldchanges) <-pair_OE$gene
      print(paste(x, y))
      if (!is_empty(go_Ensembl)){
        if (nrow(go_Ensembl)>=2){
          pdfFN(paste(output_dir, "/ENSEMBL_GO/", pair_name, "/",y,"/",x,"_EnsemblGO.pdf", sep=""), 
            dotplot(go_Ensembl, showCategory=20, orderBy = "x") + ggtitle(paste("Cluster ", x," GO Term", y, " - ", pair_name)),
            12, 8)
        }}}}
}

# Generate Volcano Plots
Idents(object = seurat_integrated) <- "sample_clust"
DefaultAssay(seurat_integrated) <- "RNA"
Mix_ccs2 = list()
WT_ccs2 = list()
WT_Mix2 = list()

for (x in 0:7){
  ml_temp <- FindMarkers(seurat_integrated, ident.1 = paste("CSS", x, sep="_"), ident.2 = paste("Mix", x, sep="_"), logfc.threshold=0.1) %>% arrange(desc(avg_log2FC))
  ml_temp$gene <- rownames(ml_temp)
  ml_temp <- left_join(ml_temp, gene_IDs, by = c("gene"="external_gene_name"))
  Mix_ccs2 [[as.character(x)]] <- ml_temp[!duplicated(ml_temp$gene),]
  #write.table(Mix_ccs2 [[as.character(x)]], file = paste("./output/", x, "_mix_css_marker2", ".txt", sep=""), sep = "\t", row.names = T, col.names = T)
  
  ml_temp <- FindMarkers(seurat_integrated, ident.1 = paste("CSS", x, sep="_"), ident.2 = paste("WT", x, sep="_"), logfc.threshold=0.1) %>% arrange(desc(avg_log2FC))
  ml_temp$gene <- rownames(ml_temp)
  ml_temp <- left_join(ml_temp, gene_IDs, by = c("gene"="external_gene_name"))
  WT_ccs2 [[as.character(x)]] <- ml_temp[!duplicated(ml_temp$gene),]
  
  ml_temp <- FindMarkers(seurat_integrated, ident.1 = paste("Mix", x, sep="_"), ident.2 = paste("WT", x, sep="_"), logfc.threshold=0.1) %>% arrange(desc(avg_log2FC))
  ml_temp$gene <- rownames(ml_temp)
  ml_temp <- left_join(ml_temp, gene_IDs, by = c("gene"="external_gene_name"))
  WT_Mix2 [[as.character(x)]] <- ml_temp[!duplicated(ml_temp$gene),]
}

# For Mix_Css2 Print 7 as well
for (x in 0:7){
  pdfVol(paste(output_dir,"/volcano/mix_css/",sep=""), Mix_ccs2[[as.character(x)]] ,x)
  if (x != 7){
    pdfVol(paste(output_dir,"/volcano/wt_css/",sep=""), WT_ccs2[[as.character(x)]] ,x)
    pdfVol(paste(output_dir,"/volcano/wt_mix/",sep=""), WT_Mix2[[as.character(x)]] ,x)
  }
}

############################ Section IV: Plots UMAPs: Fig4 ###########################################
# Plotting Definitions
DefaultAssay(seurat_integrated) <- "integrated"
Idents(object = seurat_integrated) <- "integrated_snn_res.0.1"
integrated_batch_markers <- FindAllMarkers(object = seurat_integrated, only.pos = TRUE, logfc.threshold = 0.25) 
metrics <-  c("nCount_RNA", "S.Score", "G2M.Score", "mitoRatio")

# Feature Plot - Sanity Check
pdfFN(paste(output_dir,"/featureplot.pdf",sep=""), FeaturePlot(seurat_integrated, 
                  reduction = "umap", features = metrics, pt.size = 0.4, order = TRUE, 
                  min.cutoff = 'q10', repel = TRUE, label = TRUE), 10, 8)

# Plot Single - Prep
seurat_integrated <- RenameIdents(object = seurat_integrated, "0" = "0 - Mesodermic Cells",
                                  "1" = "1 - Neural Stem Cells", "2" = "2 - Excitatory Neurons",
                                  "3" = "3 - Radial Glia", "4" = "4 - PFC Neurons - 1",
                                  "5" = "5 - Intermediate Progenitor Cells",
                                  "6" = "6 - PFC Neurons - 2", "7" = "7 - Endothelium")

hex_cols = c("#CD9600", "#7CAE00", "#00BE67", "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC",  "red")
single_plot = DimPlot(object = seurat_integrated, 
        reduction = "umap", label = TRUE,
        cols = hex_cols, label.size = 4,
        repel = TRUE) + NoLegend()
ggsave(file=paste(output_dir,"/Fig4_Seurat_Single.svg",sep=""), plot=single_plot, width=10, height=10, units = "cm")

# Plot Multiple
multi_plot = DimPlot(object = seurat_integrated, 
        reduction = "umap", label = FALSE,
        cols = hex_cols, split.by = "sample_type",
        repel = TRUE) + NoLegend()
ggsave(file=paste(output_dir,"/Seurat_Multi.svg",sep=""), plot=multi_plot, width=25, height=10, units = "cm")

# plot iendo, ineuron
seurat_integrated$sample_type <- factor(x=seurat_integrated$sample_type, levels = c("WT", "Mix", "CSS"))
Idents(object = seurat_integrated) <- "is_NGN1"
iendo = DimPlot(seurat_integrated,
        cols = c("gray70", "red", "#0c6c40"),
        reduction = "umap", label = FALSE,
        label.size = 3, split.by = "sample_type",
        ncol = 4) + NoLegend()
ggsave(file=paste(output_dir,"/iendo_ineuro.svg",sep=""), plot=iendo, width=25, height=10, units = "cm")
#show_col(hue_pal()(8))

#################################### Section 1: Find Markers for Clusters S13 ###################################
# Find heatmap markers 
DefaultAssay(seurat_integrated) <- "RNA"
Idents(object = seurat_integrated) <- "integrated_snn_res.0.1"
init_df <- FindConservedMarkers(seurat_integrated, ident.1 = 0, grouping.var = "batch",
                                only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
init_df$cluster <- 0
for (x in 1:7){
  temp_df <- FindConservedMarkers(seurat_integrated, ident.1 = x, grouping.var = "batch",
                                  only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  temp_df$cluster <- x
  init_df <- rbind (init_df, temp_df)
}

write.csv(init_df, paste(output_dir, "/FigS13_initdf_hm_markers_10.csv", sep=""))
gc()

## Scale the data before doing heatmap
DefaultAssay(seurat_integrated) <- "integrated"
init_df$gene <- rownames(init_df)
top10 <- init_df %>% group_by(cluster) %>% top_n(n = 10, wt = B1_avg_log2FC)

all.genes <- rownames(seurat_integrated)
hmp <- DoHeatmap(ScaleData(seurat_integrated, features = all.genes), features = top10$gene, slot="scale.data")
pg <- ggplot_build(hmp)
raw_hmpdata <- hmp$data
write.csv(raw_hmpdata, paste(output_dir, "/FigS13_raw_hmp_10.csv", sep=""))

pdfFN(paste(output_dir, "/FigS13_hm_markers_10.pdf", sep=""), hmp, 6 , 6)
gc()


