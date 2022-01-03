### Bioconductor and CRAN libraries used
# Note current version of Gencode is 3.8
#update.packages(ask=FALSE)
set.seed(0)  # random seed for reproducibility
rm(list = ls())
load_pkg <- rlang::quos(DESeq2, GenomicFeatures, BiocParallel, tidyverse, RColorBrewer, data.table, biomaRt,
                        DEGreport, ggplot2, ggrepel, DOSE, org.Hs.eg.db, dplyr, GOSemSim, enrichplot, tximport,
                        AnnotationDbi, clusterProfiler, KEGGREST, msigdbr, gt, limma)

invisible(lapply(lapply(load_pkg, rlang::quo_name), library, character.only = TRUE))
`%nin%` = Negate(`%in%`)

############################## Global Variables ###############################
#Set Thresholds for p-value and log2fold threshold
ont_list = c("CC", "MF", "BP")
padj.cutoff <- 0.05
lfc.cutoff <- 0.58  

output_dir = paste("../results/data_analysis_organoid_", format(Sys.Date(), "%Y%m%d"), sep="")
dir.create(output_dir)

# Plot Dot Plots
pdfFN <- function(path, object, w, h){
  pdf(path, width = w, height = h)
  print(object)
  dev.off()
}

################################ Section I: Load SALMON Datasets #######################################
## List all directories containing data_all & obtain a vector of all filenames including the path 
samples <- list.files(path = "../data/salmon", recursive = T, include.dirs = T, full.names = T, pattern="salmon$")
files <- file.path(samples, "quant.sf")
names(files) <- str_replace(samples, "../data/salmon", "") %>% str_replace("_concat.salmon", "") %>%
  str_replace("LIB052393_", "") %>% gsub(".*_", "", .)

txdb <- loadDb('../references/txdb.gencode38.sqlite')
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID","TXNAME") 
txi <- tximport(files, type="salmon", tx2gene=tx2gene, countsFromAbundance="lengthScaledTPM")
data_all <- txi$counts %>% round() %>% data.frame()

# Get ENSG Gene_IDs from Biomart
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene_IDs <- getBM(attributes = c("external_gene_name","ensembl_gene_id", "entrezgene_id"),
                  filters = "ensembl_gene_id", values =  gsub("\\..*", "",rownames(data_all)), mart = mart) 

meta <- column_to_rownames(read.csv("../meta/organoid_meta.csv", fileEncoding="UTF-8-BOM"), var="Library_ID")
neuro_elevated <- column_to_rownames(read.delim("../references/brain_elevated_proteinatlas.tsv", fileEncoding="UTF-8-BOM"), var="Ensembl")
endo_elevated <- column_to_rownames(read.delim("../references/vascular_elevated_proteinatlas.tsv", fileEncoding="UTF-8-BOM"), var="Ensembl")
colnames(data_all) <- rownames(meta)
rownames(data_all) <- gsub("\\..*", "",rownames(data_all))

meta_endo <- filter(meta, Tissue_Type == "Endo" | Tissue_Type == "Endo.WT" | Tissue_Type == "Endo.Neuron")
data_endo <- data_all[rownames(data_all) %in% rownames(neuro_elevated) == FALSE,]
data_endo <- data_endo[,colnames(data_endo) %in% rownames(meta_endo) == TRUE]

meta_neuron <- filter(meta, Tissue_Type == "Neuron" | Tissue_Type == "Neuron.WT" | Tissue_Type == "Neuron.Endo")
data_neuron <- data_all[rownames(data_all) %in% rownames(endo_elevated) == FALSE,]
data_neuron <- data_neuron[,colnames(data_neuron) %in% rownames(meta_neuron) == TRUE]

# Tissue_Types & data_all
all_Tissues <- c("Endo", "Endo.WT", "Endo.Neuron", "Neuron", "Neuron.WT", "Neuron.Endo")
endo_Tissues <- c("Endo", "Endo.WT", "Endo.Neuron")
neuron_Tissues <- c("Neuron", "Neuron.WT", "Neuron.Endo")

####################### SECTION II: Count Normalization #################

AnalysisFN <- function(data, meta, output_dir, tissues){
  # Create subdirectories
  dir.create(paste(output_dir, "/deseq_results", sep=""), recursive=TRUE)
  dir.create(paste(output_dir, "/pca", sep=""), recursive=TRUE)
  dir.create(paste(output_dir, "./sig_results/", sep = ""), recursive = TRUE)
  dir.create(paste(output_dir, "/ENSEMBL_GO/", sep = ""), recursive = TRUE)
  dir.create(paste(output_dir, "/volcano/", sep = "") , recursive = TRUE)
  
  ## Create DESeq2Dataset object.
  dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ Tissue_Type)
  dds <- estimateSizeFactors(dds)
  idx <- rowSums(counts(dds, normalized=TRUE)>=5)>=3
  dds <- dds[idx,]
  normalized_counts <- counts(dds, normalized=TRUE)
  write.table(normalized_counts, file=paste(output_dir, "/normalized_salmon_counts.txt", sep=""), sep="\t", quote=F, col.names=NA)
  
  # Run Deseq2 & make VST for PCA   # Returns list of paired DEseq results and a character list of pairs used.
  dds <- DESeq(dds)
  rld <- vst(dds, blind=TRUE)
  data_deseqResult = list()
  pairs_list = list()
  
  # PlotPCA
  pcaData <- plotPCA(rld[ , rld$Tissue_Type %in% tissues ], 
                     intgroup=c("Tissue_Type", "Replicate_Num"), returnData=TRUE)
  pdfFN(paste(output_dir, "/pca/pca_tissuetype.pdf", sep="") , 
        ggplot(pcaData, aes(PC1, PC2, color=Tissue_Type, shape=Replicate_Num)) +
          geom_point(size=3) +
          coord_fixed(), 6, 5)
  write.csv(pcaData, paste(output_dir, "/pca/all_pca.csv", sep=""), row.names = TRUE)
  
  # Find Significant Genes for GO
  for (p in 1:(length(tissues)-1)){
    ind_1 = p+1
    for (x in ind_1:length(tissues)) {
      pair_name = paste(tissues[p],tissues[x], sep="_")
      contrast_nam <- c("Tissue_Type", tissues[p], tissues[x])
      res_table_temp <- results(dds, contrast=contrast_nam, alpha = 0.05)         # Type S4
      res_table_temp <- lfcShrink(dds, contrast=contrast_nam, res=res_table_temp, type="ashr")
      res_table_temp <- setDT(res_table_temp %>% data.frame(), keep.rownames=T)   # Type data.frame
      names(res_table_temp)[1] <- "ENSG" 
      res_table_temp$ENSG_noVer <- gsub("\\.[^.]*$","", res_table_temp$ENSG)
      res_table_temp <- left_join(res_table_temp, gene_IDs, by = c("ENSG_noVer"="ensembl_gene_id"))
      res_table_temp <- subset(res_table_temp, select=c("ENSG","external_gene_name",names(res_table_temp)[2:6], "entrezgene_id"))
      res_table_temp <- res_table_temp[!duplicated(res_table_temp[,"ENSG"]),]
      data_deseqResult[[pair_name]] <- res_table_temp
      pairs_list <- append(pairs_list, paste(tissues[p],tissues[x], sep="_"))
      write.table(res_table_temp, file = paste(output_dir, "/deseq_results/", pair_name, "_deseq", ".txt", sep=""), sep = "\t", row.names = T, col.names = T)
    }
  }

  # Returns a list of tibbles for significant genes based on cutoff, filter for padj, and lfc.cutoff 
  data_sigResult  = list()
  for (x in pairs_list){
    res_tableOE_tb <- data_deseqResult[[x]] %>% as_tibble()
    data_sigResult[[x]] <- res_tableOE_tb %>% filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff) %>% arrange(desc(log2FoldChange))
    write.table(data_sigResult[[x]], file = paste(output_dir, "/sig_results/", x, "_sig", ".txt", sep=""), sep = "\t", row.names = T, col.names = T)
  }

########################### Section V: GO/KEGG Term Analysis #########################

  
  for (x in ont_list) {dir.create(paste(output_dir, "/ENSEMBL_GO/", x, "/", sep = ""))}
  screened_sigresult = list()
  screened_universeresult = list()
  GO_Ensembl_list = list()
  for (x in pairs_list){
    temp <- data_sigResult[[x]]
    screened_sigresult[[x]] <- temp[!(is.na(temp$external_gene_name)),]  %>% arrange(padj)
    temp2 <- data_deseqResult[[x]]
    screened_universeresult[[x]] <- temp2[!(is.na(temp2$external_gene_name)),]
  }
  
  # Run GO all for all paired lists using Ensembl
  for (x in pairs_list){
    allOE_genes_ensembl <- gsub("\\.[^.]*$","", as.character(screened_universeresult[[x]]$ENSG))
    sigOE_genes_ensembl <- gsub("\\.[^.]*$","",as.character(screened_sigresult[[x]]$ENSG))
    for (y in ont_list){
      ego_ensembl <- enrichGO(gene = sigOE_genes_ensembl, universe = allOE_genes_ensembl, keyType = "ENSEMBL", 
                              OrgDb = org.Hs.eg.db, ont = y, pAdjustMethod = "BH", readable = TRUE, qvalueCutoff = 0.05)
      if(!is.null(ego_ensembl)){
        if (!nrow(ego_ensembl)==0) {
          ego_ensembl <- pairwise_termsim(ego_ensembl, method = "Wang", semData = godata('org.Hs.eg.db', ont=y))
          }
        save(ego_ensembl, file=paste(output_dir, "/ENSEMBL_GO/", y, "/", x, "_ego_ensembl.rda", sep=""))
        write.csv(ego_ensembl, paste(output_dir, "/ENSEMBL_GO/", y, "/", x, "_EnsemblGO.csv", sep=""))
        ID = paste(x, y, sep="_")
        GO_Ensembl_list[[ID]] <- ego_ensembl
      }
    }
  }
  #  Filter based on Sig Results Print Volcano
  for (x in pairs_list){
    res_tableOE_tb <- data_deseqResult[[x]][complete.cases(data_deseqResult[[x]]),] %>% drop_na(entrezgene_id) # %>% filter(abs(log2FoldChange)<10)
    res_tableOE_tb <- res_tableOE_tb %>% mutate(threshold_OE = padj < padj.cutoff & abs(log2FoldChange) >= lfc.cutoff)
    res_tableOE_tb <- res_tableOE_tb[!grepl("LINC*", res_tableOE_tb$external_gene_name),]
    res_tableOE_tb <- res_tableOE_tb %>% arrange(padj) %>% mutate(genelabels = "")
    res_tableOE_tb$genelabels[1:15] <- res_tableOE_tb$external_gene_name[1:15]
    pdfFN(paste(output_dir, "/volcano/" ,x,"_volcano.pdf", sep=""),
          ggplot(res_tableOE_tb, aes(x = log2FoldChange, y = -log10(padj))) +
            geom_point(aes(colour = threshold_OE)) +
            geom_text_repel(aes(label = genelabels), max.overlaps = Inf) +
            ggtitle(paste(x, "Volcano Plot", sep=" ")) +
            xlab("log2 fold change") + 
            xlim(-8, 8) +
            ylab("-log10 adjusted p-value") +
            geom_hline(yintercept=-log10(padj.cutoff), linetype="dashed", color = "grey") +
            geom_vline(xintercept=-lfc.cutoff, linetype="dashed", color = "grey") +
            geom_vline(xintercept=lfc.cutoff, linetype="dashed", color = "grey") +
            theme(legend.position = "none",
                  plot.title = element_text(size = rel(1.5), hjust = 0.5),
                  axis.title = element_text(size = rel(1.25))) 
          ,5, 8)
  }
  
  ## EMAPS&GO&Netplot Print
  color_pal <- brewer.pal(3,"RdBu")
  
  for (x in pairs_list){
    for (y in ont_list){
      if(!is.null(GO_Ensembl_list[[paste(x,y,sep="_")]])){
        go_Ensembl <- GO_Ensembl_list[[paste(x,y,sep="_")]]
        go_Ensembl@result$ID <- factor(go_Ensembl@result$ID, levels = go_Ensembl@result$ID[order(go_Ensembl@result$p.adjust)])
        pair_OE <- screened_sigresult[[x]] # %>% filter(abs(log2FoldChange)<10)
        OE_foldchanges <- pair_OE$log2FoldChange %>% as.vector()
        names(OE_foldchanges) <-pair_OE$external_gene_name
        print(paste(x, y))
        if (nrow(go_Ensembl)>=2){
          pdfFN(paste(output_dir, "/ENSEMBL_GO/",y,"/",x,"_EnsemblGO.pdf", sep=""), dotplot(go_Ensembl, showCategory=10,  orderBy = "x") + scale_y_discrete(labels=function(x) str_wrap(x, width=45)), 8, 6) 
        }
      }
    }
  }
}
dir.create(paste(output_dir, "/all", sep=""), recursive = TRUE)
dir.create(paste(output_dir, "/endo", sep=""), recursive = TRUE)
dir.create(paste(output_dir, "/neuron", sep=""), recursive = TRUE)
AnalysisFN(data_endo, meta_endo, paste(output_dir, "/endo", sep=""), endo_Tissues)
AnalysisFN(data_neuron, meta_neuron, paste(output_dir, "/neuron", sep=""), neuron_Tissues)
AnalysisFN(data_all, meta, paste(output_dir, "/all", sep=""), all_Tissues)
