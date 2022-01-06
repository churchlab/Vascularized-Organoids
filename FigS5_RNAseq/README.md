# Vascularized_organoids
Repository for Vascularized_organoids Paper
Dr. Jeremy Y Huang - George M. Church Lab, Harvard University

Purpose: 
Flow-cytometry isolated endothelium and neuron co-culture bulk RNA sequencing. Figure S5 

./data contains salmon alignment quant.sf outputs 
./meta contains organoid metadata
./references contains .tsv files for elevated and enriched brain and vascular genes
	- needs txdb.gencode38.sqlite build from Gencode3.8 
./results contains latest build
./scripts contains de_organoid_script.R file for processing.


Most recent build

>sessionInfo()

R version 4.1.2 (2021-11-01)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19043)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252 
[2] LC_CTYPE=English_United States.1252   
[3] LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                          
[5] LC_TIME=English_United States.1252    

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils    
[7] datasets  methods   base     

other attached packages:
 [1] limma_3.48.3                gt_0.3.1                   
 [3] msigdbr_7.4.1               KEGGREST_1.32.0            
 [5] clusterProfiler_4.0.5       tximport_1.20.0            
 [7] enrichplot_1.12.3           GOSemSim_2.18.0            
 [9] org.Hs.eg.db_3.13.0         DOSE_3.18.3                
[11] ggrepel_0.9.1               DEGreport_1.28.0           
[13] biomaRt_2.48.3              data.table_1.14.2          
[15] RColorBrewer_1.1-2          forcats_0.5.1              
[17] stringr_1.4.0               dplyr_1.0.7                
[19] purrr_0.3.4                 readr_2.1.1                
[21] tidyr_1.1.4                 tibble_3.1.6               
[23] ggplot2_3.3.5               tidyverse_1.3.1            
[25] BiocParallel_1.26.1         GenomicFeatures_1.44.2     
[27] AnnotationDbi_1.54.1        DESeq2_1.32.0              
[29] SummarizedExperiment_1.22.0 Biobase_2.52.0             
[31] MatrixGenerics_1.4.3        matrixStats_0.61.0         
[33] GenomicRanges_1.44.0        GenomeInfoDb_1.28.4        
[35] IRanges_2.26.0              S4Vectors_0.30.0           
[37] BiocGenerics_0.38.0        

loaded via a namespace (and not attached):
  [1] utf8_1.2.2                  tidyselect_1.1.1           
  [3] RSQLite_2.2.9               grid_4.1.2                 
  [5] scatterpie_0.1.7            munsell_0.5.0              
  [7] codetools_0.2-18            withr_2.4.3                
  [9] colorspace_2.0-2            filelock_1.0.2             
 [11] knitr_1.37                  rstudioapi_0.13            
 [13] lasso2_1.2-22               GenomeInfoDbData_1.2.6     
 [15] mnormt_2.0.2                polyclip_1.10-0            
 [17] bit64_4.0.5                 farver_2.1.0               
 [19] pheatmap_1.0.12             downloader_0.4             
 [21] treeio_1.16.2               vctrs_0.3.8                
 [23] generics_0.1.1              xfun_0.29                  
 [25] BiocFileCache_2.0.0         R6_2.5.1                   
 [27] doParallel_1.0.16           clue_0.3-60                
 [29] graphlayouts_0.7.2          locfit_1.5-9.4             
 [31] bitops_1.0-7                cachem_1.0.6               
 [33] reshape_0.8.8               fgsea_1.18.0               
 [35] gridGraphics_0.5-1          DelayedArray_0.18.0        
 [37] assertthat_0.2.1            BiocIO_1.2.0               
 [39] scales_1.1.1                ggraph_2.0.5               
 [41] gtable_0.3.0                Cairo_1.5-14               
 [43] tidygraph_1.2.0             rlang_0.4.11               
 [45] genefilter_1.74.1           GlobalOptions_0.1.2        
 [47] splines_4.1.2               lazyeval_0.2.2             
 [49] rtracklayer_1.52.1          broom_0.7.10               
 [51] yaml_2.2.1                  reshape2_1.4.4             
 [53] modelr_0.1.8                backports_1.4.1            
 [55] qvalue_2.24.0               tools_4.1.2                
 [57] psych_2.1.9                 logging_0.10-108           
 [59] ggplotify_0.1.0             ellipsis_0.3.2             
 [61] ggdendro_0.1.22             Rcpp_1.0.7                 
 [63] plyr_1.8.6                  progress_1.2.2             
 [65] zlibbioc_1.38.0             RCurl_1.98-1.5             
 [67] prettyunits_1.1.1           GetoptLong_1.0.5           
 [69] viridis_0.6.2               cowplot_1.1.1              
 [71] haven_2.4.3                 cluster_2.1.2              
 [73] fs_1.5.2                    magrittr_2.0.1             
 [75] DO.db_2.9                   circlize_0.4.13            
 [77] reprex_2.0.1                tmvnsim_1.0-2              
 [79] patchwork_1.1.1             hms_1.1.1                  
 [81] xtable_1.8-4                XML_3.99-0.8               
 [83] readxl_1.3.1                gridExtra_2.3              
 [85] shape_1.4.6                 compiler_4.1.2             
 [87] shadowtext_0.1.0            crayon_1.4.2               
 [89] htmltools_0.5.2             ggfun_0.0.4                
 [91] tzdb_0.2.0                  geneplotter_1.70.0         
 [93] aplot_0.1.1                 lubridate_1.8.0            
 [95] DBI_1.1.2                   tweenr_1.0.2               
 [97] dbplyr_2.1.1                ComplexHeatmap_2.8.0       
 [99] MASS_7.3-54                 rappdirs_0.3.3             
[101] babelgene_21.4              Matrix_1.4-0               
[103] cli_3.1.0                   igraph_1.2.10              
[105] pkgconfig_2.0.3             GenomicAlignments_1.28.0   
[107] xml2_1.3.3                  foreach_1.5.1              
[109] ggtree_3.0.4                annotate_1.70.0            
[111] XVector_0.32.0              rvest_1.0.2                
[113] yulab.utils_0.0.4           digest_0.6.29              
[115] ConsensusClusterPlus_1.56.0 Biostrings_2.60.1          
[117] cellranger_1.1.0            fastmatch_1.1-3            
[119] tidytree_0.3.6              edgeR_3.34.1               
[121] restfulr_0.0.13             curl_4.3.2                 
[123] Rsamtools_2.8.0             rjson_0.2.20               
[125] lifecycle_1.0.1             nlme_3.1-153               
[127] jsonlite_1.7.2              viridisLite_0.4.0          
[129] fansi_0.5.0                 pillar_1.6.4               
[131] lattice_0.20-45             Nozzle.R1_1.1-1            
[133] fastmap_1.1.0               httr_1.4.2                 
[135] survival_3.2-13             GO.db_3.13.0               
[137] glue_1.6.0                  png_0.1-7                  
[139] iterators_1.0.13            bit_4.0.4                  
[141] ggforce_0.3.3               stringi_1.7.6              
[143] blob_1.2.2                  memoise_2.0.1              
[145] ape_5.6                    