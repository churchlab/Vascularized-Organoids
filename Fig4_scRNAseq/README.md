# Vascularized_organoids
Repository for Vascularized_organoids Paper
Dr. Jeremy Y Huang - George M. Church Lab, Harvard University

Purpose:
Processed Data in Seurat Object for Figure 4 and accompanying supplemental data for Skylar Scott, Huang, Lue, et.al

Scripts generates figures and data using data/meta folders. 
Runs natively in R 4.1.2 and during second submission in 4.1.0

./data folder needs Seurat dataset in .rds format (14.5 gb loaded), which can be provided from authors: 
./meta folder
  -  JH_org_sc_seurat_integrated_CCA_2batch_clustered_meta_with_NGN1.rds contains information for cells with NGN1 barcode for iNeuron identification.
  -  marker_list.csv contains known endothelial markers for supplemental heatmap plot.
  -  marker_list_neuro.csv contains known neural markers for supplemental heatmap plot.

./scripts/marker_analysis.r executes scripts.
28+ GB of RAM needed natively as Seurat file is 14.5GB in size (RNA, SCT, Integrated)
Recommend running on cluster to reproduce.

./results contains latest re-run of marker_analysis.r and outputs. concurs with published data/plots.



Confirmed Latest Run SessionInfo()

R version 4.1.2 (2021-11-01)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19043)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
[3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.1252    

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] scales_1.1.1          ggrepel_0.9.1         enrichplot_1.12.3     DOSE_3.18.3           clusterProfiler_4.0.5
 [6] org.Hs.eg.db_3.13.0   AnnotationDbi_1.54.1  IRanges_2.26.0        S4Vectors_0.30.0      Biobase_2.52.0       
[11] BiocGenerics_0.38.0   GOSemSim_2.18.0       biomaRt_2.48.3        DataCombine_0.2.21    svglite_2.0.0        
[16] SeuratObject_4.0.4    Seurat_4.0.6          heatmaply_1.3.0       viridis_0.6.2         viridisLite_0.4.0    
[21] plotly_4.10.0         pheatmap_1.0.12       RColorBrewer_1.1-2    forcats_0.5.1         stringr_1.4.0        
[26] dplyr_1.0.7           purrr_0.3.4           readr_2.1.1           tidyr_1.1.4           tibble_3.1.6         
[31] ggplot2_3.3.5         tidyverse_1.3.1      

loaded via a namespace (and not attached):
  [1] rappdirs_0.3.3         scattermore_0.7        bit64_4.0.5            multcomp_1.4-17       
  [5] irlba_2.3.5            data.table_1.14.2      rpart_4.1-15           KEGGREST_1.32.0       
  [9] RCurl_1.98-1.5         generics_0.1.1         metap_1.7              TH.data_1.1-0         
 [13] cowplot_1.1.1          RSQLite_2.2.9          shadowtext_0.1.0       RANN_2.6.1            
 [17] future_1.23.0          bit_4.0.4              tzdb_0.2.0             mutoss_0.1-12         
 [21] spatstat.data_2.1-2    webshot_0.5.2          xml2_1.3.3             lubridate_1.8.0       
 [25] httpuv_1.6.4           assertthat_0.2.1       hms_1.1.1              promises_1.2.0.1      
 [29] TSP_1.1-11             fansi_0.5.0            progress_1.2.2         dendextend_1.15.2     
 [33] dbplyr_2.1.1           readxl_1.3.1           tmvnsim_1.0-2          igraph_1.2.10         
 [37] DBI_1.1.2              htmlwidgets_1.5.4      spatstat.geom_2.3-1    ellipsis_0.3.2        
 [41] backports_1.4.1        deldir_1.0-6           vctrs_0.3.8            ROCR_1.0-11           
 [45] abind_1.4-5            cachem_1.0.6           withr_2.4.3            ggforce_0.3.3         
 [49] sctransform_0.3.2      treeio_1.16.2          prettyunits_1.1.1      goftest_1.2-3         
 [53] mnormt_2.0.2           cluster_2.1.2          ape_5.6                lazyeval_0.2.2        
 [57] crayon_1.4.2           pkgconfig_2.0.3        tweenr_1.0.2           GenomeInfoDb_1.28.4   
 [61] nlme_3.1-153           seriation_1.3.1        rlang_0.4.11           globals_0.14.0        
 [65] lifecycle_1.0.1        miniUI_0.1.1.1         sandwich_3.0-1         downloader_0.4        
 [69] registry_0.5-1         filelock_1.0.2         BiocFileCache_2.0.0    mathjaxr_1.4-0        
 [73] modelr_0.1.8           cellranger_1.1.0       polyclip_1.10-0        matrixStats_0.61.0    
 [77] lmtest_0.9-39          Matrix_1.4-0           aplot_0.1.1            zoo_1.8-9             
 [81] reprex_2.0.1           ggridges_0.5.3         png_0.1-7              bitops_1.0-7          
 [85] KernSmooth_2.23-20     Biostrings_2.60.1      blob_1.2.2             qvalue_2.24.0         
 [89] parallelly_1.30.0      gridGraphics_0.5-1     memoise_2.0.1          magrittr_2.0.1        
 [93] plyr_1.8.6             ica_1.0-2              zlibbioc_1.38.0        compiler_4.1.2        
 [97] scatterpie_0.1.7       plotrix_3.8-2          fitdistrplus_1.1-6     cli_3.1.0             
[101] XVector_0.32.0         listenv_0.8.0          patchwork_1.1.1        pbapply_1.5-0         
[105] MASS_7.3-54            mgcv_1.8-38            tidyselect_1.1.1       stringi_1.7.6         
[109] grid_4.1.2             fastmatch_1.1-3        tools_4.1.2            future.apply_1.8.1    
[113] rstudioapi_0.13        foreach_1.5.1          gridExtra_2.3          farver_2.1.0          
[117] Rtsne_0.15             ggraph_2.0.5           digest_0.6.29          shiny_1.7.1           
[121] Rcpp_1.0.7             broom_0.7.10           later_1.3.0            RcppAnnoy_0.0.19      
[125] httr_1.4.2             Rdpack_2.1.3           colorspace_2.0-2       rvest_1.0.2           
[129] XML_3.99-0.8           fs_1.5.2               tensor_1.5             reticulate_1.22       
[133] splines_4.1.2          uwot_0.1.11            yulab.utils_0.0.4      tidytree_0.3.6        
[137] sn_2.0.1               spatstat.utils_2.3-0   graphlayouts_0.7.2     multtest_2.48.0       
[141] ggplotify_0.1.0        systemfonts_1.0.3      xtable_1.8-4           jsonlite_1.7.2        
[145] ggtree_3.0.4           tidygraph_1.2.0        ggfun_0.0.4            R6_2.5.1              
[149] TFisher_0.2.0          pillar_1.6.4           htmltools_0.5.2        mime_0.12             
[153] glue_1.6.0             fastmap_1.1.0          BiocParallel_1.26.1    codetools_0.2-18      
[157] fgsea_1.18.0           mvtnorm_1.1-3          utf8_1.2.2             lattice_0.20-45       
[161] spatstat.sparse_2.1-0  numDeriv_2016.8-1.1    curl_4.3.2             leiden_0.3.9          
[165] GO.db_3.13.0           limma_3.48.3           survival_3.2-13        munsell_0.5.0         
[169] DO.db_2.9              GenomeInfoDbData_1.2.6 iterators_1.0.13       haven_2.4.3           
[173] reshape2_1.4.4         gtable_0.3.0           rbibutils_2.2.7        spatstat.core_2.3-2   