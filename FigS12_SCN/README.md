# Vascularized_organoids
Repository for Vascularized_organoids Paper
Dr. Jeremy Y Huang - George M. Church Lab, Harvard University

Purpose:
Code for cluster identification using SingleCellNet

./data needs mat_primary_2.rda file data from published human cortical brain data detailed in Bhaduri, A. et al. Nature (2020)
./meta contains .rda file metadata for annotation based on published human cortical brain data detailed in Bhaduri, A. et al. Nature (2020)
./results contains all output from SCN algorithm including heatmaps and violin plots utilized in our paper.
./scripts contains both a shell-script for cluster-based calling utilized and .r script which was called

Needs: JH_org_sc_seurat_integrated_CCA_2batch_clustered.rds to run uploaded which can be provided from authors.

RunTime Note: Allocated 100GB ram during run.
Ran on R 4.1.0 with Seurat 4.0 and SCN 