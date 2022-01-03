library(singleCellNet)
library(Seurat)

######################### Section I: Load Training Data and Seurat DATA ########################
#loading training data
stPrim <- utils_loadObject(fname = "meta_primary.rda")  #metadata
stPrim<-stPrim[!(stPrim$Subtype=="Outlier"),] # Remove outlier from calculus
expPrim <- utils_loadObject(fname = "mat_primary_2.rda")  #expression matrix
stPrim <- droplevels(stPrim)

#loading query data
clustered_seurat <- readRDS("JH_org_sc_seurat_integrated_CCA_2batch_clustered.rds")
stOrganoid = clustered_seurat@meta.data         #metadata
expOrganoid = clustered_seurat@assays$RNA@counts                        #expression
stOrganoid$cell <- rownames(stOrganoid) #Add rownames as column "cell"
genesOrganoid = rownames(expOrganoid)

# Find genes in common to the datasets and limit analysis to these
commonGenes = intersect(rownames(expPrim), genesOrganoid)
expPrim = expPrim[commonGenes,]

######################### Section II: Begin Training Dataset & QC ###############################

# Split for training and assessment, and transform training data
set.seed(100) #can be any random seed number
stList = splitCommon(sampTab=stPrim, ncells=100, dLevel="Subtype")
stTrain = stList[[1]]
expTrain = expPrim[,rownames(stTrain)]

system.time(class_info<-scn_train(stTrain = stTrain, expTrain = expTrain, nTopGenes = 75, 
	nRand = 70, nTrees = 1000, nTopGenePairs = 150, dLevel = "Subtype", colName_samp = "cell"))

#validate data
stTestList = splitCommon(sampTab=stList[[2]], ncells=100, dLevel="Subtype") #normalize validation data so that the assessment is as fair as possible
stTest = stTestList[[1]]
expTest = expPrim[commonGenes,rownames(stTest)]

#predict
classRes_val_all = scn_predict(cnProc=class_info[['cnProc']], expDat=expTest, nrand = 70)
tm_heldoutassessment = assess_comm(ct_scores = classRes_val_all, stTrain = stTrain, 
	stQuery = stTest, dLevelSID = "cell", classTrain = "Subtype", classQuery = "Subtype", nRand = 70)

# Plot Dot Plots, EMAPS
pdfFN <- function(path, object){
  pdf(path, width = 7, height = 11)
  print(object)
  dev.off()
}

#Create a name vector label used later in classification heatmap where the values are cell types/ clusters and names are the sample names
nrand = 70
sla = as.vector(stTest$Subtype)
names(sla) = as.vector(stTest$cell)
slaRand = rep("rand", nrand) 
names(slaRand) = paste("rand_", 1:nrand, sep='')
sla = append(sla, slaRand) #include in the random cells profile created

# Assess Clarifier
pdfFN("tm_heldoutassessment.pdf", plot_PRs(tm_heldoutassessment))
pdfFN("tm_heldoutassessment.pdf", plot_metrics(tm_heldoutassessment))

# Classification result heatmap
pdfFN("sc_hmClass.pdf", sc_hmClass(classMat = classRes_val_all,grps = sla, max=300, isBig=TRUE))

# Attribution Plot
pdfFN("plot_attr.pdf",plot_attr(classRes=classRes_val_all, sampTab=stTest, nrand=nrand, dLevel="Subtype", sid="cell"))

# Visualize average top pairs genes expression for training data

gpTab = compareGenePairs(query_exp = expTest, training_exp = expTrain, training_st = stTrain, classCol = "Subtype", sampleCol = "cell", RF_classifier = class_info$cnProc$classifier, numPairs = 20, trainingOnly= TRUE)
train = findAvgLabel(gpTab = gpTab, stTrain = stTrain, dLevel = "Subtype")
pdfFN("top_pairs.pdf", hm_gpa_sel(gpTab, genes = class_info$cnProc$xpairs, grps = train, maxPerGrp = 50))

############################## SECTION III: Apply Trained Datasets to Organoid Data ##########################

nqRand = 70
system.time(crOrganoidall<-scn_predict(class_info[['cnProc']], expOrganoid, nrand=nqRand))

sgrp = as.vector(stOrganoid$integrated_snn_res.0.1)
names(sgrp) = as.vector(stOrganoid$cell)
grpRand =rep("rand", nqRand)
names(grpRand) = paste("rand_", 1:nqRand, sep='')
sgrp = append(sgrp, grpRand)

# This classifies a cell with  the catgory with the highest classification score or higher than a classification score threshold of your choosing.
# The annotation result can be found in a column named category in the query sample table.
stOrganoid <- get_cate(classRes = crOrganoidall, sampTab = stOrganoid, dLevel = "integrated_snn_res.0.1", sid = "cell", nrand = nqRand)

# heatmap classification result
pdfFN("hm_Organoid_all.pdf", sc_hmClass(crOrganoidall, sgrp, max=5000, isBig=TRUE, cCol=F, font=8))
pdfFN("violin_Organoid_all.pdf", sc_violinClass(sampTab = stOrganoid, classRes = crOrganoidall, sid = "cell", dLevel = "integrated_snn_res.0.1", addRand = nqRand))

saveRDS(crOrganoidall, file="crOrganoidall.rds")
saveRDS(expTrain, file="expTrain.rds")
saveRDS(stTrain, file="stTrain.rds")
saveRDS(stOrganoid, file="stOrganoid.rds")
write.csv(t(crOrganoidall),"crOrganoidall.csv")
write.csv(stOrganoid,"stOrganoid.csv")