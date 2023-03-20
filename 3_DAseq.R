#Here, the idea is that the DAseq algorithm both resembles ours more, and it
#uses a more conservative method of defining DA regions, so it might be 
#more amenable to statistical findings actually. 
library(DAseq)
library(SingleCellExperiment)

metaData <- readRDS("Data/metaData.rds")
umapData <- readRDS("Data/umapData.rds")
logCountsRed <- readRDS("Data/rnaScaleData.rds")
harmonyData <- readRDS("Data/harmonyData.rds")

python2use <- "/Users/jakob.theorell/Library/r-miniconda/envs/r-reticulate/bin/python"
GPU <- 8

allSce <- SingleCellExperiment(assays = list(logcounts = logCountsRed), 
                               colData = metaData)

reducedDim(allSce, "Harmony") <- harmonyData
reducedDim(allSce, "UMAP") <- umapData

donInflam <- paste0(metaData$donor1, metaData$Tissue)
donInflam <- gsub("-", "", donInflam)
highLowSce <- allSce[,which(donInflam %in% c( "D1Noninflamed", "D16Noninflamed", 
                                              "D6Inflamed", "D19Noninflamed",  
                                              "D19Inflamed","D10Inflamed", 
                                              "D11Noninflamed", "D11Inflamed"))]

highLowSce$Inflam <- "Low"
highLowSce$Inflam[which(highLowSce$CHASTea > 3)] <- "High"
table(highLowSce$Inflam)
#High   Low 
#12701 15515 

label_low <- unique(highLowSce$donor[which(highLowSce$Inflam == "Low")])
label_high <- unique(highLowSce$donor[which(highLowSce$Inflam == "High")])

daCells <- getDAcells(reducedDim(highLowSce, "Harmony"), 
                      highLowSce$donor, label_low, label_high, 
                      plot.embedding = reducedDim(highLowSce, "UMAP"))

saveRDS(daCells, "Data_files/daCells.rds")

dir.create("Results/DAseq")

pdf("Results/DAseq/Prediction_plot.pdf")
daCells$pred.plot
dev.off()

daCells2 <- updateDAcells(
  X = daCells, pred.thres = c(-0.5,0.5),
  plot.embedding = reducedDim(highLowSce, "UMAP")
)

pdf("Results/DAseq/DA_cells.pdf")
daCells2$da.cells.plot
dev.off()

daRegions <- getDAregion(
  X = reducedDim(highLowSce, "Harmony"),
  da.cells = daCells2,
  cell.labels = highLowSce$donor,
  labels.1 = label_low,
  labels.2 = label_high, resolution = 0.01,
  plot.embedding = reducedDim(highLowSce, "UMAP"),
  min.cell = 50
)

#Here, we check which clusters that are significantly different between the groups
daRegions$DA.stat
#DA.score pval.wilcoxon pval.ttest
#[1,]  0.5116863    0.02857143 0.06531729
#[2,]  0.6830990    0.02857143 0.03268729
#[3,] -0.7294222    0.02857143 0.15625757
#[4,] -0.6437529    0.68571429 0.21973425

#And with FDR correction: 
p.adjust(daRegions$DA.stat[,2], method = "fdr")
#0.03809524 0.03809524 0.03809524 0.68571429

pdf("Results/DAseq/DA_cell_regions.pdf")
daRegions$da.region.plot
dev.off()

daRegions2 <- daRegions
daRegions2$da.region.label[which(daRegions$da.region.label == 4)] <- 0

pdf("Results/DAseq/DA_cell_regions_significant.pdf")
daRegions2$da.region.plot
dev.off()

highLowSce$daRegions <- daRegions$da.region.label

#And we remove label 4, as it is non-significant. 

