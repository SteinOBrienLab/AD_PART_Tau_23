## 1/14/2023 Ryan Palaganas
#read in data from GeoMx files-----
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
# Reference the main folder 'file.path' containing the sub-folders with each----
datadir <- file.path("~/Directory/NanostringTauProject/")
# the DataLocation folder should contain a dccs, pkcs, and annotation folder
# with each set of files present as needed
# automatically list files in each directory for use
DCCFiles <- dir(file.path(datadir, "dccs"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)
PKCFiles <- dir(file.path(datadir, "pkcs"), pattern = ".pkc$",
                full.names = TRUE, recursive = TRUE)
AnnotationFile <-
  dir(file.path(datadir, "annotation"), pattern = ".xlsx$",
      full.names = TRUE, recursive = TRUE)
myData <-
  readNanoStringGeoMxSet(dccFiles = DCCFiles,
                         pkcFiles = PKCFiles,
                         phenoDataFile = AnnotationFile,
                         phenoDataSheet = "Annotation",
                         phenoDataDccColName = "Sample_ID",
                         protocolDataColNames = c("aoi", "roi"),
                         experimentDataColNames = c("panel"))

saveRDS(myData,"~/Directory/NanostringTauProject/cache/myData.rds")
#check data is in place with correct pkc
library(knitr)
pkcs <- annotation(myData)
modules <- gsub(".pkc", "", pkcs)
kable(data.frame(PKCs = pkcs, modules = modules))


#looking at the data----
library(dplyr)
library(ggforce)
# select the annotations we want to show, use `` to surround column names with
# spaces or special symbols
count_mat <- count(pData(myData), `slide_name`, class, tau, `group`, Label, sex, age, PMI,collection)
#fixes error due to combination of character and double in count_mat[1] and count_mat[4]
count_mat$group <- as.character(count_mat$group)
# gather the data and plot in order: class, slide name, region, segment
test_gr <- gather_set_data(count_mat, 1:4)

#test_gr$x <- factor(test_gr$x,
                    #levels = c("class", "slide_name", "tau", "group", "Label", "sex", "age", "PMI", "collection"))
ggplot(test_gr, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = tau), alpha = 0.5, axis.width = 0.1) +
  geom_parallel_sets_axes(axis.width = 0.2) +
  geom_parallel_sets_labels(color = "white", size = 5) +
  theme_classic(base_size = 17) + 
  theme(legend.position = "bottom",
        axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_blank()) +
  scale_y_continuous(expand = expansion(0)) + 
  scale_x_discrete(expand = expansion(0)) +
  labs(x = "", y = "") 

# Shift counts to one
head(exprs(myData))
range(exprs(myData)) #0 42354


myDataS<- shiftCountsOne(myData, useDALogic = TRUE)

range(exprs(myDataS)) 

myData <- shiftCountsOne(myData, useDALogic = TRUE)

#QC process----
QC_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 75,    # Minimum % of reads aligned (80%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 1,   # Minimum negative control counts (10)
       maxNTCCount = 9000,     # Maximum counts observed in NTC well (1000)
       minNuclei = 0,         # Minimum # of nuclei estimated (100) - none in this dataset since not counting nuclei well
       minArea = 2000)         # Minimum segment area (5000)

myData <-
  setSegmentQCFlags(myData, 
                    qcCutoffs = QC_params)  
print((protocolData(myData)[["QCFlags"]]))

# Collate QC Results
QCResults <- protocolData(myData)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))

#visualize QC
library(ggplot2)
col_by <- "group"
# Graphical summaries of QC statistics plot function
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
}

QC_histogram(sData(myData), "Trimmed (%)", col_by, 80)
ggsave(paste0("~/Directory/NanostringTauProject/plots/trimmed.pdf"))
QC_histogram(sData(myData), "Stitched (%)", col_by, 80)
ggsave(paste0("~/Directory/NanostringTauProject/plots/stitched.pdf"))
QC_histogram(sData(myData), "Aligned (%)", col_by, 75)
ggsave(paste0("~/Directory/NanostringTauProject/plots/aligned.pdf"))
QC_histogram(sData(myData), "Saturated (%)", col_by, 50) +
  labs(title = "Sequencing Saturation (%)",
       x = "Sequencing Saturation (%)")
ggsave(paste0("~/Directory/NanostringTauProject/plots/sequencingsaturation.pdf"))
QC_histogram(sData(myData), "area", col_by, 1000, scale_trans = "log10")
ggsave(paste0("~/Desktop/Directory/NanostringTauProject/plots/area.pdf"))
# calculate the negative geometric means for each module
negativeGeoMeans <- 
  esBy(negativeControlSubset(myData), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       })
protocolData(myData)[["NegGeoMean"]] <- negativeGeoMeans
# explicitly copy the Negative geoMeans from sData to pData
negCols <- paste0("NegGeoMean_", modules)
pData(myData)[, negCols] <- sData(myData)[["NegGeoMean"]]

for(ann in negCols) {
  ngmgroup <- QC_histogram(pData(myData), ann, 'group', 2, scale_trans = "log10")
  print(ngmgroup)
}
ngmgroup
ggsave(paste0("~/Directory/NanostringTauProject/plots/ngmgroup.pdf"))
#look by tau
for(ann in negCols) {
  ngmtau <- QC_histogram(pData(myData), ann, 'tau', 2, scale_trans = "log10")
  print(ngmtau)
}
ngmtau
ggsave(paste0("~/Directory/NanostringTauProject/plots/ngmtau.pdf"))

# detatch neg_geomean columns ahead of aggregateCounts call
pData(myData) <- pData(myData)[, !colnames(pData(myData)) %in% negCols]

# show all NTC values, Freq = # of Segments with a given NTC count:
kable(table(NTC_Count = sData(myData)$NTC),
      col.names = c("NTC Count", "# of Segments"))
kable(QC_Summary, caption = "QC Summary Table for each Segment")
myData <- myData[, QCResults$QCStatus == "PASS"]
# Subsetting our dataset has removed samples which did not pass QC
dim(myData)

#probeQC to remove outlier NC probes
# Generally keep the qcCutoffs parameters unchanged. Set removeLocalOutliers to 
# FALSE if you do not want to remove local outliers
myData <- setBioProbeQCFlags(myData, 
                             qcCutoffs = list(minProbeRatio = 0.1,
                                              percentFailGrubbs = 20), 
                             removeLocalOutliers = TRUE)

ProbeQCResults <- fData(myData)[["QCFlags"]]

# Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))
qc_df
#Subset object to exclude all that did not pass Ratio & Global testing
ProbeQCPassed <- 
  subset(myData, 
         fData(myData)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(myData)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dim(ProbeQCPassed)

myData <- ProbeQCPassed 

#count level data assessment
# Check how many unique targets the object has
length(unique(featureData(myData)[["TargetName"]]))
# collapse to targets
target_myData <- aggregateCounts(myData)
dim(target_myData)
exprs(target_myData)[1:5, 1:2]

#LOQ determination with segment and gene filtration----
# Define LOQ SD threshold and minimum value
cutoff <- 2
minLOQ <- 2

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_myData))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(target_myData)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(target_myData)[, vars[1]] * 
             pData(target_myData)[, vars[2]] ^ cutoff)
  }
}
pData(target_myData)$LOQ <- LOQ

#filter data by segment
LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(target_myData)$Module == module
  Mat_i <- t(esApply(target_myData[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(target_myData)$TargetName, ]

# Save detection rate information to pheno data
pData(target_myData)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)
pData(target_myData)$GeneDetectionRate <-
  pData(target_myData)$GenesDetected / nrow(target_myData)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_myData)$DetectionThreshold <- 
  cut(pData(target_myData)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.02, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-2%", "2-5%", "5-10%", "10-15%", ">15%"))
# stacked bar plot of different cut points (1%, 5%, 10%, 15%)

pdf("~/Directory/NanostringTauProject/plots/GeneDetectionRatePerSegmentNo.pdf")
ggplot(pData(target_myData),
       aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = tau)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_fill_manual(values = c('negative' = 'royalblue3', 'positive' = 'red3')) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Tau")
dev.off()


# cut percent genes detected at 1, 5, 10, 15
kable(table(pData(target_myData)$DetectionThreshold,
            pData(target_myData)$group))

kable(table(pData(target_myData)$DetectionThreshold,
            pData(target_myData)$slide_name))

target_myData <-
  target_myData[, pData(target_myData)$GeneDetectionRate >= .02]

dim(target_myData)

# select the annotations we want to show, use `` to surround column names with
# spaces or special symbols
count_mat <- count(pData(target_myData), `slide_name`, class, tau, group, Label, sex, age, PMI, collection)

# gather the data and plot in order: class, slide name, region, segment
test_gr <- gather_set_data(count_mat, 1:4)
#test_gr$x <-
  #factor(test_gr$x,
         #levels = c("class", "slide_name", "tau", "group"))
# plot Sankey
ggplot(test_gr, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = tau), alpha = 0.5, axis.width = 0.1) +
  geom_parallel_sets_axes(axis.width = 0.2) +
  geom_parallel_sets_labels(color = "white", size = 5) +
  theme_classic(base_size = 17) + 
  theme(legend.position = "bottom",
        axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_blank()) +
  scale_y_continuous(expand = expansion(0)) + 
  scale_x_discrete(expand = expansion(0)) +
  labs(x = "", y = "") 
ggsave("~/Directory/NanostringTauProject/plots/sankey.pdf")
dev.off()
library(scales) # for percent

# Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(target_myData)]
fData(target_myData)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_myData)$DetectionRate <-
  fData(target_myData)$DetectedSegments / nrow(pData(target_myData))


# Gene of interest detection table
goi <- c("GRIA1", "GRIA2", "GRIA3", "GRIA4", "GRIN1", "GRIN2A","ARC",
         "GRIN2B", "BDNF", "HOMER1", "WWC1", "NPTX2", "SYN1", "SYP")
goi_df <- data.frame(
  Gene = goi,
  Number = fData(target_myData)[goi, "DetectedSegments"],
  DetectionRate = percent(fData(target_myData)[goi, "DetectionRate"]))

cellMarker <- c("SNAP25", "SLC12A5", "SYT1", "NEFM", "GABRA1","RBFOX3", 
                "GAD1","PVALB", "SST", "NEUROD6","THY1", "CAMK2A", "MBP", 
                "RANGRF", "GPAT3", "SOX10","GJC2","SLC1A2","GFAP","FGFR3",
                "GJB6","AQP4", "ALDH1L1", "TMEM119","GPR34","CD14","AIF1")
cellMarker_df <- data.frame(
  Gene = cellMarker,
  Number = fData(target_myData)[cellMarker, "DetectedSegments"],
  DetectionRate = percent(fData(target_myData)[cellMarker, "DetectionRate"]))

litMarker <- c("PRNP", "APP","CANX","NCAM1","ARHGAP35","SARAF","GRIN2B","CDK5R1",
               "ENC1","TUBB2A","SYT4","ACTG1","PLPPR4","FAM171B","HS6ST3","LAMP1",
               "OLFM1","JUN","HMP19","ATF4","FTL","STMN2","SCG2","SCG2","CPE","ATP1B1",
               "MAP1A","HSP90AA1","PSAP","NNAT","MPHOSPH8")
litMarker_df <- data.frame(
  Gene = litMarker,
  Number = fData(target_myData)[litMarker, "DetectedSegments"],
  DetectionRate = percent(fData(target_myData)[litMarker, "DetectionRate"]))

#gene filtering workflow
# Plot detection rate:
plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <-
  unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                function(x) {sum(fData(target_myData)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(target_myData))
rownames(plot_detect) <- plot_detect$Freq

ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
            vjust = 1.6, color = "black", size = 4) +
  scale_fill_gradient2(low = "orange2", mid = "lightblue",
                       high = "dodgerblue3", midpoint = 0.65,
                       limits = c(0,1),
                       labels = scales::percent) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "% of Segments",
       y = "Genes Detected, % of Panel > LOQ")
ggsave("~/Directory/NanostringTauProject/plots/detectionrate.pdf")

# Subset to target genes detected in at least 10% of the samples.
#   Also manually include the negative control probe, for downstream use
negativeProbefData <- subset(fData(target_myData), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
target_myData <- 
  target_myData[fData(target_myData)$DetectionRate >= 0.2 |
                  fData(target_myData)$TargetName %in% neg_probes, ]
dim(target_myData)

goi <- goi[goi %in% rownames(target_myData)]
cellMarker <- cellMarker[cellMarker %in% rownames(target_myData)]
litMarker <- litMarker[litMarker %in% rownames(target_myData)]

#normalizing the data----
library(reshape2)  # for melt
library(cowplot)   # for plot_grid
library(tidyverse)

# Graph Q3 value vs negGeoMean of Negatives
ann_of_interest <- "tau"
Stat_data <- 
  data.frame(row.names = colnames(exprs(target_myData)),
             Segment = colnames(exprs(target_myData)),
             Annotation = pData(target_myData)[, ann_of_interest],
             Q3 = unlist(apply(exprs(target_myData), 2,
                               quantile, 0.75, na.rm = TRUE)),
             NegProbe = exprs(target_myData)[neg_probes, ])
Stat_data_m <- melt(Stat_data, measure.vars = c("Q3", "NegProbe"),
                    variable.name = "Statistic", value.name = "Value")

plt1 <- ggplot(Stat_data_m,
               aes(x = Value, fill = Statistic)) +
  geom_histogram(bins = 40) + theme_bw() +
  scale_x_continuous(trans = "log2") +
  facet_wrap(~Annotation, nrow = 1) + 
  scale_fill_brewer(palette = 3, type = "qual") +
  labs(x = "Counts", y = "Segments, #")

plt2 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3, color = Annotation)) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
  geom_point() + guides(color = "none") + theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")

plt3 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)) +
  geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
  geom_point() + theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")

btm_row <- plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""),
                     rel_widths = c(0.43,0.57))
plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))
ggsave("~/Directory/NanostringTauProject/plots/Q3vsnegGeoMean_plotgrid.pdf")

# Q3 norm (75th percentile) for WTA/CTA  with or without custom spike-ins
target_myData <- normalize(target_myData ,
                           norm_method = "quant", 
                           desiredQuantile = .75,
                           toElt = "q_norm")


# visualize the first 10 segments with each normalization method
boxplot(exprs(target_myData)[,1:10],
        col = "#9EDAE5", main = "Raw Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Raw")
pdf(file = "~/Directory/NanostringTauProject/plots/rawsegments.pdf")
dev.off()

boxplot(assayDataElement(target_myData, elt = "q_norm"),
        col = "#2CA02C", main = "Q3 Norm Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Q3 Normalized")
pdf(file = "~/Directory/NanostringTauProject/plots/q3normsegments.pdf")
dev.off()

sampleID <- target_myData$Label
target_myData.mat <- exprs(target_myData)
target_myData.df <- as_tibble(target_myData.mat, rownames = "geneID")
colnames(target_myData.df) <- c("geneID", sampleID)
target_myData.df.pivot <- pivot_longer(target_myData.df, # dataframe to be pivoted
                                       cols = 2:ncol(target_myData.df), # column names to be stored as a SINGLE variable
                                       names_to = "samples", # name of that new variable (column)
                                       values_to = "expression") # name of new variable (column) storing all the values (data)

ggplot(target_myData.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  scale_y_log10() +
  labs(y="Expression", x = "sample",
       title="Raw Data",
       subtitle="gene and segment filtered, Raw data",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()
ggsave("~/Directory/NanostringTauProject/plots/rawexpressiondata_sample.pdf")
dev.off()

target_myDataQ3.mat <- assayDataElement(target_myData, elt = "q_norm")
target_myDataQ3.df <- as_tibble(target_myDataQ3.mat, rownames = "geneID")
colnames(target_myDataQ3.df) <- c("geneID", sampleID)
target_myDataQ3.df.pivot <- pivot_longer(target_myDataQ3.df, # dataframe to be pivoted
                                         cols = 2:ncol(target_myDataQ3.df), # column names to be stored as a SINGLE variable
                                         names_to = "samples", # name of that new variable (column)
                                         values_to = "expression") # name of new variable (column) storing all the values (data)

ggplot(target_myDataQ3.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  scale_y_log10() +
  labs(y="Expression", x = "sample",
       title="Q3 Normalized",
       subtitle="gene and segment filtered, Q3 Normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()
ggsave("~/Directory/NanostringTauProject/plots/q3normexpressiondata_sample.pdf")
dev.off()
#Limma normlization----
library(edgeR) # well known package for differential expression analysis, but we only use for the DGEList object and for normalization methods
library(matrixStats) # let's us easily calculate stats on rows or columns of a data matrix
library(cowplot)

myDGEListFiltered <- DGEList(target_myData.mat)
myDGEListNorm <- calcNormFactors(myDGEListFiltered, method = "TMM")
target_myDataTMM.mat <- cpm(myDGEListNorm, log=FALSE)
target_myDataTMM.df <- as_tibble(target_myDataTMM.mat, rownames = "geneID")
colnames(target_myDataTMM.df) <- c("geneID", sampleID)
target_myDataTMM.df.pivot <- pivot_longer(target_myDataTMM.df, # dataframe to be pivoted
                                          cols = 2:ncol(target_myDataTMM.df), # column names to be stored as a SINGLE variable
                                          names_to = "samples", # name of that new variable (column)
                                          values_to = "expression") # name of new variable (column) storing all the values (data)

ggplot(target_myDataTMM.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  scale_y_log10() +
  labs(y="Expression", x = "sample",
       title="TMM Normalized",
       subtitle="gene and segment filtered, TMM Normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()
ggsave("~/Directory/NanostringTauProject/plots/targettmmnormexpressiondata_sample.pdf")
dev.off()

#Put TMM data into Nanostring object
target_myDataTMM <- target_myData
exprs(target_myDataTMM) <- target_myDataTMM.mat
exprs(target_myData)[1:5, 1:2]
exprs(target_myDataTMM)[1:5, 1:2]


saveRDS(target_myData,"target_myData.rds")

#Plotting the clusters----
library(umap)
library(Rtsne)
target_myDataTMM$group <- factor(target_myDataTMM$group)
target_myDataTMM$sex <- factor(target_myDataTMM$sex)
target_myDataTMM$collection <- factor(target_myDataTMM$collection)

# update defaults for umap to contain a stable random_state (seed)
custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42
# run UMAP
umap_out <-
  umap(t(log2(exprs(target_myDataTMM))),  
       config = custom_umap)
pData(target_myDataTMM)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]
#visualize umap by collection and sex
ggplot(pData(target_myDataTMM),
       aes(x = UMAP1, y = UMAP2, color = collection, shape = sex)) +
  geom_point(size = 3) +
  theme_bw()
ggsave("~/Directory/NanostringTauProject/plots/umapcollectionvssex.pdf")
dev.off()

str(pData(target_myDataTMM))
#visualize umap by tau and sex
ggplot(pData(target_myDataTMM),
       aes(x = UMAP1, y = UMAP2, color = tau, shape = sex)) +
  geom_point(size = 3) +
  theme_bw()
ggsave("~/Directory/NanostringTauProject/plots/umaptauvssex.pdf")
dev.off()

saveRDS(target_myDataTMM,"target_myDataTMM.rds")

#ComBAT correction on collection------
collection_adjusted <- sva::ComBat(dat = exprs(target_myDataTMM), 
                                   batch = target_myDataTMM$collection, par.prior = TRUE)
saveRDS(collection_adjusted,"~/Directory/NanostringTauProject/cache/collection_adjusted.rds")
#load corrected data into nanostring object
exprs(target_myDataTMM) <- collection_adjusted
saveRDS(target_myDataTMM,"~/Directory/NanostringTauProject/cache/target_myDataTMM_collection_adjusted.rds")
