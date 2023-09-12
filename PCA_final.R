##--1/14/2023 Ryan Palaganas

library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(ggplot2)

#read in batch corrected data----
target_myDataTMM<-readRDS("~/Directory/NanostringTauProject/cache/target_myDataTMM_collection_slide_adjusted.rds")

#extract experiment data
data <- exprs(target_myDataTMM)
range(data)

#shifting data by minimum value plus 1 
data <- data + abs(min(data)) + 1
range(data)

# PCA on adjusted data----
comp_merged<- prcomp(data)

#Selecting top two PCAs from calculated PCAs
table <- subset(comp_merged$rotation, select = c(PC1, PC2))

#extract metadata (colDataDf) and bind it to PC1 and PC2 (table)
colDataDf <- pData(target_myDataTMM)
dtable_merged <- cbind(table, colDataDf)

#plot the PCs----
cvar<-"tau"
svar<-"class"
pcaplot <- ggplot(dtable_merged, aes_string(x = "PC1", y = "PC2", color = cvar, shape = svar)) + 
  scale_color_manual(values = c("negative" = "royalblue3", "positive"="red3")) + 
  scale_shape_manual(values=c("AD" = 1, "control" = 17, "PART" = 15)) +
  scale_size_manual(values=c("AD" = 5, "control" = 5, "PART" = 5)) +
  geom_point(size = 5) + labs(color = paste0(cvar), shape = paste0(svar))
pcaplot
ggsave(paste0("~/Directory/NanostringTauProject/plots/",cvar,".",svar,"PCplot.pdf"))

#export top genes----
head(comp_merged$x)
sort(comp_merged$x[,1])[1:10]
sort(comp_merged$x[,1],decreasing = TRUE)[1:10]

x <- data.frame(comp_merged$x[,1:2])
x$genes <- rownames(x)
write_xlsx(x,"~/Directory/NanostringTauProject/cache/pcaloadings.xlsx")


