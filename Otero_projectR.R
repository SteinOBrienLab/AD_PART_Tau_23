##-- 1/14/2023 Ryan Palaganas
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(CoGAPS)
library(SingleCellExperiment)
library(ggplot2)
library(FSA)
library(projectR)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)

target_myDataTMM<-readRDS("~/Directory/NanostringTauProject/cache/target_myDataTMM_collection_slide_adjusted.rds")

#Data from Cobos (Otero-Garcia et al) "Molecular signatures underlying neurofibrillary tangle susceptibility in Alzheimer's disease----
##Projecting patterns into excitatory and inhibitory neuron data 
###Loading rds objects into environment 
excitatory <- readRDS("~/Directory/Raw Data /Excitatory.rds")
inhibitory <- readRDS("~/Directory/Raw Data /Inhibitory.rds")
result <- readRDS("~/Directory/NanostringTauProject/cache/CoGAPSnp5.rds")
cogapsRes<-list("np5"=result)


##Extracting assay data (assay data is already log normalized)
excitatory_edat <- GetAssayData(excitatory, assay = "RNA")
range(excitatory_edat)


excitatory_featurenames <- excitatory@assays$RNA@meta.features$feature_name %>% as.character()
excitatory_featurenames <- as.character(excitatory_featurenames)
rownames(excitatory_edat) <- excitatory_featurenames
str(excitatory_edat)
#data is too large to project into, presubset 33178 features in excitatory dataset by the 1393 features the patterns are learned in
features <- rownames(cogapsRes$np5@featureLoadings)   # Specify features to keep
#index excitatory_edat by features
indx<-which(rownames(excitatory_edat) %in% features)
length(indx)
#restrict features to index
excitatory_edat<-excitatory_edat[indx,]
str(excitatory_edat)
excitatory_data <- data.frame(excitatory_edat)
subset(excitatory_edat, rownames(excitatory_edat) %in% features)
dim(excitatory_edat)
excitatory_data <- as.matrix(excitatory_edat[rownames(excitatory_edat) %in% features, ])  # Extract rows from data
data_subset                                                # Print data frame subset

#creating average expression plots
sampleID <- colnames(excitatory_data)
excitatory_data.df <- as_tibble(excitatory_data, rownames = "geneID")
colnames(excitatory_data.df) <- c("geneID", sampleID)
library(tidyr)
excitatory_data.df.pivot <- pivot_longer(excitatory_data.df, # dataframe to be pivoted
                                          cols = 2:ncol(excitatory_data.df), # column names to be stored as a SINGLE variable
                                          names_to = "samples", # name of that new variable (column)
                                          values_to = "expression") # name of new variable (column) storing all the values (data)

#averaging the data for same group by tau status and slide----

excitatory_data.wide <- pivot_wider(
  excitatory_data.df.pivot, 
  id_cols = samples, 
  names_from = geneID, 
  values_from = expression
)

Sort <- excitatory@meta.data$SORT
Tau <- excitatory_pdata$Tau
Cell.Types <- excitatory_pdata$Cell.Type
sampleid <- excitatory_pdata$Sample.ID

excitatory_data.wide <- excitatory_data.wide %>%
  mutate(Sort = Sort, Tau = Tau, Cell.Types = Cell.Types, sampleid = sampleid)

excitatory_data.long <- excitatory_data.wide %>%
  pivot_longer(
    !(samples|Sort|Tau|Cell.Types|sampleid),
    names_to = "geneID",
    values_to = "expression"
  )

excitatory_data.longAVG <- excitatory_data.long %>%
  group_by(Sort, geneID, Tau, Cell.Types, sampleid) %>%
  summarize(meanExprs=mean(expression))

excitatory_data.wideAVG <- excitatory_data.longAVG %>%
  pivot_wider(
    names_from = geneID, 
    values_from = meanExprs)

#removing cell type groups that were not analyzed in otero et. al.
excitatory_data.wideAVG <- subset(excitatory_data.wideAVG, 
                               Cell.Types!="Ex09_FEZF2-ADRA1A (L5b)" & 
                                 Cell.Types!="Ex12_FEZF2-SYT6 (L6)" & 
                                 Cell.Types!="Ex11_THEMIS-NTNG2 (L6)")
excitatory_data.wideAVG$Cell.Types <- as.character(excitatory_data.wideAVG$Cell.Types)
excitatory_data.wideAVG[excitatory_data.wideAVG == "Ex04_RORB-GABRG1 (L4-L5)"] <- "Ex04_05_RORB-GABRG1-ADGRLA (L4-L5)"
excitatory_data.wideAVG[excitatory_data.wideAVG == "Ex05_RORB-ADGRLA (L5)"] <- "Ex04_05_RORB-GABRG1-ADGRLA (L4-L5)"

SYT1 <- ggplot(excitatory_data.wideAVG) +
  aes(x=Sort, y=log2(SYT1), fill=Tau) +
  geom_boxplot() +
  scale_fill_manual(values = c('Negative' = 'royalblue3', 'Positive' = 'red3')) +
  labs(y = "SYT1 Expression (log2 CPM)") +
  scale_y_continuous() +
  theme_bw() + 
  scale_x_discrete(limits = c("MAP2control", "MAP2", "AT8")) 
SYT1
ggsave(paste0("~/Directory/NanostringTauProject/plots/SYT1_otero.pdf"))

CALM3 <- ggplot(excitatory_data.wideAVG) +
  aes(x=Sort, y=log2(CALM3), fill=Tau) +
  geom_boxplot() +
  scale_fill_manual(values = c('Negative' = 'royalblue3', 'Positive' = 'red3')) +
  labs(y = "CALM3 Expression (log2 CPM)") +
  scale_y_continuous() +
  theme_bw() + 
  scale_x_discrete(limits = c("MAP2control", "MAP2", "AT8")) 
CALM3
ggsave(paste0("~/Directory/NanostringTauProject/plots/CALM3_otero.pdf"))

PRNP <- ggplot(excitatory_data.wideAVG) +
  aes(x=Sort, y=log2(PRNP), fill=Tau) +
  geom_boxplot() +
  scale_fill_manual(values = c('Negative' = 'royalblue3', 'Positive' = 'red3')) +
  labs(y = "PRNP Expression (log2 CPM)") +
  scale_y_continuous() +
  theme_bw() + 
  scale_x_discrete(limits = c("MAP2control", "MAP2", "AT8")) 
PRNP
ggsave(paste0("~/Directory/NanostringTauProject/plots/PRNP_otero.pdf"))

NRGN <- ggplot(excitatory_data.wideAVG) +
  aes(x=Sort, y=log2(NRGN), fill=Tau) +
  geom_boxplot() +
  scale_fill_manual(values = c('Negative' = 'royalblue3', 'Positive' = 'red3')) +
  labs(y = "NRGN Expression (log2 CPM)") +
  scale_y_continuous() +
  theme_bw() + 
  scale_x_discrete(limits = c("MAP2control", "MAP2", "AT8")) 
NRGN
ggsave(paste0("~/Directory/NanostringTauProject/plots/NRGN_otero.pdf"))

SYN3 <- ggplot(excitatory_data.wideAVG) +
  aes(x=Sort, y=log2(SYN3), fill=Tau) +
  geom_boxplot() +
  scale_fill_manual(values = c('Negative' = 'royalblue3', 'Positive' = 'red3')) +
  labs(y = "SYN3 Expression (log2 CPM)") +
  scale_y_continuous() +
  theme_bw() + 
  scale_x_discrete(limits = c("MAP2control", "MAP2", "AT8")) 
SYN3
ggsave(paste0("~/Directory/NanostringTauProject/plots/SYN3_otero.pdf"))

APP <- ggplot(excitatory_data.wideAVG) +
  aes(x=Sort, y=log2(APP), fill=Tau) +
  geom_boxplot() +
  scale_fill_manual(values = c('Negative' = 'royalblue3', 'Positive' = 'red3')) +
  labs(y = "APP Expression (log2 CPM)") +
  scale_y_continuous() +
  theme_bw() + 
  scale_x_discrete(limits = c("MAP2control", "MAP2", "AT8")) 
APP
ggsave(paste0("~/Directory/NanostringTauProject/plots/APP_otero.pdf"))
library(stringr)
#extract meta data
excitatory_pdata <- excitatory@meta.data
excitatory_pdata <- excitatory_pdata %>% mutate(Tau = NA)
excitatory_pdata$Tau[str_detect(excitatory_pdata$Sample.ID, "MAP2")] <- "Negative"
excitatory_pdata$Tau[str_detect(excitatory_pdata$Sample.ID, "AT8")] <- "Positive"
excitatory_pdata$Tau[str_detect(excitatory_pdata$Sample.ID, "CTRL")] <- "Negative"
#bind metadata to experiment data
excitatory_colDataDf <- cbind((p_excitatory_data), excitatory_pdata)
#exclude cell.type groups 9, 12 and 11, merge groups 4 and 5 as in paper
excitatory_colDataDf <- subset(excitatory_colDataDf, 
                                 Cell.Types!="Ex09_FEZF2-ADRA1A (L5b)" & 
                                 Cell.Types!="Ex12_FEZF2-SYT6 (L6)" & 
                                 Cell.Types!="Ex11_THEMIS-NTNG2 (L6)")
excitatory_colDataDf$Cell.Types <- as.character(excitatory_colDataDf$Cell.Types)
excitatory_colDataDf[excitatory_colDataDf == "Ex04_RORB-GABRG1 (L4-L5)"] <- "Ex04_05_RORB-GABRG1-ADGRLA (L4-L5)"
excitatory_colDataDf[excitatory_colDataDf == "Ex05_RORB-ADGRLA (L5)"] <- "Ex04_05_RORB-GABRG1-ADGRLA (L4-L5)"
#reformatting cell type metadata
excitatory_colDataDf <- excitatory_colDataDf %>% mutate(groups = NA)
excitatory_colDataDf$groups[str_detect(excitatory_colDataDf$Cell.Types, "Ex01_")] <- "Ex01 (L2L3)"
excitatory_colDataDf$groups[str_detect(excitatory_colDataDf$Cell.Types, "Ex02_")] <- "Ex02 (L2L4)"
excitatory_colDataDf$groups[str_detect(excitatory_colDataDf$Cell.Types, "Ex03_")] <- "Ex03 (L4L5)"
excitatory_colDataDf$groups[str_detect(excitatory_colDataDf$Cell.Types, "Ex04_05_")] <- "Ex04/05 (L4L5)"
excitatory_colDataDf$groups[str_detect(excitatory_colDataDf$Cell.Types, "Ex06_")] <- "Ex06 (L5)"
excitatory_colDataDf$groups[str_detect(excitatory_colDataDf$Cell.Types, "Ex07_")] <- "Ex07 (L5)"
excitatory_colDataDf$groups[str_detect(excitatory_colDataDf$Cell.Types, "Ex08_")] <- "Ex08 (L5)"
excitatory_colDataDf$groups[str_detect(excitatory_colDataDf$Cell.Types, "Ex10_")] <- "Ex10 (L5L6)"
excitatory_colDataDf$groups[str_detect(excitatory_colDataDf$Cell.Types, "Ex13_")] <- "Ex13 (L6b)"
saveRDS(excitatory_colDataDf, '/Users/ryanpalaganas/Desktop/Directory/NanostringTauProject/cache/excitatory_colDataDf.rds')
excitatory_colDataDf <- readRDS('/Users/ryanpalaganas/Desktop/Directory/NanostringTauProject/cache/excitatory_colDataDf.rds')
saveRDS(test, '/Users/ryanpalaganas/Desktop/Directory/NanostringTauProject/cache/test.rds')
#projecting patterns n=5----
p_excitatory_data <-t(projectR(excitatory_data,cogapsRes$np5@featureLoadings))
saveRDS(p_excitatory_data, '/Users/ryanpalaganas/Desktop/Directory/NanostringTauProject/cache/oteroprojectr.rds')

#exporting genes in each pattern to csv 
gene_patterns <- as.data.frame(cogapsRes$np5@featureLoadings)
write.xlsx(gene_patterns, file = '~/Directory/NanostringTauProject/plots/gene_patterns.xlsx',
           colNames=TRUE, rowNames=TRUE)

library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

#visualize 5 patterns in excitatory neurons----
for(i in 1:length(cogapsRes)) {
  # make metadata table
  pdf(paste0("/Users/ryanpalaganas/Desktop/Directory/NanostringTauProject/plots/",names(cogapsRes[i]), "_boxplots_excitatorytest",".pdf"))
  names(cogapsRes[i])
  for(k in colnames(p_excitatory_data)){
    boxplot_list <- list()
    categoricalVariables <- c("Tau")
    variableTitles <- c("Tau")
    for(l in 1:length(categoricalVariables)){
      boxplot_list[[categoricalVariables[l]]] <- ggplot(excitatory_colDataDf, aes_string(fill = paste(categoricalVariables[l]), x = 'SORT', y = paste(k))) +
        geom_boxplot(alpha = 0.80) + 
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))+
        labs(color = paste(variableTitles[l]), x = paste(variableTitles[l])) +
        #facet_grid(cols = vars(groups), 
                   #scales = "free_x", space = "free_x",
                   #labeller = labeller(groups = label_wrap_gen(width = 6)))+
        scale_fill_manual(values = c('Negative' = 'royalblue3', 'Positive' = 'red3')) +
        scale_x_discrete(limits = c("MAP2control", "MAP2", "AT8")) 
    }
    print(do.call(grid.arrange, boxplot_list))
  }
  dev.off()
}

#students t test pattern 4 positive vs negative 

positive <- filter(excitatory_colDataDf,Tau == "Positive")
positive <- positive$Pattern_4
negative <- filter(excitatory_colDataDf,Tau == "Negative")
negative <- negative$Pattern_4

t.test(positive, negative, var.equal = TRUE) #t = 70.419, df = 91949, p-value < 2.2e-16

#students t test pattern 5 positive vs negative 

positive <- filter(excitatory_colDataDf,Tau == "Positive")
positive <- positive$Pattern_5
negative <- filter(excitatory_colDataDf,Tau == "Negative")
negative <- negative$Pattern_5

t.test(positive, negative, var.equal = TRUE) #t = -36.641, df = 91949, p-value < 2.2e-16

#t.test of every combination of excitatory cell type 
# Transform the data into long format
# Put all variables in the same column except `Cell.type`, the grouping variable
t.data <- excitatory_colDataDf[,c(4,20,40)]
library(tidyverse)
library(rstatix)
library(ggpubr)

#average pattern matrix by sample ID
p_excitatoryDf <- excitatory_colDataDf[,c("Pattern_1", "Pattern_2", "Pattern_3", "Pattern_4", "Pattern_5")]
p_excitatoryDf <- p_excitatoryDf %>% mutate(uniqueid = rownames(excitatory_colDataDf)) %>% select(uniqueid, everything())
p_excitatoryDf.pivot <- pivot_longer(p_excitatoryDf, # dataframe to be pivoted
                                         cols = 2:ncol(p_excitatoryDf), # column names to be stored as a SINGLE variable
                                         names_to = "patterns", # name of that new variable (column)
                                         values_to = "expression")
p_excitatoryDf.wide <- pivot_wider(
  p_excitatoryDf.pivot, 
  id_cols = uniqueid, 
  names_from = patterns, 
  values_from = expression
) 


Tau <- excitatory_colDataDf$Tau
Cell.Types <- excitatory_colDataDf$Cell.Type
sampleid <- excitatory_colDataDf$Sample.ID

p_excitatoryDf.wide <- p_excitatoryDf.wide %>%
  mutate(Tau = Tau, Cell.Types = Cell.Types, sampleid = sampleid)


p_excitatoryDf.long <- p_excitatoryDf.wide %>%
  pivot_longer(
    !(uniqueid|Tau|Cell.Types|sampleid),
    names_to = "patternID",
    values_to = "expression"
  )

p_excitatoryDf.longAVG <- p_excitatoryDf.long %>%
  group_by(patternID, Tau, sampleid, Cell.Types) %>%
  summarize(meanExprs=mean(expression))

p_excitatoryDf.wideAVG <- p_excitatoryDf.longAVG %>%
  pivot_wider(
    names_from = patternID, 
    values_from = meanExprs)

test <- excitatory_colDataDf %>% group_by(donor_id, Cell.Types, Tau) %>% 
  summarize(avgp1 = mean(Pattern_1),
            avgp2 = mean(Pattern_2),
            avgp3 = mean(Pattern_3),
            avgp4 = mean(Pattern_4),
            avgp5 = mean(Pattern_5))

stat.test.p4 <- test %>%
  group_by(Cell.Types) %>%
  t_test(avgp4 ~ Tau) %>%
  add_significance()
write.xlsx(stat.test.p4, file = '/Users/ryanpalaganas/Desktop/Directory/NanostringTauProject/plots/stat.test.p4noadj.xlsx',
           colNames=TRUE, rowNames=TRUE)

stat.test.p5 <- test %>%
  group_by(Cell.Types) %>%
  t_test(avgp5 ~ Tau) %>%
  add_significance()
write.xlsx(stat.test.p5, file = '/Users/ryanpalaganas/Desktop/Directory/NanostringTauProject/plots/stat.test.p5.xlsx',
           colNames=TRUE, rowNames=TRUE)




#anova
two.way <- aov(Pattern_4 ~ Tau * Cell.Types, data = excitatory_colDataDf)
summary(two.way)
two.way <- aov(Pattern_5 ~ Tau * Cell.Types, data = excitatory_colDataDf)
summary(two.way)

#displaying gene lists that contribute to patterns below 
Pattern3 <- cogapsRes[["np5"]]@featureLoadings[,3]
Pattern4 <- cogapsRes[["np5"]]@featureLoadings[,4]
Pattern5 <- cogapsRes[["np5"]]@featureLoadings[,5]
head(sort(Pattern5, decreasing = TRUE), n=20)

#visualize 10 patterns in excitatory neurons
for(i in 1:length(cogapsRes10)) {
  # extract pattern matrix
  p_excitatory_data_10 <- t(p_excitatory_data_10)
  # make metadata table
  p_excitatory_colDataDf_10 <- cbind(p_excitatory_data_10, excitatory_pdata)
  pdf(paste0("/Users/ryanpalaganas/Desktop/Directory/NanostringTauProject/plots/",names(cogapsRes10[i]), "_boxplots_10pattern.pdf"))
  names(cogapsRes10[i])
  for(k in colnames(p_excitatory_data_10)){
    boxplot_list <- list()
    # extract the categorical variables
    categoricalVariables <- c("Tau")
    variableTitles <- c("Tau")
    for(l in 1:length(categoricalVariables)){
      boxplot_list[[categoricalVariables[l]]] <- ggplot(p_excitatory_colDataDf_10, aes_string(color = paste(categoricalVariables[l]), x = paste(categoricalVariables[l]), y = paste(k))) +
        geom_boxplot() + 
        geom_jitter(shape=16, position=position_jitter(0.2)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))+
        labs(color = paste(variableTitles[l]), x = paste(variableTitles[l])) #+
        #facet_grid(rows = "subclusterAssignment")
    }
    print(do.call(grid.arrange, boxplot_list))
  }
  dev.off()
}

#visualize 10 patterns in excitatory neurons with stdloadings
for(i in 1:length(cogapsRes10)) {
  # extract pattern matrix
  p_excitatory_data_10_loadings <- t(p_excitatory_data_10_loadings)
  # make metadata table
  p_excitatory_colDataDf_10_l <- cbind(p_excitatory_data_10_loadings, excitatory_pdata)
  pdf(paste0("/Users/ryanpalaganas/Desktop/Directory/NanostringTauProject/plots/",names(cogapsRes10[i]), "_boxplots_10pattern_loadings.pdf"))
  names(cogapsRes10[i])
  for(k in colnames(p_excitatory_data_10_loadings)){
    boxplot_list <- list()
    # extract the categorical variables
    categoricalVariables <- c("Tau")
    variableTitles <- c("Tau")
    for(l in 1:length(categoricalVariables)){
      boxplot_list[[categoricalVariables[l]]] <- ggplot(p_excitatory_colDataDf_10_l, aes_string(color = paste(categoricalVariables[l]), x = paste(categoricalVariables[l]), y = paste(k))) +
        geom_boxplot() + 
        geom_jitter(shape=16, position=position_jitter(0.2)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))+
        labs(color = paste(variableTitles[l]), x = paste(variableTitles[l])) #+
      #facet_grid(rows = "subclusterAssignment")
    }
    print(do.call(grid.arrange, boxplot_list))
  }
  dev.off()
}





##i think this can be deleted
library("readxl")
genes <- read_excel("/Users/ryanpalaganas/Downloads/SynGoAnalysisOverall.xlsx")
genes <- genes$`gene symbol`
str(genes)
genes <- genes[0:353]

measures <- list("genes"=genes)
genes

cogapsRes<-unlist(cogapsRes)
cogapsRes

calcCoGAPSStat(
  cogapsRes,
  sets = measures,
  whichMatrix = "featureLoadings",
  numPerm = 1000
  #GStoGenes = genes
)

class(cogapsRes)<-"CogapsResult"

