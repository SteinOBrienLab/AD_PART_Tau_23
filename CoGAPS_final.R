##-- 1/14/2023 Ryan Palaganas
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(CoGAPS)
library(SingleCellExperiment)
library(ggplot2)
library(gridExtra)
#read in preprocessed, batch (collection) corrected data----
target_myDataTMM<-readRDS("/Users/ryanpalaganas/Desktop/Directory/NanostringTauProject/cache/target_myDataTMM_collection_adjusted.rds")

#Applying batch correction for slide (slide and collection)----
slide_adjusted <- sva::ComBat(dat = exprs(target_myDataTMM), 
                              batch =target_myDataTMM$slide_name, par.prior = TRUE)


#load batch corrected data into nanostring object
exprs(target_myDataTMM) <- slide_adjusted
saveRDS(target_myDataTMM,"(~/Directory/NanostringTauProject/cache/target_myDataTMM_collection_slide_adjusted.rds")


#extract experiment data----
data <- exprs(target_myDataTMM)
range(data)

#shifting data by minimum value plus 1 
data <- data + abs(min(data)) + 1
range(data)

#learn patterns and create CoGAPS object----
params <- new("CogapsParams")
params <- setParam(params, "nPatterns", 5)
params <- setParam(params, "seed",834)
params <- setParam(params,"geneNames", rownames(data))

#create single cell experiment object
sce<-SingleCellExperiment(list(counts=data))

#create CoGAPS object
result <-CoGAPS(sce, params, outputFrequency=500)
saveRDS(result,"~/Directory/NanostringTauProject/cache/CoGAPSnp5.rds")
result <- readRDS("/Users/ryanpalaganas/Desktop/Directory/NanostringTauProject/cache/CoGAPSnp5.rds")
cogapsRes<-list("np5"=result)

str(cogapsRes)

#pattern visualization----
for(i in 1:length(cogapsRes)) {
  # extract pattern matrix
  tempDf <- as.data.frame(cogapsRes[[i]]@sampleFactors)
  # make metadata table
  colDataDf <- cbind(tempDf, pData(target_myDataTMM))
  pdf(paste0("/Users/ryanpalaganas/Desktop/Directory/NanostringTauProject/plots/",names(cogapsRes[i]), "_boxplotst.pdf"))
  names(cogapsRes[i])
  for(k in colnames(tempDf)){
    boxplot_list <- list()
    # extract the categorical variables
    categoricalVariables <- c("tau","collection","class","slide_name","group")
    variableTitles <- c("tau","collection","class","slide_name","group")
    for(l in 1:length(categoricalVariables)){
      boxplot_list[[categoricalVariables[l]]] <- ggplot(colDataDf, aes_string(color = paste(categoricalVariables[l]), x = paste(categoricalVariables[l]), y = paste(k))) +
        geom_boxplot() + 
        geom_jitter(shape=16, position=position_jitter(0.2)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))+
        labs(color = paste(variableTitles[l]), x = paste(variableTitles[l])) +
        facet_grid()
    }
    print(do.call(grid.arrange, boxplot_list))
  }
  dev.off()
}

#visualize patterns by disease and tau
for(i in 1:length(cogapsRes)) {
  # extract pattern matrix
  tempDf <- as.data.frame(cogapsRes[[i]]@sampleFactors)
  # make metadata table
  colDataDf <- cbind(tempDf, pData(target_myDataTMM))
  pdf(paste0("/Users/ryanpalaganas/Desktop/Directory/NanostringTauProject/plots/",names(cogapsRes[i]), "_boxplotst.pdf"))
  names(cogapsRes[i])
  for(k in colnames(tempDf)){
    boxplot_list <- list()
    # extract the categorical variables
    categoricalVariables <- c("tau")
    variableTitles <- c("tau")
    for(l in 1:length(categoricalVariables)){
      boxplot_list[[categoricalVariables[l]]] <- ggplot(colDataDf, aes_string(fill = paste(categoricalVariables[l]), x = 'class', y = paste(k))) +
        geom_boxplot(alpha = 0.80, outlier.shape = NA) + 
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))+
        labs(color = paste(variableTitles[l]), x = paste(variableTitles[l])) +
        scale_fill_manual(values = c('negative' = 'royalblue3', 'positive' = 'red3')) +
        scale_x_discrete(limits = c("control", "PART", "AD")) 
    }
    print(do.call(grid.arrange, boxplot_list))
  }
  dev.off()
}

#students t test pattern 4 positive vs negative----
tempDf <- as.data.frame(cogapsRes[[1]]@sampleFactors)
t_colDataDf <- cbind(tempDf, pData(target_myDataTMM))
t_colDataDf <- as.factor(t_colDataDf$class)

#create a data obect consisting of only tau positive and average by slide name
positive <- filter(t_colDataDf,tau == "positive") %>% as.data.table()

p <- positive[, mean(Pattern_4), by = slide_name]

#create a data object consisting of only tau negative and average by slide name
negative <- filter(t_colDataDf,tau == "negative") %>% as.data.table()

n <- negative[, mean(Pattern_4), by = slide_name]
#perform students t-test
t.test(p$V1, n$V1, var.equal = TRUE)

#create a data obect consisting of only tau positive and average by slide name
positive <- filter(t_colDataDf,tau == "positive") %>% as.data.table()

p <- positive[, mean(Pattern_5), by = slide_name]

#create a data object consisting of only tau negative and average by slide name
negative <- filter(t_colDataDf,tau == "negative") %>% as.data.table()

n <- negative[, mean(Pattern_5), by = slide_name]
#perform students t-test
t.test(p$V1, n$V1, var.equal = TRUE)

#anova
mway <- aov(Pattern_4 ~ tau + sex + class, data = t_colDataDf)
summary(mway)
mway <- aov(Pattern_5 ~ tau + sex + class, data = t_colDataDf)
summary(mway)
#disease vs control anova: i am aggregating both PART and AD into 'disease'
t_colDataDf <- t_colDataDf %>% 
  mutate(disease = case_when(
    class == "control" ~ "control",
    class == "AD" ~ "disease",
    class == "PART" ~ "disease"
  ))
t_colDataDf$disease <- as.factor(t_colDataDf$disease)
t_colDataDf$tau <- as.factor(t_colDataDf$tau)
t_colDataDf <- t_colDataDf %>% 
  mutate(control = case_when(
    class == "control" ~ 1,
    class == "AD" ~ 0,
    class == "PART" ~ 0
  ))
t_colDataDf <- t_colDataDf %>% 
  mutate(AD = case_when(
    class == "control" ~ 0,
    class == "AD" ~ 1,
    class == "PART" ~ 0
  ))
t_colDataDf <- t_colDataDf %>% 
  mutate(PART = case_when(
    class == "control" ~ 0,
    class == "AD" ~ 0,
    class == "PART" ~ 1
  ))


anvtst <- t_colDataDf %>% anova_test(Pattern_4 ~ tau*AD*PART + sex)

#4 t-tests: [tau+/AD vs tau+/PART], [tau-/AD vs tau-/PART], [tau+/AD vs tau-/AD], [tau+/PART vs tau-/PART]
taupos.AD <- filter(t_colDataDf, tau == "positive"  & class=="AD") %>% as.data.table() #tau+/AD
taupos.PART <- filter(t_colDataDf, tau == "positive"  & class=="PART") %>% as.data.table() #tau+/PART
tauneg.AD <- filter(t_colDataDf, tau == "negative"  & class=="AD") %>% as.data.table() #tau-/AD
tauneg.PART <- filter(t_colDataDf, tau == "negative"  & class=="PART") %>% as.data.table() #tau-/PART
statdf <- data.frame()
#t.test [tau+/AD vs tau+/PART]
#pattern4
t.test(taupos.AD$Pattern_4, taupos.PART$Pattern_4, var.equal = TRUE)
#pattern5
t.test(taupos.AD$Pattern_5, taupos.PART$Pattern_5, var.equal = TRUE)

#t.test [tau-/AD vs tau-/PART]
#pattern 4
t.test(tauneg.AD$Pattern_4, tauneg.PART$Pattern_4, var.equal = TRUE)
#pattern 5
t.test(tauneg.AD$Pattern_5, tauneg.PART$Pattern_5, var.equal = TRUE)

#t.test [tau+/AD vs tau-/AD]
t.test(tposAD$V1, tnegAD$V1, var.equal = TRUE)
#pattern 4
t.test(taupos.AD$Pattern_4, tauneg.AD$Pattern_4, var.equal = TRUE)
#pattern 5
t.test(taupos.AD$Pattern_5, tauneg.AD$Pattern_5, var.equal = TRUE)

#t.test [tau+/PART vs tau-/PART]
x <- t.test(tposPART$V1, tnegPART$V1, var.equal = TRUE)
#pattern 4 
t.test(taupos.PART$Pattern_4, tauneg.PART$Pattern_4, var.equal = TRUE)
#pattern 5 
t.test(taupos.PART$Pattern_5, tauneg.PART$Pattern_5, var.equal = TRUE)


mway <- aov(Pattern_4 ~ tau*disease + sex, data = t_colDataDf)
summary(mway)


library(knitr)
library(rstatix)
summary(mway)

res.aov <- t_colDataDf %>% anova_test(Pattern_4 ~ tau + slide_name + sex + class)
kable(summary(res.aov))
#learning 10 patterns----
params10 <- new("CogapsParams")
params10 <- setParam(params10, "nPatterns", 10)
params10 <- setParam(params10,"geneNames", rownames(data))

sce10<-SingleCellExperiment(list(counts=data))

exprs(sce10)

result10 <-CoGAPS(sce10, params10, outputFrequency=500) #,  geneNames = )
saveRDS(result10,"~/Directory/NanostringTauProject/cache/CoGAPSnp10.rds")

cogapsRes10<-list("np10"=result10)
str(cogapsRes10)

#visualize 10 patterns----
for(i in 1:length(cogapsRes10)) {
  # extract pattern matrix
  tempDf <- as.data.frame(cogapsRes10[[i]]@sampleFactors)
  # make metadata table
  colDataDf <- cbind(tempDf, pData(target_myDataTMM))
  pdf(paste0("~/Directory/NanostringTauProject/plots/",names(cogapsRes10[i]), "_boxplots.pdf"))
  names(cogapsRes10[i])
  for(k in colnames(tempDf)){
    boxplot_list10 <- list()
    # extract the categorical variables
    categoricalVariables <- c("tau", "sex","collection","class","slide_name","group")
    variableTitles <- c("tau", "sex","collection","class","slide_name","group")
    for(l in 1:length(categoricalVariables)){
      boxplot_list10[[categoricalVariables[l]]] <- ggplot(colDataDf, aes_string(color = paste(categoricalVariables[l]), x = paste(categoricalVariables[l]), y = paste(k))) +
        geom_boxplot() + 
        geom_jitter(shape=16, position=position_jitter(0.2)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))+
        labs(color = paste(variableTitles[l]), x = paste(variableTitles[l])) 
    }
    print(do.call(grid.arrange, boxplot_list10))
  }
  dev.off()
}

##heatmaps----
#read in CoGAPSnp5 object
result <- readRDS("~/Directory/NanostringTauProject/cache/CoGAPSnp5.rds")

cogapsRes<-list("np5"=result)

#extract patternMarkers----
pm <- patternMarkers(cogapsRes$np5, threshold = "all")

saveRDS(pm,"~/Directory/NanostringTauProject/cache/patternmarkers.rds")
#export patternmarkers, patternmarker ranks and patternmarkerscores as csv
library(stringi)
patternmarker <- list(patternmarkers = data.frame(stri_list2matrix(pm$PatternMarkers, byrow=FALSE)), 
                      patternmarkerranks = data.frame(pm$PatternMarkerRanks, 'genes' = rownames(pm[["PatternMarkerRanks"]])),
                      patternmarkerscores = data.frame(pm$PatternMarkerScores, 'genes' = rownames(pm[["PatternMarkerRanks"]])))
colnames(patternmarker$patternmarkers) = c("Pattern1","Pattern2","Pattern3","Pattern4","Pattern5")

write_xlsx(patternmarker,"~/Directory/NanostringTauProject/cache/patternmarkers.xlsx")

library('ComplexHeatmap')
library(circlize)
#extract metadata and drop control group
colDataDf <- pData(target_myDataTMM)%>% subset(class!='control')
names <- rownames(colDataDf)

#extract top 20 pattern marker genes from patterns 4 and 5 (pattern 4 is tau up, pattern 5 is tau down)
GOI4 <- pm$PatternMarkers$Pattern_4[1:21]
#extract top 21 genes for pattern 5 because RUFY3 will be excluded as an outlier
GOI5 <- pm$PatternMarkers$Pattern_5[1:22]

#subset expression data by top 20 pattern 4 and 5 marker genes, restrict to experimental groups
p4hm <- subset(exprs(target_myDataTMM), rownames(exprs(target_myDataTMM)) %in% GOI4, colnames(exprs(target_myDataTMM)) %in% names)
p5hm <- subset(exprs(target_myDataTMM), rownames(exprs(target_myDataTMM)) %in% GOI5, colnames(exprs(target_myDataTMM)) %in% names)

#row mean normalization function----
rowmeannormalize <- function(x){
  (x-mean(x))/(max(x)-mean(x))
}
#create normalized data object using normalization function for tau up and tau down patterns
data_norm4 <- t(apply(p4hm, 1, rowmeannormalize))
data_norm5 <- t(apply(p5hm, 1, rowmeannormalize))
#setting up complex heatmap variables----
Tau = sample(paste0("Tau", 1:21), nrow(data_norm4), replace = TRUE)
Tau_level = intersect(paste0("Tau", 1:21), Tau)
column_ha = HeatmapAnnotation(Tau = colDataDf$tau,
                              Class = colDataDf$class,
                              col = list(Tau = c("negative" = "royalblue3", "positive" = "red3"),
                                         Class = c("AD" = "red3", "PART" = "yellow2")),
                              annotation_legend_param = list(Tau = list()))

#color spectrum for heatmap
col_fun = colorRamp2(c(-1.5, 0, 1.5), c("purple", "white", "orange"))
col_fun(seq(-3, 3))

#create and export heatmaps----
#pattern tau up
pdf("~/Directory/NanostringTauProject/plots/pattern4t20.pdf")
Heatmap(data_norm4, 
        top_annotation = column_ha, 
        show_row_names = TRUE, 
        show_column_names = FALSE, 
        col = col_fun)
dev.off()
#tau down, RUFY3 was excluded as an outlier
data_norm5 <- data_norm5[!row.names(data_norm5) %in% 'RUFY3',]
#reordering top annotation dendrogram for pattern5 heatmap
library(dendsort)
library(dendextend)
#column dendrogram re-ordering for tau
col_dend5 = as.dendrogram(hclust(dist(t(data_norm5)))) %>% sort()
pdf("~/Directory/NanostringTauProject/plots/pattern5t20.pdf")
Heatmap(data_norm5, 
        name = 'mat', 
        top_annotation = column_ha, 
        show_row_names = TRUE, 
        show_column_names = FALSE, 
        col = col_fun,
        cluster_columns = col_dend5)
dev.off()
