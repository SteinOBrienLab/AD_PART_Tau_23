## 1/14/2023 Ryan Palaganas
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(ggplot2)
#read in collection corrected data object----
target_myDataTMM <- readRDS("~/Directory/NanostringTauProject/cache/target_myDataTMM_collection_adjusted.rds")
#averaging the data for same group by tau status and slide----

#create a vector of slide names
sampleID <- target_myDataTMM$Label
#create an expression matrix (genes and slide names)
target_myDataTMM.df <- as_tibble(exprs(target_myDataTMM), rownames = "geneID")
colnames(target_myDataTMM.df) <- c("geneID", sampleID)
#convert data to long format with gene, slide name and expression column
target_myDataTMM.df.pivot <- pivot_longer(target_myDataTMM.df, # dataframe to be pivoted
                                          cols = 2:ncol(target_myDataTMM.df), # column names to be stored as a SINGLE variable
                                          names_to = "samples", # name of that new variable (column)
                                          values_to = "expression") # name of new variable (column) storing all the values (data)
#transpose expression matrix (target_myDataTMM.df) to slide names x genes
target_myDataTMM.wide <- pivot_wider(
  target_myDataTMM.df.pivot, 
  id_cols = samples, 
  names_from = geneID, 
  values_from = expression
)

#create metadata character strings for Tau, slide and group. Bind these to the expression matrix
Tau <- target_myDataTMM$tau
Slide <- target_myDataTMM$slide_name
Group <- target_myDataTMM$group

target_myDataTMM.wide <- target_myDataTMM.wide %>%
  mutate(tau=Tau, slide=Slide, group=Group)

#convert data to long format and group expression by slide and sample name, tau status and group
target_myDataTMM.long <- target_myDataTMM.wide %>%
  pivot_longer(
    !(samples|tau|slide|group),
    names_to = "geneID",
    values_to = "expression"
  )

#average gene expression by tau, slide and group
target_myDataTMM.longAVG <- target_myDataTMM.long %>%
  group_by(tau,slide,geneID,group) %>%
  summarize(meanExprs=mean(expression))

#convert to wide format
target_myDataTMM.wideAVG <- target_myDataTMM.longAVG %>%
  pivot_wider(
    names_from = geneID, 
    values_from = meanExprs)
#factor by metadata of interest
target_myDataTMM.wideAVG$tau <- factor(target_myDataTMM.wideAVG$tau)
target_myDataTMM.wideAVG$group <- factor(target_myDataTMM.wideAVG$group)
target_myDataTMM.wideAVG$slide <- factor(target_myDataTMM.wideAVG$slide)
TMMAVGwide.mat <- as.matrix(target_myDataTMM.wideAVG[ ,4:ncol(target_myDataTMM.wideAVG)])

#Plotting the AVG sample clusters----
library(umap)

# update defaults for umap to contain a stable random_state (seed)
custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42
# run UMAP
umap_out2 <-
  umap(TMMAVGwide.mat,  
       config = custom_umap)

umap_resultsTMM <- data.frame(umap_out2$layout[, c(1,2)])

ggplot(umap_resultsTMM,
       aes(x = X1, y = X2, color=target_myDataTMM.wideAVG$group)) +
  geom_point(size = 3) +
  theme_bw()
ggsave("~/Directory/NanostringTauProject/plots/umapresults.pdf")
dev.off()
#create a condition column based on group ID
x <- mutate(target_myDataTMM.wideAVG, condition = case_when(group == 1 ~ 'control',
                                                            group == 2 ~ 'PART',
                                                            group == 3 ~ 'PART',
                                                            group == 4 ~ 'AD',
                                                            group == 5 ~ 'AD'))
#plotting average sample expression data for specific genes by group----

APP <- ggplot(target_myDataTMM.wideAVG) +
  aes(x=group, y=log2(APP), fill=tau) +
  geom_boxplot() +
  scale_fill_manual(values = c('negative' = 'royalblue3', 'positive' = 'red3')) +
  labs(y = "APP Expression (log2 CPM)") +
  scale_y_continuous() +
  theme_bw()
APP
ggsave(paste0("~/Directory/NanostringTauProject/plots/APP.pdf"))

PRNP <- ggplot(target_myDataTMM.wideAVG) +
  aes(x=group, y=log2(PRNP), fill=tau) +
  geom_boxplot() +
  scale_fill_manual(values = c('negative' = 'royalblue3', 'positive' = 'red3')) +
  labs(y = "PRNP Expression (log2 CPM)") +
  scale_y_continuous() +
  theme_bw()
PRNP
ggsave(paste0("~/Directory/NanostringTauProject/plots/PRNP.pdf"))

SYN3 <- ggplot(target_myDataTMM.wideAVG) +
  aes(x=group, y=log2(SYN3), fill=tau) +
  geom_boxplot() +
  scale_fill_manual(values = c('negative' = 'royalblue3', 'positive' = 'red3')) +
  labs(y = "SYN3 Expression (log2 CPM)") +
  scale_y_continuous() +
  theme_bw()
SYN3
ggsave(paste0("~/Directory/NanostringTauProject/plots/SYN3.pdf"))

NRGN <- ggplot(target_myDataTMM.wideAVG) +
  aes(x=group, y=log2(NRGN), fill=tau) +
  geom_boxplot() +
  scale_fill_manual(values = c('negative' = 'royalblue3', 'positive' = 'red3')) +
  labs(y = "NRGN Expression (log2 CPM)") +
  scale_y_continuous() +
  theme_bw()
NRGN
ggsave(paste0("~/Directory/NanostringTauProject/plots/NRGN.pdf"))

CALM3 <- ggplot(target_myDataTMM.wideAVG) +
  aes(x=group, y=log2(CALM3), fill=tau) +
  geom_boxplot() +
  scale_fill_manual(values = c('negative' = 'royalblue3', 'positive' = 'red3')) +
  labs(y = "CALM3 Expression (log2 CPM)") +
  scale_y_continuous() +
  theme_bw()
CALM3
ggsave(paste0("~/Directory/NanostringTauProject/plots/CALM3.pdf"))

SYT1 <- ggplot(target_myDataTMM.wideAVG) +
  aes(x=group, y=log2(SYT1), fill=tau) +
  geom_boxplot() +
  scale_fill_manual(values = c('negative' = 'royalblue3', 'positive' = 'red3')) +
  labs(y = "SYT1 Expression (log2 CPM)") +
  scale_y_continuous() +
  theme_bw()
SYT1
ggsave(paste0("~/Directory/NanostringTauProject/plots/SYT1.pdf"))

#plotting average sample expression data for specific genes by slide----
g1 <- filter(target_myDataTMM.wideAVG, group == 1)
APP <- ggplot(x) +
  aes(x=slide, y=log2(APP), fill=tau) +
  geom_boxplot() +
  scale_fill_manual(values = c('negative' = 'royalblue3', 'positive' = 'red3')) +
  labs(y = "APP Expression (log2 CPM)") +
  scale_y_continuous() +
  theme_bw() +
  facet_wrap(~condition, ncol = 5)
APP <- APP + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste0("~/Directory/NanostringTauProject/plots/APPslide.pdf"))

PRNP <- ggplot(x) +
  aes(x=slide, y=log2(PRNP)) +
  geom_boxplot() +
  scale_fill_manual(values = c('negative' = 'royalblue3', 'positive' = 'red3')) +
  labs(y = "PRNP Expression (log2 CPM)") +
  scale_y_continuous() +
  theme_bw() +
  facet_wrap(~condition, ncol = 5)
PRNP <- PRNP + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste0("~/Directory/NanostringTauProject/plots/PRNPslide.pdf"))

PALM3 <- ggplot(x) +
  aes(x=slide, y=log2(PRNP)) +
  geom_boxplot() +
  scale_fill_manual(values = c('negative' = 'royalblue3', 'positive' = 'red3')) +
  labs(y = "PRNP Expression (log2 CPM)") +
  scale_y_continuous() +
  theme_bw() +
  facet_wrap(~condition, ncol = 5)
PRNP <- PRNP + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste0("~/Directory/NanostringTauProject/plots/PRNPslide.pdf"))

#preparing for mixed effects linear models and graphing volcano plots----
# convert test variables to factors
pData(target_myDataTMM)$testTau <- 
  factor(pData(target_myDataTMM)$tau, c("positive", "negative"))
pData(target_myDataTMM)[["slide"]] <- 
  factor(pData(target_myDataTMM)[["slide_name"]])
assayDataElement(object = target_myDataTMM, elt = "log_q") <-
  assayDataApply(target_myDataTMM, 2, FUN = log, base = 2)
library(ggrepel) 

# run LMM within slide by Tau----
# formula follows conventions defined by the lme4 package
results <- c()
for(status in c("AD", "PART")) {
  ind <- pData(target_myDataTMM)$class == status
  mixedOutmc <-
    mixedModelDE(target_myDataTMM[, ind],
                 elt = "log_q",
                 modelFormula = ~ testTau + (1 + testTau | slide),
                 groupVar = "testTau",
                 nCores = parallel::detectCores(),
                 multiCore = FALSE)
  
  # format results as data.frame
  r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- as.data.frame(r_test)
  r_test$Contrast <- tests
  
  # use lapply in case you have multiple levels of your test factor to
  # correctly associate gene name with it's row in the results table
  r_test$Gene <- 
    unlist(lapply(colnames(mixedOutmc),
                  rep, nrow(mixedOutmc["lsmeans", ][[1]])))
  r_test$Subset <- status
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate", 
                       "Pr(>|t|)", "FDR")]
  results <- rbind(results, r_test)
}
library(knitr)
kable(tail(results), digits = 3,
      caption = "DE results",
      align = "lc", row.names = FALSE)

# Categorize Results based on P-value & FDR for plotting
results$Color <- "NS or FC < 0.2"
results$Color[results$`Pr(>|t|)` < 0.05] <- "P < 0.05"
results$Color[results$FDR < 0.05] <- "FDR < 0.05"
results$Color[results$FDR < 0.001] <- "FDR < 0.001"
results$Color[abs(results$Estimate) < 0.2] <- "NS or FC < 0.2"
results$Color <- factor(results$Color,
                        levels = c("NS or FC < 0.2", "P < 0.05",
                                   "FDR < 0.05", "FDR < 0.001"))

# pick top genes for either side of volcano to label
# order genes for convenience:
results$invert_P <- (-log10(results$`Pr(>|t|)`)) * sign(results$Estimate)
top_g <- c()
for(cond in c("AD", "PART")) {
  ind <- results$Subset == cond
  top_g <- c(top_g,
             results[ind, 'Gene'][
               order(results[ind, 'invert_P'], decreasing = TRUE)[1:15]],
             results[ind, 'Gene'][
               order(results[ind, 'invert_P'], decreasing = FALSE)[1:15]])
}
top_g <- unique(top_g)
results <- results[, -1*ncol(results)] # remove invert_P from matrix
range(results$Estimate)
# Graph results
tauvolcano <- ggplot(results,
                     aes(x = Estimate, y = -log10(`Pr(>|t|)`),
                         color = Color, label = Gene)) +
  geom_vline(xintercept = c(0.2, -0.2), lty = "dashed") +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point() +
  labs(x = "Enriched in Tau Negative <- log2(FC) -> Enriched in Tau Positive",
       y = "Significance, -log10(P)",
       color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
                                `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",
                                `NS or FC < 0.2` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(results, Gene %in% top_g & FDR < 0.05),
                  size = 4, point.padding = 0.15, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom") +
  facet_wrap(~Subset, scales = "free_y")
tauvolcano
ggsave(paste0("~/Directory/NanostringTauProject/plots/tauvolcano.pdf"))

#export results as excel
library(writexl)
write_xlsx(results,"~/Directory/NanostringTauProject/plots/TauComparison.xlsx")

# run LMM across slides for diseases----
# convert test variables to factors
pData(target_myDataTMM)$testClass <-
  factor(pData(target_myDataTMM)$class, c("control", "PART","AD"))

# formula follows conventions defined by the lme4 package
resultsClass <- c()
mixedOutmcClass <-
  mixedModelDE(target_myDataTMM,
               elt = "log_q",
               modelFormula =  ~  testClass + (1 | slide),
               groupVar = "testClass",
               nCores = parallel::detectCores(),
               multiCore = FALSE)
# format results as data.frame
r_test2 <- do.call(rbind, mixedOutmcClass["lsmeans", ])
tests2 <- rownames(r_test2)
r_test2 <- as.data.frame(r_test2)
r_test2$Contrast <- tests2

# use lapply in case you have multiple levels of your test factor to
# correctly associate gene name with it's row in the results table
r_test2$Gene <- 
  unlist(lapply(colnames(mixedOutmcClass),
                rep, nrow(mixedOutmcClass["lsmeans", ][[1]])))
r_test2$FDR <- p.adjust(r_test2$`Pr(>|t|)`, method = "fdr")
r_test2 <- r_test2[, c("Gene", "Contrast", "Estimate", 
                       "Pr(>|t|)", "FDR")]
resultsClass <- rbind(resultsClass, r_test2)
resultsClass$Contrast <- factor(resultsClass$Contrast)

resultsClass$Color <- "NS or FC < 0.2"
resultsClass$Color[resultsClass$`Pr(>|t|)` < 0.05] <- "P < 0.05"
resultsClass$Color[resultsClass$FDR < 0.05] <- "FDR < 0.05"
resultsClass$Color[resultsClass$FDR < 0.001] <- "FDR < 0.001"
resultsClass$Color[abs(resultsClass$Estimate) < 0.2] <- "NS or FC < 0.2"
resultsClass$Color <- factor(resultsClass$Color,
                             levels = c("NS or FC < 0.2", "P < 0.05",
                                        "FDR < 0.05", "FDR < 0.001"))

# pick top genes for either side of volcano to label
# order genes for convenience:
resultsClass$invert_P <- (-log10(resultsClass$`Pr(>|t|)`)) * sign(resultsClass$Estimate)

resultsClass <- resultsClass[, -1*ncol(resultsClass)]# remove invert_P from matrix and multiply estimate by -1
resultsClass$Estimate <- resultsClass$Estimate*-1
# Graph results (multiply 'Estimate' by -1)
classcomparison <- ggplot(resultsClass,
                          aes(x = Estimate, y = -log10(`Pr(>|t|)`),
                              color = Color, label = Gene)) +
  geom_vline(xintercept = c(0.2, -0.2), lty = "dashed") +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point() +
  labs(x = "Enriched in Control <- log2(FC) -> Enriched in Disease",
       y = "Significance, -log10(P)",
       color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
                                `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",
                                `NS or FC < 0.2` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(resultsClass, FDR < 0.05),
                  size = 4, point.padding = 0.15, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom") +
  facet_wrap(~Contrast, scales = "free_y")
classcomparison
ggsave(paste0("~/Directory/NanostringTauProject/plots/ClassComparison.pdf"))

#export resultsclass as csv
write_xlsx(resultsClass,"~/Directory/NanostringTauProject/plots/ClassComparison.xlsx")
