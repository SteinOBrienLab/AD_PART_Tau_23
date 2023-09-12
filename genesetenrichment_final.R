##-- 1/14/2023 Ryan Palaganas
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(CoGAPS)
library(SingleCellExperiment)
library(ggplot2)
library(org.Hs.eg.db)
library(AnnotationDbi)

#load in CoGAPS object----
result <- readRDS("~/Directory/NanostringTauProject/cache/CoGAPSnp5.rds")
cogapsRes<-list("np5"=result)

#creating ranked gene lists for each pattern---- 
#extract feature loading matrix
patterns <- as.data.frame(cogapsRes[["np5"]]@featureLoadings)
#attach entrezid column to the pattern matrix
entrezids <- mapIds(org.Hs.eg.db, rownames(patterns), 'ENTREZID', 'SYMBOL')
patterns <- cbind(patterns, entrezids)
#remove genes that do not have an entrezid map
patterns <- na.omit(patterns)

#create rank vectors for each pattern to input into fgsea
rank1 <- setNames(patterns[,1], as.character(patterns$entrezids))
rank2 <- setNames(patterns[,2], as.character(patterns$entrezids))
rank3 <- setNames(patterns[,3], as.character(patterns$entrezids))
rank4 <- setNames(patterns[,4], as.character(patterns$entrezids))
rank5 <- setNames(patterns[,5], as.character(patterns$entrezids))

library(readxl)
library(fgsea)
#load in synapse enrichment genes, this can also be used to introduce any gene set of interest outside of what is available
#from fgsea 'examplepathways' gene set 
#bring examplepathways into environment
data(examplePathways)
#load synaptic pathway
Syngo <- read_excel("~/Directory/SynGoAnalysisOverall (1).xlsx")
synapsepathway <- Syngo$`gene symbol`
#attach entrezid column
Syngo$ENTREZID <- mapIds(org.Hs.eg.db, synapsepathway, 'ENTREZID','SYMBOL')
length(Syngo$ENTREZID)

#create a list of genes
synapticgenes <- Syngo$ENTREZID %>% as.character()
pathways <- list(synapticgenes)
#load synaptic pathway into example pathways
examplePathways$syngo_genes <- synapticgenes

#create fgseaRes object based on pattern 4 (tau up)
fgseaRes <- fgsea(pathways = examplePathways, 
                  stats = rank4,
                  eps = 0,
                  minSize=15,
                  maxSize=500,
                  scoreType = "pos")

class(fgseaRes)

# top 6 enriched pathways
head(fgseaRes[order(pval), ])

# number of significant pathways at padj < 0.01
sum(fgseaRes[, padj < 0.01])

#plot significantly enriched pathway----
pdf("~/Directory/NanostringTauProject/plots/synapticgeneenrichment.pdf",width=10,height=5)
plotEnrichment(examplePathways[["syngo_genes"]],
               rank4) + labs(title="syngo_genes")
dev.off()
