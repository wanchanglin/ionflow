#' wl-03-07-2020, Fri: Load packages here. 
#' wl-06-07-2020, Mon: debug functions PreProcessing and ExploratoryAnalysis
library(IonFlow)

## ==== General settings ====
rm(list = ls(all = T))
setwd("~/my_galaxy/ionflow")

#' pkgs <- c("reshape","knitr","Matrix","gridExtra",
#'            "network", "igraph","psych","ggrepel","dplyr",)

#' wl-03-07-2020, Fri: qgraph loads plent of R packages
Packages <- 
  c("data.table","reshape2","tidyr", "ggplot2", 
    "corrplot","gplots","pheatmap", "factoextra","ggfortify","mixOmics",
    "GGally", "intergraph", "sna", "qgraph", 
    "org.Sc.sgd.db","GO.db","GOstats")
suppressWarnings(invisible(lapply(Packages, library, character.only = TRUE)))

pkgs <- c("data.table","reshape2","tidyr", "ggplot2","ggrepel", 
          "corrplot","gplots","pheatmap", "mixOmics", "qgraph")
invisible(lapply(pkgs, library,character.only = TRUE))

source("all_IonFlow.R")

#' IonData <- read.table("./test-data/IonData.txt", header = T, sep = " ", 
#'                       stringsAsFactors = F)
#' save(IonData,file="./test-data/IonData.rdata")
load(file="./test-data/IonData.rdata")

## ==== Pre-processing ====
# data=IonData
# stdev=NULL
pre_proc <- PreProcessing(data = IonData)

#' pre_proc <- PreProcessing(data = IonData, stdev = pre_defined_sd)
# stats
pre_proc$stats.raw_data
pre_proc$stats.outliers
pre_proc$stats.median_batch_corrected_data
pre_proc$stats.standardised_data
# plots
pre_proc$plot.logConcentration_by_batch
pre_proc$plot.logConcentration_z_scores
# data
head(pre_proc$dataR.long)
head(pre_proc$data.long)
head(pre_proc$data.wide)
head(pre_proc$data.wide_Symb)

# save(pre_proc,file="./test-data/pre_proc.rdata")

## ==== Load Pre-proceesed data from GitHub ====

#' load data set from github 
data <- read.csv("./test-data/data.wide.csv", stringsAsFactors = F)
data <- data[,-1]
data_Symb <- read.csv("./test-data/data.wide_Symb.csv", stringsAsFactors = F)
data_Symb <- data_Symb[,-1]

#' Use saved pre-processing results
#' load(file="./test-data/pre_proc.rdata")
load(file="./doc/pre_proc.rdata")

#' data = pre_proc$data.wide 
#' data_Symb = pre_proc$data.wide_Symb

## ==== Exploratory analysis ====
exp_anal <- ExploratoryAnalysis(data = data)

exp_anal <- ExploratoryAnalysis(data = pre_proc$data.wide)

# plots
exp_anal$plot.Pearson_correlation
exp_anal$plot.PCA_Individual
exp_anal$plot.heatmap
exp_anal$plot.pairwise_correlation_map
exp_anal$plot.regularized_partial_correlation_network
# data
head(exp_anal$data.PCA_loadings)

#' save(exp_anal,file="./doc/exp_anal.rdata")

## ==== Gene Clustering ====

gene_clust <- GeneClustering(data = data, data_Symb = data_Symb)

gene_clust <- GeneClustering(data = pre_proc$data.wide, 
                             data_Symb = pre_proc$data.wide_Symb)
# stats
gene_clust$stats.clusters
gene_clust$stats.Kegg_Goslim_annotation
gene_clust$stats.Goterms_enrichment
# plots
gene_clust$plot.profiles

#' save(gene_clust,file="./doc/gene_clust.rdata")

## ==== Gene Network ====

gene_net <- GeneNetwork(data=data, data_Symb=data_Symb)

gene_net <- GeneNetwork(data = pre_proc$data.wide, 
                        data_Symb = pre_proc$data.wide_Symb)
# stats
gene_net$stats.impact_betweeness
gene_net$stats.impact_betweeness_by_cluster
# plots
gene_net$plot.pnet
gene_net$plot.impact_betweenees

