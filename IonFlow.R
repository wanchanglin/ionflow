#' wl-03-07-2020, Fri: Load packages here. 
#' wl-06-07-2020, Mon: debug functions PreProcessing and ExploratoryAnalysis
#' wl-06-07-2020, Mon: R packages 
#' - Not used: "intergraph", "factoextra", "ggfortify", "knitr", "reshape"
#' - attached: "Matrix", "network", "igraph","psych", "data.table",
#'             "gridExtra","GGally"
#' - option: "mixOmics" 
#' wl-13-07-2020, Mon: NAs and PCA
#' wl-14-07-2020, Tue: Fix a bug in network analysis
#' wl-23-07-2020, Thu: test PreProcessing with NAs removal

## ==== General settings ====
rm(list = ls(all = T))
#' setwd("~/my_galaxy/ionflow")
setwd("C:/R_lwc/my_galaxy/ionflow")

#' wl-03-07-2020, Fri: qgraph loads plent of R packages
pkgs <- c("reshape2", "plyr", "dplyr", "tidyr", "ggplot2", "ggrepel",
          "corrplot", "gplots", "pheatmap", 
          "network", "sna", "GGally", "qgraph",
          "org.Sc.sgd.db", "GO.db", "GOstats")
suppressWarnings(invisible(lapply(pkgs, library, character.only = TRUE)))

if (F) {
  library(IonFlow) 
} else {
  source("all_IonFlow.R")
  load(file="./test-data/IonData.rdata")
  load(file="./test-data/pre_defined_sd.rdata")
  load(file="./test-data/data_GOslim.rdata")
  load(file="./test-data/data_ORF2KEGG.rdata")
}

## ==== Pre-processing ====
#' data = IonData
#' stdev = pre_defined_sd

#' pre_proc <- PreProcessing(data = IonData, stdev = NULL)
pre_proc <- PreProcessing(data = IonData, stdev = pre_defined_sd)

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
head(pre_proc$data.wide_symb)

#' save(pre_proc,   file="./test-data/pre_proc.rdata")

## ==== Load Pre-proceesed data ====

#' load(file="./test-data/pre_proc.rdata")
data      <- pre_proc$data.wide 
data_symb <- pre_proc$data.wide_symb
sum(is.na(data))

## ==== Exploratory analysis ====
exp_anal <- ExploratoryAnalysis(data = pre_proc$data.wide)
# plots
exp_anal$plot.Pearson_correlation
exp_anal$plot.PCA_Individual
exp_anal$plot.heatmap
exp_anal$plot.pairwise_correlation_map
exp_anal$plot.regularized_partial_correlation_network
# data
head(exp_anal$data.PCA_loadings)

## ==== Gene Clustering ====
gene_clust <- GeneClustering(data = pre_proc$data.wide, 
                             data_symb = pre_proc$data.wide_symb)
# stats
gene_clust$stats.clusters
gene_clust$stats.Kegg_Goslim_annotation
gene_clust$stats.Goterms_enrichment
# plots
gene_clust$plot.profiles

## ==== Gene Network ====
gene_net <- GeneNetwork(data = pre_proc$data.wide, 
                        data_symb = pre_proc$data.wide_symb)
# stats
gene_net$stats.impact_betweeness
gene_net$stats.impact_betweeness_by_cluster
# plots
gene_net$plot.pnet
gene_net$plot.impact_betweenees

