#' wl-11-10-2020, Sun: test PreProcessing

## ==== General settings ====
rm(list = ls(all = T))

tool_dir <- "~/my_galaxy/ionflow/"       #' wl: change here

setwd(tool_dir)
pkgs <- c("optparse", "reshape2", "plyr", "dplyr", "tidyr", "ggplot2", 
          "ggrepel", "corrplot", "gplots", "network", "sna", "GGally",
          "org.Sc.sgd.db", "GO.db", "GOstats")
invisible(lapply(pkgs, library, character.only = TRUE))
source("funcs_ionflow_ji.R")

## ==== Data preparation ====

ion_data <- read.table("./test-data/iondata.tsv", header = T, sep = "\t")
#' ion_data <- read.table("./test-data/ionome_oe.tsv", header = T, sep = "\t")
#' ion_data <- read.table("./test-data/ionome_ko.tsv", header = T, sep = "\t")

## ==== Pre-processing ====
pre_proc <- PreProcessing(data = ion_data, 
                          var_id = 1, batch_id = 2, data_id = 3,
                          method_norm = "median",
                          control_lines = NULL,
                          control_use = "control",
                          method_outliers = "mad",
                          n_thrs = 4,
                          stand_method = "std",
                          stdev = NULL, symb_thr = 4, p_symb = 0.5)


pre_proc$stats.raw_data
pre_proc$stats.outliers
pre_proc$stats.batch_data
pre_proc$stats.stand_data
head(pre_proc$data.long_bat)
head(pre_proc$data.long)
head(pre_proc$data.wide)
head(pre_proc$data.wide_symb)
pre_proc$plot.dot
pre_proc$plot.hist

## ==== Exploratory analysis ====
exp_anal <- ExploratoryAnalysis(data = pre_proc$data.wide)

exp_anal$plot.Pearson_correlation
exp_anal$plot.PCA_Individual
exp_anal$plot.heatmap
exp_anal$plot.pairwise_correlation_map
exp_anal$plot.correlation_network
head(exp_anal$data.PCA_loadings)

## ==== Gene Clustering ====
#' data for annotations
lib_dir <- paste0(tool_dir, "libraries/")
data_GOslim <- read.table("./libraries/data_GOslim.tsv", sep = "\t", header = T)
data_ORF2KEGG <- read.table("./libraries/data_ORF2KEGG.tsv", sep = "\t", header = T)

gene_clust <- GeneClustering(data = pre_proc$data.wide,
                             data_symb = pre_proc$data.wide_symb,
                             thres_clus = 10, thres_anno = 5)
gene_clust$plot.profiles
gene_clust$stats.clusters
gene_clust$stats.Kegg_Goslim_annotation
gene_clust$stats.Goterms_enrichment

## ==== Gene Network ====
gene_net <- GeneNetwork(data = pre_proc$data.wide,
                        data_symb = pre_proc$data.wide_symb,
                        thres_clus = 10, thres_cor = 0.75)

gene_net$plot.pnet
gene_net$plot.impact_betweenness
gene_net$stats.impact_betweenness
gene_net$stats.impact_betweenness_tab
