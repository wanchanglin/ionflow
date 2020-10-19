#' wl-11-10-2020, Sun: test PreProcessing
#' wl-19-10-2020, Mon: test and debug jacopo's changes. Fix some bugs.

## ==== General settings ====
rm(list = ls(all = T))

tool_dir <- "~/my_galaxy/ionflow/"       #' wl: change here

setwd(tool_dir)
pkgs <- c("optparse", "reshape2", "plyr", "dplyr", "tidyr", "ggplot2", 
          "ggrepel", "corrplot", "gplots", "network", "sna", "GGally",
          "org.Sc.sgd.db", "GO.db", "GOstats","pracma")
invisible(lapply(pkgs, library, character.only = TRUE))
source("funcs_ionflow.R")

## ==== Data preparation ====

#' ion_data <- read.table("./test-data/iondata.tsv", header = T, sep = "\t")
#' ion_data <- read.table("./test-data/ionome_oe.tsv", header = T, sep = "\t")
ion_data <- read.table("~/R_lwc/r_data/icl/test-data/ionome_ko.tsv", header = T, sep = "\t")

## ==== Pre-processing ====
pre_proc <- PreProcessing(data = ion_data, 
                          var_id = 1, batch_id = 3, data_id = 5,
                          method_norm = "median",
                          control_lines =  "BY4741",
                          control_use = "control",
                          method_outliers = "mad",
                          n_thrs = 3,
                          stand_method = "std",
                          stdev = NULL, symb_thr = 4, p_symb = 0.5)

#' names(pre_proc)
#' [1] "stats.raw_data"    "stats.outliers"    "stats.batch_data"
#' [4] "data.long"         "data.gene.logFC"   "data.gene.zscores"
#' [7] "data.gene.symb"    "plot.dot"          "plot.hist"

pre_proc$stats.raw_data
pre_proc$stats.outliers
pre_proc$stats.batch_data
head(pre_proc$data.long)
head(pre_proc$data.gene.logFC)
head(pre_proc$data.gene.zscores)
head(pre_proc$data.gene.symb)
#' pre_proc$plot.dot             #' wl-19-10-2020, Mon: problem. 
#' pre_proc$plot.hist            #' wl-19-10-2020, Mon: problem. 

## ==== Exploratory analysis ====
exp_anal <- ExploratoryAnalysis(data = pre_proc$data.gene.zscores)

exp_anal$plot.Pearson_correlation
exp_anal$plot.PCA_Individual
exp_anal$plot.heatmap
exp_anal$plot.pairwise_correlation_map
exp_anal$plot.correlation_network
head(exp_anal$data.PCA_loadings)

## ==== Gene Network ====
#' data = pre_proc$data.gene.zscores
#' data_symb = pre_proc$data.gene.symb
#' method_corr = "pearson"
#' min_clust_size = 10
#' thres_corr = 0.6

gene_net <- GeneNetwork(data = pre_proc$data.gene.zscores,
                        data_symb = pre_proc$data.gene.symb,
                        min_clust_size = 10, thres_corr = 0.6)

gene_net$plot.pnet
gene_net$plot.impact_betweenness
gene_net$stats.impact_betweenness
gene_net$stats.impact_betweenness_tab

## ==== Gene Clustering ====
#' data for annotations
lib_dir <- paste0(tool_dir, "libraries/")
data_GOslim <- read.table("./libraries/data_GOslim.tsv", sep = "\t", header = T)
data_ORF2KEGG <- read.table("./libraries/data_ORF2KEGG.tsv", sep = "\t", header = T)

gene_clust <- GeneClustering(data = pre_proc$data.gene.zscores,
                             data_symb = pre_proc$data.gene.symb,
                             thres_clus = 10, thres_anno = 5)
gene_clust$plot.profiles
gene_clust$stats.clusters
gene_clust$stats.Kegg_Goslim_annotation
gene_clust$stats.Goterms_enrichment

