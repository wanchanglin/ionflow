#' wl-11-10-2020, Sun: test PreProcessing
#' wl-19-10-2020, Mon: test and debug jacopo's changes. Fix some bugs.
#' wl-21-10-2020, Wed: test enrichment analysis
#' wl-26-10-2020, Mon: test GeneNetwork with cosM
#' wl-30-10-2020, Fri: test some changes

## ==== General settings ====
rm(list = ls(all = T)) 

tool_dir <- "~/my_galaxy/ionflow/"
#' tool_dir <- "C:/R_lwc/my_galaxy/ionflow/"

setwd(tool_dir)
pkgs <- c("optparse", "reshape2", "plyr", "dplyr", "tidyr", "ggplot2",
          "ggrepel", "corrplot", "gplots", "network", "sna", "GGally",
          "org.Sc.sgd.db", "GO.db", "GOstats", "pheatmap") #, "pracma")
invisible(lapply(pkgs, library, character.only = TRUE))
source("funcs_ionflow.R")

## ==== Data preparation ====

#' ion_data <- read.table("./test-data/iondata_test.tsv", header = T, sep = "\t")
#' ion_data <- read.table("./test-data/iondata.tsv", header = T, sep = "\t")
#' ion_data <- read.table("C:/R_lwc/r_data/icl/test-data/ionome_ko.tsv", header = T, sep = "\t")
ion_data <- read.table("~/R_lwc/r_data/icl/test-data/ionome_ko_test.tsv", header = T, sep = "\t")
#' ion_data <- read.table("~/R_lwc/r_data/icl/test-data/ionome_oe.tsv", header = T, sep = "\t")
#' ion_data <- read.table("~/R_lwc/r_data/icl/test-data/ionome_oe_test.tsv", header = T, sep = "\t")

## ==== Pre-processing ====
pre_proc <- PreProcessing(data = ion_data,
                          var_id = 1, batch_id = 4, data_id = 5,
                          method_norm = "median",
                          method_outliers = "IQR",
                          n_thrs = 3,
                          stand_method = "std",
                          stdev = NULL, symb_thr = 3)

head(pre_proc$data.gene.zscores)
head(pre_proc$data.gene.symb)
#' pre_proc$stats.raw_data
#' pre_proc$stats.outliers
#' pre_proc$stats.batch_data
#' head(pre_proc$data.long)
#' head(pre_proc$data.gene.logFC)
#' pre_proc$plot.dot
#' pre_proc$plot.hist

#' ==== Filter data set ====
data      <- pre_proc$data.gene.zscores
data_symb <- pre_proc$data.gene.symb
dim(data)

#' Select phenotypes of interest
idx       <- rowSums(abs(data_symb[, -1])) > 0
data      <- data[idx, ]
data_symb <- data_symb[idx, ]
dim(data)

## ==== Exploratory analysis ====
exp_anal <- ExploratoryAnalysis(data = data)
#' exp_anal$plot.Pearson_correlation
#' exp_anal$plot.PCA_Individual
#' exp_anal$plot.heatmap
#' exp_anal$plot.pairwise_correlation_map
#' exp_anal$plot.correlation_network
#' head(exp_anal$data.PCA_loadings)

## ==== Gene Network ====

gene_net <- GeneNetwork(data = data,
                        data_symb = data_symb,
                        min_clust_size = 5, thres_corr = 0.6,
                        #' method_corr = "cosine")
                        method_corr = "pearson")
                        #' method_corr = "hybrid_mahal_cosine")
                        #' method_corr = "mahal_cosine")
gene_net$plot.pnet
gene_net$plot.impact_betweenness
gene_net$stats.impact_betweenness
gene_net$stats.impact_betweenness_tab

#' ==== GO/KEGG enrichment analysis ====
kegg_en <- kegg_enrich(data = data_symb, min_clust_size = 10,
                       pval = 0.05)
kegg_en

go_en  <- go_enrich(data = data_symb, min_clust_size = 10,
                    pval = 0.05, ont = "BP")
go_en

## ==== Gene Clustering ====
#' data for annotations
lib_dir <- paste0(tool_dir, "libraries/")
data_GOslim <- read.table("./libraries/data_GOslim.tsv", sep = "\t", header = T)
data_ORF2KEGG <- read.table("./libraries/data_ORF2KEGG.tsv", sep = "\t", header = T)

gene_clust <- GeneClustering(data = data,
                             data_symb = data_symb,
                             min_clust_size = 10, thres_anno = 5)
gene_clust$plot.profiles
gene_clust$stats.clusters
gene_clust$stats.Kegg_Goslim_annotation
gene_clust$stats.Goterms_enrichment

## =========================================================================
#' wl-21-10-2020, Wed: DEBUG stuff

if (F) {

  #' for PreProcessing
  data = ion_data
  var_id = 1
  batch_id = 2
  data_id = 3
  method_norm = "median"
  #' control_lines =  "BY4741"    #' only for inome_ko'
  control_lines = NULL
  control_use = "control"
  method_outliers = "mad"
  n_thrs = 3
  stand_method = "std"
  stdev = NULL
  symb_thr = 4

  #' For GeneNetwork
  data = pre_proc$data.gene.zscores
  data_symb = pre_proc$data.gene.symb
  min_clust_size = 10
  thres_corr = 0.6
  method_corr = c("pearson", "spearman", "kendall", "cosine",
                  "mahal_cosine", "hybrid_mahal_cosine")

  #' for enrichment analysis
  data            <- pre_proc$data.gene.symb
  min_clust_size  <- 5
  annot_pkg       <- "org.Sc.sgd.db"
  pval            <- 0.05
  ont             <- "BP"
  f_max           <- F
  x               <- data[, -1]

  gene_clus(x, min_clust_size = 5, f_max = F)
}
