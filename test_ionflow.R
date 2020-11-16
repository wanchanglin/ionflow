#' wl-11-10-2020, Sun: test PreProcessing
#' wl-26-10-2020, Mon: test cosM
#' wl-07-11-2020, Sat: test enrichment
#' wl-09-11-2020, Mon: test GeneNetwork

## ==== General settings ====
rm(list = ls(all = T))

tool_dir <- "~/my_galaxy/ionflow/"

setwd(tool_dir)
pkgs <- c("optparse", "reshape2", "plyr", "dplyr", "tidyr", "ggplot2",
          "ggrepel", "corrplot", "gplots", "network", "sna", "GGally",
          "org.Sc.sgd.db", "GO.db", "GOstats", "KEGG.db", "pheatmap")
invisible(lapply(pkgs, library, character.only = TRUE))
source("funcs_ionflow.R")

## ==== Data preparation ====

#' ion_data <- read.table("./test-data/iondata_test.tsv", header = T, sep = "\t")
#' ion_data <- read.table("./test-data/iondata.tsv", header = T, sep = "\t")
#' ion_data <- read.table("~/R_lwc/r_data/icl/test-data/ionome_ko_test.tsv", header = T, sep = "\t")
ion_data <- read.table("~/R_lwc/r_data/icl/test-data/ionome_oe_test.tsv", header = T, sep = "\t")

#' Test for batch control
#' idx <- ion_data[, 1] %in% "BY4741"
#' sum(idx)
#' control_lines =  "BY4741"    #' only for inome_ko'
#' control_use = "all",


#' Test for cumtom standardisation
#' stdev <- read.table("./test-data/user_std.tsv", header = T, sep = "\t")
#' stand_method = "custom"
#' stdev = stdev

## ==== Pre-processing ====
pre <- PreProcessing(data = ion_data,
                     var_id = 1, batch_id = 4, data_id = 5,
                     method_norm = "median",
                     control_lines = NULL,
                     method_outliers = "IQR",
                     thres_outl = 3,
                     stand_method = "std",
                     stdev = NULL,
                     thres_symb = 2)

head(pre$data.gene.zscores)
head(pre$data.gene.symb)

#' ==== Filter data set ====
dat      <- pre$data.gene.zscores
dat_symb <- pre$data.gene.symb
dim(dat)

#' Select phenotypes of interest
idx      <- rowSums(abs(dat_symb[, -1])) > 0
dat      <- dat[idx, ]
dat_symb <- dat_symb[idx, ]
dim(dat)

## ==== Gene Network ====
gene_net <- GeneNetwork(data = dat,
                        data_symb = dat_symb,
                        min_clust_size = 10,
                        thres_corr = 0.75,
                        method_corr = "pearson")
                        #' method_corr = "cosine")
                        #' method_corr = "hybrid_mahal_cosine")
                        #' method_corr = "mahal_cosine")

gene_net$plot.pnet1    #' symbolic pheno (hclust)
X11()
gene_net$plot.pnet2    #' network (comminuty detection)

gene_net$plot.impact_betweenness
gene_net$stats.impact_betweenness
gene_net$stats.impact_betweenness_tab

#' ==== GO/KEGG enrichment analysis ====
kegg_en <- kegg_enrich(data = dat_symb, min_clust_size = 10, pval = 0.05,
                       annot_pkg =  "org.Sc.sgd.db")
kegg_en

go_en  <- go_enrich(data = dat_symb, min_clust_size = 10, pval = 0.05,
                    ont = "BP", annot_pkg =  "org.Sc.sgd.db")
go_en

## ==== Exploratory analysis ====
exp_anal <- ExploratoryAnalysis(data = dat)
#' exp_anal$plot.Pearson_correlation
#' exp_anal$plot.PCA_Individual
#' exp_anal$plot.heatmap
#' exp_anal$plot.pairwise_correlation_map
#' exp_anal$plot.correlation_network
#' head(exp_anal$data.PCA_loadings)

