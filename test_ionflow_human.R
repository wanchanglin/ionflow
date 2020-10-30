#' wl-11-10-2020, Sun: test PreProcessing

## ==== General settings ====
rm(list = ls(all = T)) 

tool_dir <- "~/my_galaxy/ionflow/"

setwd(tool_dir)
pkgs <- c("optparse", "reshape2", "plyr", "dplyr", "tidyr", "ggplot2",
          "ggrepel", "corrplot", "gplots", "network", "sna", "GGally",
          "org.Sc.sgd.db", "GO.db", "GOstats", "pheatmap")
invisible(lapply(pkgs, library, character.only = TRUE))
source("funcs_ionflow.R")

## ==== Data preparation ====

data <- read.table("./test-data/human_z_profiles.csv", header = T, sep = ",")
data <- data[!duplicated(data$Gene), ]
colnames(data)[1] <- "Line"
dim(data)

data_symb <- symbol_data(x = data, symb_thr = 2)

## ==== Exploratory analysis ====
exp_anal <- ExploratoryAnalysis(data = data)

## ==== Gene Network ====

gene_net <- GeneNetwork(data = data,
                        data_symb = data_symb,
                        min_clust_size = 5, thres_corr = 0.6,
                        method_corr = "cosine")
                        #' method_corr = "pearson")
                        #' method_corr = "hybrid_mahal_cosine")
                        #' method_corr = "mahal_cosine")
gene_net$plot.pnet

#' ==== GO/KEGG enrichment analysis ====
library("org.Hs.eg.db")
kegg_en <- kegg_enrich(data = data_symb, min_clust_size = 5,
                       pval = 0.05, annot_pkg =  "org.Hs.eg.db")
kegg_en

go_en  <- go_enrich(data = data_symb, min_clust_size = 10,
                    pval = 0.05, ont = "BP", annot_pkg =  "org.Hs.eg.db")
go_en
