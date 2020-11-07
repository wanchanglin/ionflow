#' wl-30-10-2020, Fri: human data
#' wl-06-11-2020, Fri: test enrichment analysis for human database

## ==== General settings ====
rm(list = ls(all = T)) 

#' tool_dir <- "~/my_galaxy/ionflow/"
tool_dir <- "C:/R_lwc/my_galaxy/ionflow/"

setwd(tool_dir)
pkgs <- c("optparse", "reshape2", "plyr", "dplyr", "tidyr", "ggplot2",
          "ggrepel", "corrplot", "gplots", "network", "sna", "GGally",
          "org.Hs.eg.db", "GO.db", "GOstats", "pheatmap")
invisible(lapply(pkgs, library, character.only = TRUE))
source("funcs_ionflow.R")

## ==== Data preparation ====
dat <- read.table("./test-data/human_z_profiles.csv", header = T, sep = ",")
dat <- dat[!duplicated(dat$Gene), ]
colnames(dat)[1] <- "Line"
dat_symb <- symbol_data(x = dat, symb_thr = 2)

#' filtering: Select phenotypes of interest
if (F) {
  idx      <- rowSums(abs(dat_symb[, -1])) > 0
  dat      <- dat[idx, ]
  dat_symb <- dat_symb[idx, ]
}

dim(dat)
dim(dat_symb)

## ==== Exploratory analysis ====
exp_anal <- ExploratoryAnalysis(data = dat)

## ==== Gene Network ====
gene_net <- GeneNetwork(data = dat,
                        data_symb = dat_symb,
                        min_clust_size = 5, thres_corr = 0.6,
                        method_corr = "pearson")
                        #' method_corr = "cosine")
                        #' method_corr = "hybrid_mahal_cosine")
                        #' method_corr = "mahal_cosine")
gene_net$plot.pnet

#' ==== GO/KEGG enrichment analysis ====
kegg_en <- kegg_enrich(data = dat_symb, min_clust_size = 5,
                       pval = 0.05, annot_pkg = "org.Hs.eg.db")
kegg_en

go_en  <- go_enrich(data = data_symb, min_clust_size = 5,
                    pval = 0.05, ont = "BP", annot_pkg = "org.Hs.eg.db")
go_en

