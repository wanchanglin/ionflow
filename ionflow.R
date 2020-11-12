#' wl-03-07-2020, Fri: Load packages here.
#' wl-06-07-2020, Mon: debug functions PreProcessing and ExploratoryAnalysis
#' wl-07-08-2020, Fri: re-write for galaxy
#' wl-09-08-2020, Sun: recordPlot has problem for command line mode.
#' wl-02-09-2020, Wed: change for PreProcessing
#' wl-12-11-2020, Thu: change for version 2

## ==== General settings ====
rm(list = ls(all = T))

#' flag for command-line use or not. If false, only for debug interactively.
com_f <- F

#' galaxy will stop even if R has warning message
options(warn = -1) #' disable R warning. Turn back: options(warn=0)

#' ------------------------------------------------------------------------
#' Setup R error handling to go to stderr
#' options( show.error.messages=F, error = function (){
#'   cat( geterrmessage(), file=stderr() )
#'   q( "no", 1, F )
#' })

#' we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

#' wl-28-08-2018, Tue: Convert a string separated by comma into character vector
str_vec <- function(x) {
  x <- unlist(strsplit(x, ","))
  x <- gsub("^[ \t]+|[ \t]+$", "", x) #' trim white spaces
}

pkgs <- c("optparse", "reshape2", "plyr", "dplyr", "tidyr", "ggplot2",
          "ggrepel", "corrplot", "gplots", "network", "sna", "GGally",
          "org.Sc.sgd.db", "GO.db", "GOstats", "pheatmap")
suppressPackageStartupMessages(invisible(lapply(pkgs, library,
                                                character.only = TRUE)))

## ==== Command line or interactive setting ====
if (com_f) {

  #' -----------------------------------------------------------------------
  #' Setup home directory
  #' wl-24-11-2017, Fri: A dummy function for the base directory. The reason
  #' to write such a function is to keep the returned values by
  #' 'commandArgs' with 'trailingOnly = FALSE' in a local environment
  #' otherwise 'parse_args' will use the results of
  #' 'commandArgs(trailingOnly = FALSE)' even with 'args =
  #' commandArgs(trailingOnly = TRUE)' in its argument area.
  func <- function() {
    argv <- commandArgs(trailingOnly = FALSE)
    path <- sub("--file=", "", argv[grep("--file=", argv)])
  }
  #' prog_name <- basename(func())
  tool_dir <- paste0(dirname(func()), "/")

  option_list <-
    list(
      make_option(c("-v", "--verbose"),
        action = "store_true", default = TRUE,
        help = "Print extra output [default]"
      ),
      make_option(c("-q", "--quietly"),
        action = "store_false",
        dest = "verbose", help = "Print little output"
      ),

      #' input
      make_option("--ion_file",
        type = "character",
        help = "ion concentration file in tabular format"
      ),
      make_option("--var_id", type = "integer", default = 1,
                  help = "Column index of variable"),
      make_option("--batch_id", type = "integer", default = 2,
                  help = "Column index of batch ID"),
      make_option("--data_id", type = "integer", default = 3,
                  help = "Start column index of data matrix"),
      make_option("--std_file_sel",
        type = "character", default = "no",
        help = "Load user defined std file or not"
      ),
      make_option("--std_file",
        type = "character",
        help = "user predifined std file with respect to ions"
      ),

      #' Clustering and network analysis
      make_option("--thres_clus",
        type = "double", default = 10.0,
        help = "Clustering threshold for clustering and network. Clustering
                centres large than threshold will be kept."
      ),
      make_option("--thres_anno",
        type = "double", default = 5.0,
        help = "Percentage threshold for annotation (0 - 100).
                Features large than threshold will be kept."
      ),
      make_option("--thres_corr",
        type = "double", default = 0.60,
        help = "Correlation threshold for network analysis (0 - 1).
                Features large than threshold will be kept."
      ),

      #' output: pre-processing
      make_option("--pre_proc_pdf",
        type = "character", default = "pre_proc.pdf",
        help = "Save plots from pre-processing"
      ),
      make_option("--df_stats_out",
        type = "character", default = "df_stats.tsv",
        help = "Save stats summary of raw, batch corrected and
                standardised data"
      ),
      make_option("--outl_out",
        type = "character", default = "outl.tsv",
        help = "Save outliers summary"
      ),
      make_option("--data_wide_out",
        type = "character", default = "data_wide.tsv",
        help = "Save pre-processed data in wide format"
      ),
      make_option("--data_wide_symb_out",
        type = "character", default = "data_wide_symb.tsv",
        help = "Save pre-processed data Symbolization in wide format"
      ),

      #' output: exploratory analysis
      make_option("--exp_anal_pdf",
        type = "character", default = "exp_anal.pdf",
        help = "Save plots from exploratory analysis"
      ),

      #' output: gene clustering
      make_option("--gene_clus_pdf",
        type = "character", default = "gene_clus.pdf",
        help = "Save plots from gene clustering"
      ),
      make_option("--clus_out",
        type = "character", default = "clus.tsv",
        help = "Save clustering stats tableq
      ),
      make_option("--anno_out",
        type = "character", default = "kegg_go_anno.tsv",
        help = "Save Kegg and GO annotation table"
      ),
      make_option("--enri_out",
        type = "character", default = "go_enri.tsv",
        help = "Save GO terms enrichment table"
      ),

      #' output: gene network
      make_option("--gene_net_pdf",
        type = "character", default = "gene_net.pdf",
        help = "Save plots from gene network"
      ),
      make_option("--imbe_out",
        type = "character", default = "impact_betweenness.tsv",
        help = "Save impact and betweenness table"
      ),
      make_option("--imbe_tab_out",
        type = "character", default = "impact_betweenness_tab.tsv",
        help = "Save impact and betweenness contingency table"
      )
    )

  opt <- parse_args(
    object = OptionParser(option_list = option_list),
    args = commandArgs(trailingOnly = TRUE)
  )
} else {
  #' tool_dir <- "C:/R_lwc/my_galaxy/ionflow/"
  tool_dir <- "~/my_galaxy/ionflow/"

  opt <- list(

    #' Input
    ion_file = paste0(tool_dir, "test-data/iondata_test.tsv"),
    var_id = 1,
    batch_id = 2,
    data_id = 3,
    method_norm = "median",
    control_lines = NULL,
    control_use = "control",
    method_outliers = "IQR",
    thres_outl = 3,
    stand_method = "std",
    thres_symb = 2,
    #' stdev = NULL,
    std_file_sel = "no",
    std_file = paste0(tool_dir, "test-data/user_std.tsv"),

    #' network and enrichment analysis
    min_clust_size = 5.0,
    thres_corr = 0.6,
    method_corr = "pearson",
    pval = 0.05,
    ont = "BP", 
    annot_pkg =  "org.Sc.sgd.db",

    #' output: pre-processing
    pre_proc_pdf       = paste0(tool_dir, "test-data/res/pre_proc.pdf"),
    df_stats_out       = paste0(tool_dir, "test-data/res/df_stats.tsv"),
    outl_out           = paste0(tool_dir, "test-data/res/outl.tsv"),
    data_wide_out      = paste0(tool_dir, "test-data/res/data_wide.tsv"),
    data_wide_symb_out = paste0(tool_dir, "test-data/res/data_wide_symb.tsv"),

    #' output: exploratory analysis
    exp_anal_pdf  = paste0(tool_dir, "test-data/res/exp_anal.pdf"),

    #' output: gene network
    gene_net_pdf = paste0(tool_dir, "test-data/res/gene_net.pdf"),
    imbe_out     = paste0(tool_dir, "test-data/res/impact_betweenness.tsv"),
    imbe_tab_out = paste0(tool_dir, "test-data/res/impact_betweenness_tab.tsv"),

    #' output: enrichment analysis
    kegg_en_out = paste0(tool_dir, "test-data/res/kegg_en.tsv"),
    go_en_out   = paste0(tool_dir, "test-data/res/go_en.tsv")
  )
}
print(opt)

suppressPackageStartupMessages({
  source(paste0(tool_dir, "funcs_ionflow.R"))
})

## ==== Data preparation ====

#' Load data set
ion_data <- read.table(opt$ion_file, header = T, sep = "\t")

if (opt$std_file_sel == "yes") {
  std_data <- read.table(opt$std_file, header = T, sep = "\t")
} else {
  std_data <- NULL
}

## ==== Pre-processing ====

pre_proc <- PreProcessing(data            = ion_data,
                          var_id          = opt$var_id,
                          batch_id        = opt$batch_id,
                          data_id         = opt$data_id,
                          method_norm     = opt$method_norm,
                          control_lines   = opt$control_lines,
                          control_use     = opt$control_use,
                          method_outliers = opt$method_outliers,
                          thres_outl      = opt$thres_outl,
                          stand_method    = opt$stand_method,
                          stdev           = std_data,
                          thres_symb      = opt$thres_symb)

#' save plot in pdf
pdf(file = opt$pre_proc_pdf, onefile = T, width=15, height=10)
plot(pre_proc$plot.dot)
plot(pre_proc$plot.hist)
dev.off()

#' bind stats
df_stats <- list(raw_data = pre_proc$stats.raw_data,
                 bat_data = pre_proc$stats.batch_data)
df_stats <- dplyr::bind_rows(df_stats, .id = "Data_Set")
row.names(df_stats) = NULL

#' save tables
write.table(df_stats, file = opt$df_stats_out, sep = "\t", row.names = F)
write.table(pre_proc$stats.outliers, file = opt$outl_out, sep = "\t",
            row.names = F)
write.table(pre_proc$data.gene.zscores, file = opt$data_wide_out, sep = "\t",
            row.names = F)
write.table(pre_proc$data.gene.symb, file = opt$data_wide_symb_out,
            sep = "\t", row.names = F)

#' ==== Filter data set ====
dat      <- pre_proc$data.gene.zscores
dat_symb <- pre_proc$data.gene.symb

#' Select phenotypes of interest
idx      <- rowSums(abs(dat_symb[, -1])) > 0
dat      <- dat[idx, ]
dat_symb <- dat_symb[idx, ]

## ==== Exploratory analysis ====
pdf(file = opt$exp_anal_pdf, onefile = T) # ,width=15, height=10)
exp_anal <- ExploratoryAnalysis(data = dat)
exp_anal$plot.correlation_network
exp_anal$plot.PCA_Individual
dev.off()

## ==== Gene Network ====
gene_net <- GeneNetwork(data           = dat,
                        data_symb      = dat_symb,
                        min_clust_size = opt$min_clust_size,
                        thres_corr     = opt$thres_corr,
                        method_corr    = opt$method_corr)

pdf(file = opt$gene_net_pdf, onefile = T, width=15, height=10)
gene_net$plot.pnet1
gene_net$plot.pnet2
gene_net$plot.impact_betweenness
dev.off()

write.table(gene_net$stats.impact_betweenness, file = opt$imbe_out,
            sep = "\t", row.names = FALSE)
write.table(gene_net$stats.impact_betweenness_tab, file = opt$imbe_tab_out,
            sep = "\t", row.names = FALSE)

#' ==== GO/KEGG enrichment analysis ====
kegg_en  <- kegg_enrich(data           = dat_symb,
                        min_clust_size = opt$min_clust_size,
                        pval           = opt$pval,
                        annot_pkg      = opt$annot_pkg)

if (nrow(kegg_en) > 0) {
  write.table(kegg_en, file = opt$kegg_en_out, sep = "\t", row.names = FALSE)
}

go_en  <- go_enrich(data           = dat_symb,
                    min_clust_size = opt$min_clust_size,
                    pval           = opt$pval,
                    ont            = opt$ont,
                    annot_pkg      = opt$annot_pkg)

if (nrow(go_en) > 0) {
  write.table(go_en, file = opt$go_en_out, sep = "\t", row.names = FALSE)
}
