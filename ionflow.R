#' wl-03-07-2020, Fri: Load packages here.
#' wl-03-07-2020, Fri: qgraph loads plenty of R packages
#' wl-06-07-2020, Mon: debug functions PreProcessing and ExploratoryAnalysis
#' wl-30-07-2020, Thu: test on subset of IonData
#' wl-07-08-2020, Fri: re-write for galaxy
#' wl-08-08-2020, Sat: works for both in comand line and interactive mode
#' wl-09-08-2020, Sun: recordPlot has problem for command line mode. 
#' wl-10-08-2020, Mon: sort out the base graphics saving via non-interactive
#'  mode.

## ==== General settings ====
rm(list = ls(all = T))

#' flag for command-line use or not. If false, only for debug interactively.
com_f <- T

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

pkgs <- c("optparse", "reshape2", "plyr", "dplyr", "tidyr", "ggplot2", "ggrepel",
          "corrplot", "gplots", "network", "sna", "GGally", 
          "org.Sc.sgd.db", "GO.db", "GOstats")
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
                centers large than threshold will be kept."
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
        help = "Save clustering stats table"
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
        type = "character", default = "impact_betweeness.tsv",
        help = "Save impact and betweeness table"
      ),
      make_option("--imbe_tab_out",
        type = "character", default = "impact_betweeness_tab.tsv",
        help = "Save impact and betweeness contigency table"
      )
    )

  opt <- parse_args(
    object = OptionParser(option_list = option_list),
    args = commandArgs(trailingOnly = TRUE)
  )
} else {
  #' tool_dir <- "C:/R_lwc/my_galaxy/ionflow/"         #' for windows
  #' tool_dir <- "~/my_galaxy/ionflow/"   #' for linux. Must be case-sensitive
  tool_dir <- "~/R_lwc/r_data/icl/"       #' for linux. Must be case-sensitive

  opt <- list(

    #' Input
    ion_file = paste0(tool_dir, "test-data/iondata_test.tsv"),
    std_file_sel = "no",
    std_file = paste0(tool_dir, "test-data/user_std.tsv"),

    #' Clustering and network analysis
    thres_clus = 10.0, 
    thres_anno = 5.0, 
    thres_corr = 0.6, 

    #' output: pre-processing
    pre_proc_pdf       = paste0(tool_dir, "test-data/res/pre_proc.pdf"),
    df_stats_out       = paste0(tool_dir, "test-data/res/df_stats.tsv"),
    outl_out           = paste0(tool_dir, "test-data/res/outl.tsv"),
    data_wide_out      = paste0(tool_dir, "test-data/res/data_wide.tsv"),
    data_wide_symb_out = paste0(tool_dir, "test-data/res/data_wide_symb.tsv"),

    #' output: exploratory analysis
    exp_anal_pdf  = paste0(tool_dir, "test-data/res/exp_anal.pdf"),

    #' output: gene clustering
    gene_clus_pdf = paste0(tool_dir, "test-data/res/gene_clus.pdf"),
    clus_out      = paste0(tool_dir, "test-data/res/clus.tsv"),
    anno_out      = paste0(tool_dir, "test-data/res/kegg_go_anno.tsv"),
    enri_out      = paste0(tool_dir, "test-data/res/go_enri.tsv"),

    #' output: gene network
    gene_net_pdf = paste0(tool_dir, "test-data/res/gene_net.pdf"),
    imbe_out     = paste0(tool_dir, "test-data/res/impact_betweeness.tsv"),
    imbe_tab_out = paste0(tool_dir, "test-data/res/impact_betweeness_tab.tsv")
  )
}
#' print(opt)

suppressPackageStartupMessages({
  source(paste0(tool_dir, "funcs_ionflow.R"))
})

## ==== Data preparation ====

#' data for annotations
lib_dir <- paste0(tool_dir, "libraries/")
data_GOslim <- read.table(paste(lib_dir, "data_GOslim.tsv", sep = "/"), 
                          sep = "\t", header = T)
data_ORF2KEGG <- read.table(paste(lib_dir, "data_ORF2KEGG.tsv", sep = "/"), 
                            sep = "\t", header = T)

#' Load data set
ion_data <- read.table(opt$ion_file, header = T, sep = "\t")

if (opt$std_file_sel == "yes") {
  std_data <- read.table(opt$std_file, header = T, sep = "\t")
} else {
  std_data <- NULL
}

## ==== Pre-processing ====

pre_proc <- PreProcessing(data = ion_data, stdev = std_data)

#' save plot in pdf
pdf(file = opt$pre_proc_pdf, onefile = T, width=15, height=10)
plot(pre_proc$plot.logConcentration_by_batch)
plot(pre_proc$plot.logConcentration_z_scores)
dev.off()

#' bind stats
df_stats <- list(raw_data = pre_proc$stats.raw_data,
                 bat_corr_data = pre_proc$stats.median_batch_corrected_data,
                 std_data = pre_proc$stats.standardised_data)
df_stats <- dplyr::bind_rows(df_stats, .id = "Data_Set")
row.names(df_stats) = NULL 

#' save tables
write.table(df_stats, file = opt$df_stats_out, sep = "\t", row.names = F)
write.table(pre_proc$stats.outliers, file = opt$outl_out, sep = "\t", 
            row.names = F)
write.table(pre_proc$data.wide, file = opt$data_wide_out, sep = "\t", 
            row.names = F)
write.table(pre_proc$data.wide_symb, file = opt$data_wide_symb_out, 
            sep = "\t", row.names = F)

## ==== Exploratory analysis ====

#' wl-10-08-2020, Mon: Base graphics saving does not work for 
#' non-interactive mode. use this dirt trick.
pdf(file = opt$exp_anal_pdf, onefile = T) # ,width=15, height=10)
exp_anal <- ExploratoryAnalysis(data = pre_proc$data.wide)

## dev.control(displaylist="enable")
## exp_anal$plot.Pearson_correlation
## exp_anal$plot.heatmap
## exp_anal$plot.pairwise_correlation_map
exp_anal$plot.correlation_network
exp_anal$plot.PCA_Individual
dev.off()

## ==== Gene Clustering ====
gene_clus <- GeneClustering(data = pre_proc$data.wide,
                            data_symb = pre_proc$data.wide_symb,
                            thres_clus = opt$thres_clus, 
                            thres_anno = opt$thres_anno)

pdf(file = opt$gene_clus_pdf, onefile = T, width=15, height=10)
gene_clus$plot.profiles
dev.off()

write.table(gene_clus$stats.clusters, file = opt$clus_out, 
            sep = "\t", row.names = FALSE)
write.table(gene_clus$stats.Kegg_Goslim_annotation, file = opt$anno_out, 
            sep = "\t", row.names = FALSE)
write.table(gene_clus$stats.Goterms_enrichment, file = opt$enri_out, 
            sep = "\t", row.names = FALSE)

## ==== Gene Network ====
gene_net <- GeneNetwork(data = pre_proc$data.wide,
                        data_symb = pre_proc$data.wide_symb,
                        thres_clus = opt$thres_clus, 
                        thres_corr = opt$thres_corr)

pdf(file = opt$gene_net_pdf, onefile = T, width=15, height=10)
gene_net$plot.pnet
gene_net$plot.impact_betweeness
dev.off()

write.table(gene_net$stats.impact_betweeness, file = opt$imbe_out, 
            sep = "\t", row.names = FALSE)
write.table(gene_net$stats.impact_betweeness_tab, file = opt$imbe_tab_out, 
            sep = "\t", row.names = FALSE)

