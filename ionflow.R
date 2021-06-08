#' wl-07-06-2021, Mon: The fourth version: based on Jacopo's new changes in 
#' 'ionflow_funcs.R' and new pipeline 'tutorial_galaxy_ionflow.R'
#' wl-08-06-2021, Tue: finalise

## ==== General settings ====
rm(list = ls(all = T))

#' flag for command-line use or not. If false, only for debug interactively.
com_f <- T

#' galaxy will stop even if R has warning message
options(warn = -1) #' disable R warning. Turn back: options(warn=0)

#' ------------------------------------------------------------------------
#' Setup R error handling to go to stderr
#' options( show.error.messages=F, error = function () {
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
          "org.Sc.sgd.db","org.Hs.eg.db","GO.db", "GOstats", "KEGG.db",
          "pheatmap")
suppressPackageStartupMessages(invisible(lapply(pkgs, library,
                                                character.only = TRUE)))

## ==== Command line or interactive setting ====
if (com_f) {

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

      #' Data pre-processing
      make_option("--ion_file", type = "character",
                  help = "ion concentration file in tabular format"),
      make_option("--var_id", type = "integer", default = 1,
                  help = "Column index of variable"),
      make_option("--batch_id", type = "integer", default = 3,
                  help = "Column index of batch ID"),
      make_option("--data_id", type = "integer", default = 5,
                  help = "Start column index of data matrix"),
      make_option("--method_norm", type = "character", default = "median",
                  help = "Batch correction methods.  Support: median,
                          median+std and none"),
      make_option("--batch_control", type = "character", default = "yes",
                  help = "Use control lines for batch correction or not"),
      make_option("--control_lines", type = "character", default = "BY4741",
                  help = "Batch control lines"),
      make_option("--control_use", type = "character", default = "all",
                  help = "Select lines used for batch correction control.
                          Three selection: control, all and control.out"),
      make_option("--method_outliers", type = "character",
                  default = "log.FC.dist",
                  help = "Outlier detection method. Currently support:
                          mad, IQR, log.FC.dist and none."),
      make_option("--thres_outl", type = "double", default = 3.0,
                  help = "Outlier detection threshold"),
      make_option("--stand_method", type = "character", default = "std",
                  help = "Standardisation method. Currently support:
                          std, mad and custom."),
      make_option("--std_file", type = "character",
                  help = "user predifined std file with respect to ions"),
      make_option("--thres_symb", type = "double", default = 2.0,
                  help = "Symbolisation threshold"),

      #' Exploratory analysis
      make_option("--thres_ion_corr", type = "double", default = 0.15,
                  help = "Threshold for Ion correlation (0 - 1)"),

      #' Clustering analysis
      make_option("--min_clust_size", type = "double", default = 10.0,
                  help = "Minimal cluster size."),
      make_option("--h_tree", type = "double", default = 0.0,
                  help = "Cutting height for hierarchical clustering."),
      make_option("--filter_zero_string", type = "logical", default = TRUE,
                  help = "Filter the zero string or not"),

      #' Enrichment analysis
      make_option("--pval", type = "double", default = 0.05,
                  help = "P-values for enrichment analysis."),
      make_option("--min_count", type = "double", default = 3.0,
                  help = "Minimal count number for enrichment analysis."),
      make_option("--ont", type = "character", default = "BP",
                  help = "Ontology method: BP, MF and CC."),
      make_option("--annot_pkg", type = "character", default = "org.Sc.sgd.db",
                  help = "Annotation package"),

      #' Network analysis
      make_option("--method_corr", type = "character", default = "cosine",
                  help = "Similarity measure method. Currently support:
                          pearson, spearman, kendall, cosine, mahal_cosine,
                          hybrid_mahal_cosine"),
      make_option("--thres_corr", type = "double", default = 0.70,
                  help = "Similarity threshold for network analysis (0 - 1).
                          Features large than threshold will be kept."),

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
      make_option("--expl_anal_pdf",
        type = "character", default = "expl_anal.pdf",
        help = "Save plots from exploratory analysis"
      ),

      #' output: clustering analysis
      make_option("--clus_anal_pdf",
        type = "character", default = "clus_anal.pdf",
        help = "Save plots from clustering analysis"
      ),

      #' output: enrichment analysis
      make_option("--go_en_out",
        type = "character", default = "go_en.tsv",
        help = "Save GO enrichment table"
      ),

      #' output: network analysis
      make_option("--gene_net_pdf",
        type = "character", default = "gene_net.pdf",
        help = "Save plots from gene network"
      ),
      make_option("--imbe_out",
        type = "character", default = "impact_betweenness.tsv",
        help = "Save impact and betweenness table"
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
    ion_file = paste0(tool_dir, "test-data/Dataset_IonFlow_Ionome_KO_short.csv"),
    var_id = 1,
    batch_id = 3,
    data_id = 5,
    method_norm = "median",
    batch_control = "yes",
    control_lines = "BY4741",
    control_use = "all",
    method_outliers = "log.FC.dist",
    thres_outl = 3.0,
    stand_method = "std",
    thres_symb = 2,

    #' Exploratory analysis
    thres_ion_corr = 0.15,

    #' Clustering analysis
    min_clust_size = 10.0,
    h_tree = 0.0,
    filter_zero_string = TRUE,

    #' Enrichment analysis
    pval = 0.05,
    min_count = 3,
    ont = "BP",
    annot_pkg =  "org.Sc.sgd.db",

    #' Network analysis
    method_corr = "cosine",
    thres_corr = 0.7,

    #' output: pre-processing
    pre_proc_pdf       = paste0(tool_dir, "test-data/res/pre_proc.pdf"),
    df_stats_out       = paste0(tool_dir, "test-data/res/df_stats.tsv"),
    outl_out           = paste0(tool_dir, "test-data/res/outl.tsv"),
    data_wide_out      = paste0(tool_dir, "test-data/res/data_wide.tsv"),
    data_wide_symb_out = paste0(tool_dir, "test-data/res/data_wide_symb.tsv"),

    #' output: exploratory analysis
    expl_anal_pdf = paste0(tool_dir, "test-data/res/expl_anal.pdf"),

    #' output: clustering analysis
    clus_anal_pdf = paste0(tool_dir, "test-data/res/clus_anal.pdf"),

    #' output: enrichment analysis
    go_en_out = paste0(tool_dir, "test-data/res/go_en.tsv"),

    #' output: network analysis
    gene_net_pdf = paste0(tool_dir, "test-data/res/gene_net.pdf"),
    imbe_out     = paste0(tool_dir, "test-data/res/impact_betweenness.tsv")
  )
}
#' print(opt)

suppressPackageStartupMessages({
  source(paste0(tool_dir, "ionflow_funcs.R"))
})

## ==== Data preparation ====
ion_data <- read.table(opt$ion_file, header = T, sep = ",")

if (opt$batch_control == "yes") {
  control_lines <- opt$control_line
} else {
  control_lines <- NULL
}

if (opt$stand_method == "custom") { #' if (lenth(opt$std_file) > 0) {
  stdev <- read.table(opt$std_file, header = T, sep = "\t")
} else {
  stdev <- NULL
}

## ==== Pre-processing ====
pre <- PreProcessing(data            = ion_data,
                     var_id          = opt$var_id,
                     batch_id        = opt$batch_id,
                     data_id         = opt$data_id,
                     method_norm     = opt$method_norm,
                     control_lines   = control_lines,
                     control_use     = opt$control_use,
                     method_outliers = opt$method_outliers,
                     thres_outl      = opt$thres_outl,
                     stand_method    = opt$stand_method,
                     stdev           = stdev,
                     thres_symb      = opt$thres_symb)

#' save plot in pdf
pdf(file = opt$pre_proc_pdf, onefile = T) # width = 15, height = 10
plot(pre$plot.hist)
plot(pre$plot.overview)
plot(pre$plot.medians)
plot(pre$plot.CV)
plot(pre$plot.change.stat)
plot(pre$plot.change.dir)
dev.off()

#' combine stats
df_stats <- list(raw_data = pre$stats.raw.data,
                 bat_data = pre$stats.batches)
df_stats <- dplyr::bind_rows(df_stats, .id = "Data_Set")
row.names(df_stats) = NULL

#' save tables
write.table(df_stats, file = opt$df_stats_out, sep = "\t", row.names = F)
write.table(pre$stats.outliers, file = opt$outl_out, sep = "\t",
            row.names = F)
write.table(pre$data.line.zscores, file = opt$data_wide_out, sep = "\t",
            row.names = F)
write.table(pre$data.line.symb, file = opt$data_wide_symb_out,
            sep = "\t", row.names = F)

## ==== Exploratory analysis ====

pdf(file = opt$expl_anal_pdf, onefile = T)
expl <- IonAnalysis(data = pre$data.line.zscores,
                    thres_ion_corr = opt$thres_ion_corr)
plot(expl$plot.pca)
plot(expl$plot.net)
dev.off()

## ==== Clustering analysis ====

gcl <- ProfileClustering(pre$data.line.symb,
                         min_clust_size = opt$min_clust_size,
                         h_tree = opt$h_tree,
                         filter_zero_string = opt$filter_zero_string)

#' select larger clusters
cluster_vector <-
  gcl$clusters.vector[gcl$clusters.vector$Cluster %in%
                      gcl$tab.clusters.subset$Cluster, ]

#' extract symbolic and z-score prifiles for lines in selected clusters
symbol_profiles <- pre$data.line.symb
symbol_profiles$Cluster <-
  cluster_vector$Cluster[match(symbol_profiles$Line, cluster_vector$Line)]

zscore_profiles <- pre$data.line.zscores
zscore_profiles$Cluster <-
  cluster_vector$Cluster[match(zscore_profiles$Line, cluster_vector$Line)]

#' remove lines showing no phenotype
symbol_profiles <- symbol_profiles[!is.na(symbol_profiles$Cluster),]
zscore_profiles <- zscore_profiles[!is.na(zscore_profiles$Cluster),]

mat_long <- reshape2::melt(zscore_profiles, id = c("Line", "Cluster"),
                 variable.name = "Ion", value.name = "zscore")

mat_long$n.genes <-
  gcl$tab.clusters.subset$Number.of.genes[match(mat_long$Cluster,
                                                gcl$tab.clusters.subset$Cluster)]
mat_long$title <- paste0('Cluster ', mat_long$Cluster,' (',
                         mat_long$n.genes, ' genes)')

p_gcl <-
  ggplot(data = mat_long, aes(x = Ion, y = zscore, group = Line), color = "gray") +
  geom_line() +
  stat_summary(fun.data = "mean_se", color = "red", geom = "line", group = 1) +
  labs(x = "", y = "z-score") +
  coord_cartesian(ylim = c(-8, 8)) +
  facet_wrap(~title) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text = element_text(size = 10))

pdf(file = opt$clus_anal_pdf, onefile = T)
plot(p_gcl)
dev.off()

## ==== Enrichment analysis ====

ge <- GOEnricher(cluster_vector,
                 pval = opt$pval,
                 min_count = opt$min_count,
                 annot_pkg = opt$annot_pkg,
                 ont = opt$ont,
                 gene_uni = as.character(pre$data.line.zscores$Line))

if (nrow(ge$enrichment.summary) > 0) {
  write.table(ge$enrichment.summary, file = opt$go_en_out, sep = "\t",
              row.names = FALSE)
}

## ==== Network analysis ====

gn <- GeneticNetwork(data                 = zscore_profiles,
                     method_corr          = opt$method_corr,
                     thres_corr           = opt$thres_corr,
                     network_modules      = "input",
                     cluster_vector       = cluster_vector,
                     cluster_label_vector = NULL)

pdf(file = opt$gene_net_pdf, onefile = T) # width = 15, height = 10
plot(gn$plot.network)
plot(gn$plot.impact_betweenness)
dev.off()

write.table(gn$stats.impact_betweenness, file = opt$imbe_out,
            sep = "\t", row.names = FALSE)
