#' wl-07-06-2021, Mon: Use jacopo's 'tutorial_galaxy_ionflow.R'

## ==== General settings ====

rm(list = ls(all = T))
tool_dir <- "~/my_galaxy/ionflow/"
setwd(tool_dir)
pkgs <- c("optparse", "reshape2", "plyr", "dplyr", "tidyr", "ggplot2",
          "ggrepel", "corrplot", "gplots", "network", "sna", "GGally",
          "org.Sc.sgd.db","org.Hs.eg.db","GO.db", "GOstats", "KEGG.db",
          "pheatmap")
invisible(lapply(pkgs, library, character.only = TRUE))
source("ionflow_funcs.R")

## ==== S.Cerevisia Ionome ====

ion_data <- read.table("./test-data/Dataset_IonFlow_Ionome_KO_short.csv", 
                       header = T, sep = ",")

## ==== Data preprocessing ====

pre <- PreProcessing(data = ion_data,
                     var_id = 1, batch_id = 3, data_id = 5,
                     method_norm = "median",
                     control_lines = "BY4741",
                     control_use = "all",
                     method_outliers = "log.FC.dist",
                     thres_outl = 3,
                     stand_method = "std",
                     stdev = NULL,
                     thres_symb = 2)
## names(pre)
##  [1] "stats.raw.data"    "stats.outliers"    "stats.batches"
##  [4] "stats.std"         "data.long"         "data.line.logFC"
##  [7] "data.line.zscores" "data.line.symb"    "plot.overview"
## [10] "plot.hist"         "plot.change.stat"  "plot.change.dir"
## [13] "plot.medians"      "plot.CV

pre$stats.raw.data
pre$stats.outliers
pre$stats.batches
pre$stats.std
head(pre$data.line.zscores)
head(pre$data.line.symb)
head(pre$data.long)

# z-score distributions 
pre$plot.hist
# processed samples overview 
pre$plot.overview 
# batch median log-concentrations 
pre$plot.medians  
# CV batch median log-concentrations
pre$plot.CV  
# Histogram Number of Elements Changed
pre$plot.change.stat  
# Histogram Number of Changes per Element
pre$plot.change.dir  

## ==== Exploratory Analysis Plot ====

expl <- IonAnalysis(data = pre$data.line.zscores, thres_ion_corr = 0.15)
## names(expl)
## [1] "data.pca.load" "plot.pca" "plot.corr" "plot.net" "plot.heat"

# PCA
expl$plot.pca
# Ion-ion correlation network
expl$plot.net
# Ions heatmap
expl$plot.heat

## ==== Clustering ====

gcl <- ProfileClustering(pre$data.line.symb, min_clust_size = 10,
                         h_tree = 0, filter_zero_string = TRUE)
## names(gcl)
## [1] "clusters.vector" "tab.clusters" "tab.clusters.subset"
          
head(gcl$clusters.vector)
head(gcl$tab.clusters)
head(gcl$tab.clusters.subset)

# select larger clusters
cluster_vector <- 
  gcl$clusters.vector[gcl$clusters.vector$Cluster %in%
                      gcl$tab.clusters.subset$Cluster, ]

# extract symbolic and z-score prifiles for lines in selected clusters
symbol_profiles <- pre$data.line.symb
symbol_profiles$Cluster <- 
  cluster_vector$Cluster[match(symbol_profiles$Line, cluster_vector$Line)]

zscore_profiles <- pre$data.line.zscores
zscore_profiles$Cluster <- 
  cluster_vector$Cluster[match(zscore_profiles$Line, cluster_vector$Line)]

# remove lines showing no phenotype 
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
  stat_summary(fun.data = "mean_se", color = "red", geom = "line", group=1) +
  labs(x = "", y = "z-score") +
  coord_cartesian(ylim = c(-8, 8)) +
  facet_wrap(~title) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text = element_text(size = 10))
p_gcl

## ==== Enrichment Analysis ====

# Input ORFs for yeast and ENTREZIDs for mouse or human

# BP
ge <- GOEnricher(cluster_vector, pval = 0.05, min_count = 3, 
                 annot_pkg = "org.Sc.sgd.db", ont = "BP",
                 gene_uni = as.character(pre$data.line.zscores$Line)) 
## names(ge)
## [1] "enrichment.summary" "enrichment.full.results"
ge$enrichment.summary

## ==== Network Analysis ====

gn <- GeneticNetwork(data = zscore_profiles, 
                     method_corr = "cosine",
                     thres_corr = 0.7,
                     network_modules = "input",
                     cluster_vector = cluster_vector,
                     cluster_label_vector = NULL)
## names(gn)
## [1] "network" "network.modules"
## [3] "plot.network" "plot.impact_betweenness"
## [5] "stats.impact_betweenness"

# Plot Network 
gn$plot.network
# Plot Impact vs Betwennes
gn$plot.impact_betweenness
