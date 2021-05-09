## ==== General settings ====
rm(list = ls(all = T))

tool_dir <- "~/ionflow-master/"

setwd(tool_dir)
pkgs <- c("optparse", "reshape2", "plyr", "dplyr", "tidyr", "ggplot2",
          "ggrepel", "corrplot", "gplots", "network", "sna", "GGally",
          "org.Sc.sgd.db","org.Hs.eg.db","GO.db", "GOstats", "KEGG.db",
          "pheatmap")

invisible(lapply(pkgs, library, character.only = TRUE))

source("ionflow_funcs.R")



## ==== S.Cerevisia Ionome ====
# Data preprocessing 
ion_data <- read.csv("./test-data/Dataset_IonFlow_Ionome_KO_short.csv", header = T)

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

## Exploratory Analysis Plot 
expl <- IonAnalysis(data = pre$data.line.zscores, thres_ion_corr = 0.15);

# PCA
expl$plot.pca

# Ion-ion correlation network
expl$plot.net

# Ions heatmap
expl$plot.heat

## Clustering 
gcl <- ProfileClustering(pre$data.line.symb, min_clust_size = 10, h_tree = 0, filter_zero_string = TRUE)

head(gcl$clusters.vector)
head(gcl$tab.clusters)
head(gcl$tab.clusters.subset)

# select larger clusters
cluster_vector <- gcl$clusters.vector[gcl$clusters.vector$Cluster %in% gcl$tab.clusters.subset$Cluster, ]

# extract symbolic and z-score prifiles for lines in selected clusters
symbol_profiles <- pre$data.line.symb
symbol_profiles$Cluster <- cluster_vector$Cluster[match(symbol_profiles$Line, cluster_vector$Line)]

zscore_profiles <- pre$data.line.zscores
zscore_profiles$Cluster <- cluster_vector$Cluster[match(zscore_profiles$Line, cluster_vector$Line)]

# remove lines showing no phenotype 
symbol_profiles <- symbol_profiles[!is.na(symbol_profiles$Cluster),]
zscore_profiles <- zscore_profiles[!is.na(zscore_profiles$Cluster),]

mat_long <- reshape2::melt(zscore_profiles,
                 id = c("Line", "Cluster"),
                 variable.name = "Ion",
                 value.name = "zscore"
                 )

mat_long$n.genes <- gcl$tab.clusters.subset$Number.of.genes[match(mat_long$Cluster, gcl$tab.clusters.subset$Cluster)]
mat_long$title <- paste0('Cluster ', mat_long$Cluster,' (', mat_long$n.genes, ' genes)')

ggplot(
    data = mat_long,
    aes(x = Ion, y = zscore, group = Line), 
    color = "gray") +
  facet_wrap(~title) +
  geom_line() +
  stat_summary(fun.data = "mean_se", color = "red", geom = "line", group=1) +
  labs(x = "", y = "z-score") +
  coord_cartesian(ylim = c(-8, 8)) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text = element_text(size = 10)
  )


## Enrichment Analysis 
# Input ORFs for yeast and ENTREZIDs for mouse or human

# BP
ge <- GOEnricher(cluster_vector, pval = 0.05, min_count = 3, 
                 annot_pkg = "org.Sc.sgd.db", ont = "BP",
                 gene_uni = as.character(pre$data.line.zscores$Line)) 
ge$enrichment.summary


## Network Analysis 
gn <- GeneticNetwork(data = zscore_profiles, 
               method_corr = "cosine",
               thres_corr = 0.7,
               network_modules = "input",
               cluster_vector = cluster_vector,
               cluster_label_vector = NULL);


# Plot Network 
gn$plot.network

# Plot Impact vs Betwennes
gn$plot.impact_betweenness



