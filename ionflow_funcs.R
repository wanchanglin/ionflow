#' =======================================================================
#'
PreProcessing <- function(data = NULL, var_id = 1, batch_id = 3, data_id = 5,
                          method_norm = c("median", "median+std", "none"),
                          control_lines = NULL,
                          control_use = c("control", "all", "control.out"),
                          method_outliers = c("mad", "IQR", "log.FC.dist", "none"),
                          thres_outl = 3,
                          stand_method = c("std", "mad", "custom"),
                          stdev = NULL, thres_symb = 2) {


  #' -------------------> Import data
  #' ji: get sample id
  data$sample_id <- rownames(data)

  data <- data[, c(var_id, ncol(data), batch_id, data_id:(ncol(data) - 1))]
  names(data)[1:3] <- c("Line", "Sample_ID", "Batch_ID")
  mat <- data[, -c(1:3)]

  #' ji: Remove samples with zeros or negative or missing values
  data <- data[!is.na(sum(mat)), ]
  data <- data[!(sum(mat <= 0) > 0), ]

  #' get summary stats
  res <- as.data.frame(t(sapply(mat, function(x) {
    c(round(summary(x), 3), round(var(x), 3))
  })))
  names(res)[ncol(res)] <- "Variance"
  res <- cbind(Ion = names(mat), res)
  rownames(res) <- NULL
  df_raw <- res

  #' data long format
  data_long <- reshape2::melt(data,
    id = c("Line", "Sample_ID", "Batch_ID"), variable.name = "Ion",
    value.name = "Concentration"
  )

  #' convert to factors before using levels function.
  data_long$Line <- factor(data_long$Line)
  data_long$Ion <- factor(data_long$Ion)
  ion_name <- levels(data_long$Ion)

  #' -------------------> Median batch correction
  data_long$control <- rep(1, length(data_long$Concentration))
  data_long$log <- log(data_long$Concentration)

  #' ji: control use in batch correction
  if (length(control_lines) > 0) {
    if (control_use == "all") { }
    if (control_use == "control") {
      data_long$control[!(data_long$Line %in% control_lines)] <- 0
    }
    if (control_use == "control.out") {
      data_long$control[data_long$Line %in% control_lines] <- 0
    }
  }

  #' batch correction methods
  if (method_norm == "median") {
    data_long <- plyr::ddply(data_long, "Ion", function(x) {
      res <- plyr::ddply(x, "Batch_ID", function(y) {
        med <- median(y$log[y$control == 1])
        y$log_corr <- y$log - med
        y$med <- med
        y
      })
    })
  }
  if (method_norm == "median+std") {
    data_long <- plyr::ddply(data_long, "Ion", function(x) {
      res <- plyr::ddply(x, "Batch_ID", function(y) {
        med <- median(y$log[y$control == 1])
        st_de <- sd(y$log[y$control == 1])
        y$log_corr <- y$log - med - st_de
        y$med <- med
        y
      })
    })
  }
  if (method_norm == "none") {
    data_long$log_corr <- data_long$log
    data_long$med <- NULL
  }

  #' get correction stats
  res <- plyr::ddply(data_long, "Ion", function(x) {
    c(round(summary(x$log_corr), 3), round(var(x$log_corr), 3))
  })
  names(res)[ncol(res)] <- "Variance"
  df_bat <- res

  #' -------------------> Outlier detection

  #' ji: outlier detection methods
  if (method_outliers == "IQR") {
    data_long <- plyr::ddply(data_long, "Ion", function(x) {
      res <- plyr::ddply(x, "Batch_ID", function(y) {
        lowerq <- quantile(y$log_corr, na.rm = T)[2]
        upperq <- quantile(y$log_corr, na.rm = T)[4]
        iqr <- upperq - lowerq
        extreme.t.upper <- (iqr * thres_outl) + upperq
        extreme.t.lower <- lowerq - (iqr * thres_outl)
        y$Outlier <- ifelse((y$log_corr > extreme.t.upper) |
          (y$log_corr < extreme.t.lower), 1, 0)
        return(y)
      })
    })
  }
  if (method_outliers == "mad") {
    data_long <- plyr::ddply(data_long, "Ion", function(x) {
      res <- plyr::ddply(x, "Batch_ID", function(y) {
        med_dev <- mad(y$log_corr)
        extreme.t.upper <- (med_dev * thres_outl)
        extreme.t.lower <- -(med_dev * thres_outl)
        y$Outlier <- ifelse((y$log_corr > extreme.t.upper) |
          (y$log_corr < extreme.t.lower), 1, 0)
        return(y)
      })
    })
  }
  if (method_outliers == "log.FC.dist") {
    data_long <- plyr::ddply(data_long, "Ion", function(x) {
      res <- plyr::ddply(x, "Batch_ID", function(y) {
        extreme.t.upper <- thres_outl
        extreme.t.lower <- -thres_outl
        y$Outlier <- ifelse((y$log_corr > extreme.t.upper) |
          (y$log_corr < extreme.t.lower), 1, 0)
        return(y)
      })
    })
  }

  if (method_outliers == "none") {
    data_long$Outlier <- rep(0, length(data_long$Line))
    df_outlier <- data.frame()
  } else {
    df_outlier <- data.frame(cbind(
      levels(data_long$Ion),
      table(data_long$Ion, data_long$Outlier),
      round(table(data_long$Ion, data_long$Outlier)[, 2] /
            dim(data_long)[1] * 100, 2)
    ))
    rownames(df_outlier) <- c()
    colnames(df_outlier) <- c("Ion", "no_outlier", "outlier", "outlier(%)")
  }

  #' ji: remove samples with at least an outlier value
  samples_to_exclude <- unique(data_long$Sample_ID[data_long$Outlier == 1])
  data_long <- data_long[!(data_long$Sample_ID %in% samples_to_exclude), ]

  #' -------------------> Standardisation

  #' ji: standardisation methods
  if (stand_method == "std") {
    sds <- plyr::ddply(data_long, "Ion", function(x) sd(x$log_corr))
    nam <- sds[, 1]
    sds <- as.numeric(as.vector(sds[, 2]))
    names(sds) <- nam
  }
  if (stand_method == "mad") {
    sds <- plyr::ddply(data_long, "Ion", function(x) mad(x$log_corr))
    nam <- sds[, 1]
    sds <- as.numeric(as.vector(sds[, 2]))
    names(sds) <- nam
  }
  if (stand_method == "custom") {
    # specific 2-columns format for vector of sd
    sds <- stdev
    nam <- sds[, 1]
    sds <- as.numeric(as.vector(sds[, 2]))
    names(sds) <- nam
  }

  #' ji: aggregate measurements at line level (median)
  data_wide_line_log_norm <- reshape2::dcast(data_long, Line ~ Ion,
    fun.aggregate = median,
    value.var = "log_corr"
  )

  #' ji: normalise by stds
  data_wide_line_z_score <- data_wide_line_log_norm
  data_wide_line_z_score[, 2:ncol(data_wide_line_z_score)] <- data_wide_line_z_score[, 2:ncol(data_wide_line_z_score)] / sds

  #' -------------------> Symbolisation
  symb_profiles <- data_wide_line_z_score[, 2:ncol(data_wide_line_z_score)]
  symb_profiles[(symb_profiles > -thres_symb) & (symb_profiles < thres_symb)] <- 0
  symb_profiles[symb_profiles >= thres_symb] <- 1
  symb_profiles[symb_profiles <= -thres_symb] <- -1

  #' ji: save sybolic profiles
  data_wide_line_symb <- cbind(Line = data_wide_line_z_score$Line, symb_profiles)

  #' plot z-score distributions
  p1 <-
    ggplot(
      data = data_long,
      aes(x = factor(Batch_ID), y = log_corr, col = factor(Batch_ID))
    ) +
    geom_point(shape = 1) +
    facet_wrap(~Ion) +
    xlab("Batch.ID") +
    ylab("Log FC to batch median") +
    theme(legend.position = "none") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  # plot overview processed samples
  dat <- reshape2::melt(data_wide_line_z_score, id = "Line")
  p2 <-
    ggplot(data = dat, aes(x = value)) +
    geom_histogram(binwidth = .1) +
    facet_wrap(~variable) +
    xlab("Concentration (z-score)") +
    ylab("Counts") +
    xlim(-10, 10) +
    geom_vline(xintercept = c(-thres_symb, thres_symb), col = "red")

  if (method_norm != "none"){
  # plot batch median log-concentrations
  plot.data.long <- data_long[data_long["Outlier"]==0, ] 
  batch_median_ions <- plot.data.long[!duplicated(plot.data.long[,"med"]), c("Batch_ID","Ion","med"), drop = FALSE];
  
  ion.batch.mean.median <- batch_median_ions %>%
    dplyr::group_by(Ion) %>%
    dplyr::summarize(mean = mean(med, na.rm = TRUE)) 
  
  ion.order <- as.character(ion.batch.mean.median$Ion[sort.list(ion.batch.mean.median$mean)],  decreasing = TRUE)
  
  batch_median_ions$Ion <- factor(batch_median_ions$Ion, levels = ion.order)
  
 p3 <- ggplot(batch_median_ions, aes(x=Ion, y=med, group = Batch_ID)) +
    geom_line(aes(), color = "gray") +
    theme(text = element_text(size=20), axis.title.x=element_blank()) +
    ylab("Log-concentration") +
    stat_summary(aes(y = med, group=1), fun = mean, colour="green1", geom="line", group=1) + 
    stat_summary(fun.min = function(x) mean(x) - sd(x), fun.max = function(x) mean(x) + sd(x), position ='dodge', 
                 colour="green1", geom = 'errorbar', aes(group = 1), width=.1)+
    ggtitle("Median log-concentrations of batches");
  
  # plot CV batch median log-concentrations
  ion.batch.cv.median <- batch_median_ions %>%
    dplyr::group_by(Ion) %>%
    dplyr::summarize(CV = sd(med, na.rm = TRUE)/mean(med, na.rm = TRUE)) 
  
  ion.batch.cv.median <- data.frame(ion.batch.cv.median)
  
  ion.batch.cv.median <- ion.batch.cv.median[order(ion.batch.cv.median$CV),]
  
  ion.batch.cv.median$Ion <- factor(ion.batch.cv.median$Ion, levels = ion.batch.cv.median$Ion)
  
  p4 <- ggplot(ion.batch.cv.median, aes(x=Ion, y=abs(CV), group =1)) +
    ylab(" Absolute CV") +
    geom_line(aes(), color = "red") +
    theme(text = element_text(size=20), axis.title.x=element_blank()) +
    scale_y_log10() +
    ggtitle("CV median log-concentrations across batches");
  
  }
  
  #' -------------------> Output
  res                    <- list()
  res$stats.raw.data     <- df_raw                    # raw data
  res$stats.outliers     <- df_outlier                # outliers
  res$stats.batches      <- df_bat                    # batch corrected data
  res$stats.std          <- sds                       # standard deviations
  res$data.long          <- data_long                 # with Batch_ID
  res$data.line.logFC    <- data_wide_line_log_norm
  res$data.line.zscores  <- data_wide_line_z_score
  res$data.line.symb     <- data_wide_line_symb
  res$plot.overview      <- p1
  res$plot.hist          <- p2
  res$plot.medians       <- p3
  res$plot.CV            <- p4
  return(res)
}

#' =======================================================================
#'
IonAnalysis <- function(data = NULL, thres_ion_corr = 0.15, method_ion_corr = "pearson"){
  
  #' -------------------> PCA
  ionProfile.PCA = NULL
  ionProfile.PCA$pr.y = prcomp(data[, -1],scale. = F)
  ionProfile.PCA$y = cbind(data.frame(Line = data$Line), as.data.frame(ionProfile.PCA$pr.y$x))
  
  unit.norm = function(x)(x / sqrt(sum(x^2)))
  ionProfile.PCA$y[,grepl("^PC",colnames(ionProfile.PCA$y))] = apply(ionProfile.PCA$y[,grepl("^PC",colnames(ionProfile.PCA$y))],2,unit.norm)
  
  p.PCA.data = ionProfile.PCA$y %>%
    dplyr::select(Line, PC1, PC2 )
  
  my.annotation = tbl_df(ionProfile.PCA$pr.y$rotation) %>%
    dplyr::select(PC1, PC2) %>%
    dplyr::mutate(Ion = row.names(ionProfile.PCA$pr.y$rotation)) %>%
    `colnames<-`(c("x", "y", "Ion")) 
  
  p_pca <- p.PCA.data %>%
    ggplot(aes(x = PC1, y = PC2)) +
    geom_segment(inherit.aes = F, data = my.annotation,
                 aes(x = 0, y = 0, xend = x/2, yend = y/2),
                 color = "blue",
                 arrow = arrow(length = unit(0.02, "npc")))+
    geom_text(inherit.aes = F, data = my.annotation,
              aes(x = 0.51*x, y = 0.51*y, label = Ion),
              color = "red", size = 4) + 
    geom_point(size = 1, alpha = 1/20) + 
    theme(aspect.ratio = 1) +
    scale_shape(solid = FALSE) +
    labs(x = paste0("PC1",
                    " (",round(summary(ionProfile.PCA$pr.y)$importance[2,"PC1"],2)*100,"%)"),
         y = paste0("PC2"," (",
                    round(summary(ionProfile.PCA$pr.y)$importance[2,"PC2"],2)*100,"%)"),
         color = "changed")  
  
  #' PCA loading
  pca_loadings <- data.frame(ionProfile.PCA$pr.y$rotation)
  rownames(pca_loadings) <- colnames(data[, -1])
  pca_loadings <- pca_loadings[, 1:2]
  
  
  #' -------------------> Correlation
  col3 <- colorRampPalette(c("steelblue4", "white", "firebrick"))
  corrplot.mixed(cor(data[, -1], use = "complete.obs", method = method_ion_corr),
    number.cex = .7,
    lower.col = "black", upper.col = col3(100)
  )
  p_corr <- recordPlot()


  #' -------------------> Heatmap with dendrogram
  pheatmap(data[, -1], show_rownames = F, cluster_cols = T, cluster_rows = T,
           legend = T, fontsize = 15, clustering_method = "ward.D",
           scale = "row")
  p_heat <- recordPlot()

  #' -------------------> Correlation network
  #' wl-13-08-2020, Thu: there is no 'qgraph' in conda forge and bio conda.
  if (T) {
    #' library(glasso)
    #' corr_reg <- glasso(corr, rho = 0.01)
    #' net <- network::network(corr_reg$w, directed = FALSE)
    corr <- cor(na.omit(data[, -1]), method = method_ion_corr)
    corr[abs(corr) < thres_ion_corr] <- 0
    net <- network::network(corr[], directed = FALSE)

    #' set edge attributes
    net %e% "weight" <- corr
    net %e% "weight_abs" <- abs(corr) * 6
    net %e% "color" <- ifelse(net %e% "weight" > 0, "lightgreen", "coral")
    p_net <-
      ggnet2(net,
        label = TRUE, mode = "fruchtermanreingold",
        node.size = 10, edge.size = "weight_abs", edge.color = "color"
      )
  } else {
    #' wl-06-07-2020, Mon: 'cor_auto' is from package qgraph(lavaan)
    #' wl-28-07-2020, Tue: cad and corr are the same
    cad <- cor_auto(data[, -1])
    suppressWarnings(qgraph(cad,
      graph = "glasso", layout = "spring",
      tuning = 0.25, sampleSize = nrow(data[, -1])
    ))
    Graph_lasso <- recordPlot()
  }

  res <- list()
  res$data.pca.load <- pca_loadings
  res$plot.pca  <- p_pca
  res$plot.corr <- p_corr
  res$plot.net <- p_net
  res$plot.heat <- p_heat
  return(res)

}



#' =======================================================================
#' wl-04-10-2020, Sun: Hierarchical clustering
#' ji-10-01-2020, Sun: Input, Method, Output simplified

ProfileClustering <- function(symbol_profiles, min_clust_size = 10, h_tree = 0, filter_zero_string = TRUE) {
  
  if (filter_zero_string){
    symbol_profiles <- symbol_profiles[rowSums(symbol_profiles[, -1])!=0, ]  
  }
  x <- symbol_profiles[, -1] 
  dis <- stats::dist(x, method = "manhattan")
  hc <- hclust(d = dis, method = "single")
  clus <- data.frame(Line = symbol_profiles$Line, Cluster = cutree(hc, h = h_tree))
  
  tab <- as.data.frame(table(clus$Cluster), stringsAsFactors = F)
  names(tab) <- c("Cluster", "Number.of.genes")
  
  tab_subset <- tab[tab$Number.of.genes > min_clust_size, ]
  tab_subset <- tab_subset[order(tab_subset$Number.of.genes, decreasing = T), ]
  idx_subset <- clus %in% tab_subset$Cluster
  
  
  res <- list()
  res$clusters.vector <- clus
  res$tab.clusters <- tab
  res$tab.clusters.subset <- tab_subset

  return(res)
  
}


#' =======================================================================
#' wl-06-11-2020, Fri: Get ENTREZID  from SYMBOL
#'
get_entrez_id <- function(gene_list, annot_pkg = "org.Sc.sgd.db", key_type = "ORF") {
  
  res <- AnnotationDbi::select(get(annot_pkg), keys = gene_list,
                               columns = "ENTREZID", keytype = key_type)
  res <- res[,2,drop = T]
  res <- res[!is.na(res)]
  res <- res[!duplicated(res)]

  return(res)
}

#' =======================================================================
#' wl-03-10-2020, Sat: KEGG enrichment analysis for symbolization data
#' wl-06-11-2020, Fri: The first column of data must be ORF for
#'  org.Sc.sgd.db or SYMBOL for any other annotation packages.
#' ji-11-01-2020, Mon: Input, Method, Output simplified
#' 

KeggEnricher <- function(cluster_vector, pval = 0.05, 
                        min_count = 3, annot_pkg = "org.Sc.sgd.db", 
                        gene_uni = NULL) {

  if(is.null(gene_uni)){
    gene_uni <- as.character(cluster_vector$Line)
  }
  
  clusters_set <- unique(cluster_vector$Cluster)
  
  #' geneIds can be ORF or ENTREZID
  enrich <- lapply(clusters_set, function(x) {
    params <- new("KEGGHyperGParams",
                  geneIds = gene_uni[cluster_vector$Cluster == x],
                  universeGeneIds = gene_uni,
                  annotation = annot_pkg,
                  categoryName = "KEGG",
                  pvalueCutoff = pval,
                  testDirection = "over")

    over <- hyperGTest(params)
  })

  # name list with original clusters ID 
  names(enrich) <- clusters_set
  
  #' There is no explicit methods for getting manual summary table.
  summ <- lapply(enrich, function(x) {
    tmp <- summary(x)
    if (nrow(tmp) == 0) tmp <- NULL #' wl-04-10-2020, Sun: it happens very often.
    return(tmp)
  })
  summ <- summ[!sapply(summ,is.null)]

  #' binding and filtering
  summ <- lapply(summ, "[", -c(3, 4)) %>%
    dplyr::bind_rows(.id = "Cluster") %>%
    dplyr::filter(Pvalue <= pval & Count >= min_count)

  res <- list()
  res$enrichment.summary <- summ
  res$enrichment.full.results <- enrich
  return(res)
  
  
}

#' =======================================================================
#' wl-03-10-2020, Sat: GO enrichment analysis for symbolization data
#' wl-06-11-2020, Fri: The first column of data must be ORF for
#'  org.Sc.sgd.db or SYMBOL for any other annotation packages.
#'
GOEnricher <- function(cluster_vector, pval = 0.05, ont = c("BP","CC","MF"),
                      min_count = 3, annot_pkg = "org.Sc.sgd.db",
                      gene_uni = NULL) {
  
  if(is.null(gene_uni)){
    gene_uni <- as.character(cluster_vector$Line)
  }
  
  ont <- match.arg(ont)
  
  genes_not_annotated <- NULL
  
  if ( annot_pkg == "org.Sc.sgd.db"){
    genes_not_annotated <- unlist(mget(
      mappedkeys(org.Sc.sgdGO2ALLORFS)[which(is.na(Term(mappedkeys(org.Sc.sgdGO2ALLORFS))))], org.Sc.sgdGO2ALLORFS)
      )
    gene_uni <- mappedkeys(org.Sc.sgdGO)[mappedkeys(org.Sc.sgdGO) %in% gene_uni]
  }
  if ( annot_pkg == "org.Hs.eg.db"){
    genes_not_annotated <- unlist(mget(
      mappedkeys(org.Hs.egGO2ALLEGS)[which(is.na(Term(mappedkeys(org.Hs.egGO2ALLEGS))))], org.Hs.egGO2ALLEGS)
    )
    #gene_uni <- mappedkeys(org.Hs.egGO)
    gene_uni <- mappedkeys(org.Hs.egGO)[mappedkeys(org.Hs.egGO) %in% gene_uni]
  }
  if ( annot_pkg == "org.Mm.eg.db"){
    genes_not_annotated <- unlist(mget(
      mappedkeys(org.Mm.egGO2ALLEGS)[which(is.na(Term(mappedkeys(org.Mm.egGO2ALLEGS))))], org.Mm.egGO2ALLEGS)
    )
    gene_uni <- mappedkeys(org.Mm.egGO)[mappedkeys(org.Mm.egGO) %in% cluster_vector$Line]
  }  
  
  gene_uni <- as.character(gene_uni[!gene_uni %in% genes_not_annotated])
    
  cluster_vector <- cluster_vector[cluster_vector$Line %in% gene_uni,]
  clusters_set <- unique(cluster_vector$Cluster)
  
  #' geneIds can be ORF or ENTREZID
  enrich <- lapply(clusters_set, function(x) { 
    gene_set <- as.character(cluster_vector$Line[cluster_vector$Cluster == x])
    params <- new("GOHyperGParams",
                  geneIds = gene_set,
                  universeGeneIds = gene_uni,
                  annotation = annot_pkg,
                  categoryName = "GO",
                  ontology = ont,
                  pvalueCutoff = pval,
                  conditional = T,
                  testDirection = "over")

    over <- hyperGTest(params)
  })

  # name list with original clusters ID 
  names(enrich) <- clusters_set
  
  summ <- lapply(enrich, function(x) {
    Pvalue        <- round(pvalues(x), digit = 4)
    ID            <- names(Pvalue)
    Description   <- Term(ID)
    OddsRatio     <- oddsRatios(x)
    ExpCount      <- expectedCounts(x)
    Count         <- geneCounts(x)
    CountUniverse <- universeCounts(x)

    tab <- cbind(ID, Description, Pvalue, OddsRatio, ExpCount, Count, CountUniverse)
    rownames(tab) <- NULL
    tab <- na.omit(tab)   #' wl-03-10-2020, Sat: this is why summary fails.
    tab <- data.frame(tab, Ontology = ont)
  })

  summ <- lapply(summ, "[", -c(4, 5)) %>%
            dplyr::bind_rows(.id = "Cluster") %>%
            dplyr::filter(Pvalue <= pval & Count >= min_count)

  res <- list()
  res$enrichment.summary <- summ
  res$enrichment.full.results <- enrich
  return(res)
}




GeneticNetwork <- function(data = NULL, 
                        method_corr = c("pearson", "spearman", "kendall",
                          "cosine", "mahal_cosine", "hybrid_mahal_cosine"),
                        network_modules = c("louvain", "input"),
                        thres_corr = 0.7,
                        cluster_vector = NULL,
                        cluster_label_vector = NULL,
                        n_labels = 3
                        ) {

  
  mat <- as.matrix(data[,!colnames(data)%in%c("Line","Cluster")])
  

  #' Compute similarity matrix
  if (method_corr == "pearson" ||
      method_corr == "spearman" ||
      method_corr == "kendall") {
    corrGenes <- cor(t(as.matrix(mat)), 
                     method = method_corr,
                     use = "pairwise.complete.obs")
  } else if (method_corr == "cosine") {
    corrGenes <- cosine(t(as.matrix(mat)))
  } else if (method_corr == "mahal_cosine") {
    corrGenes <- cosM(mat, mode = "normal")
  } else if (method_corr == "hybrid_mahal_cosine") {
    corrGenes <- cosM(mat, mode = "hybrid")
  }

  #' Adjacency matrix
  A <- corrGenes
  diag(A) <- 0
  A <- (A > thres_corr)
  A <- ifelse(A == TRUE, 1, 0)
  
  if (network_modules == "input"){
   
    if (length(cluster_label_vector) == 0){
      cluster_vector <- as.character(cluster_vector$Cluster)
    }else{
      mm <- match(cluster_vector$Cluster, cluster_label_vector$Cluster)
      cluster_vector <- cluster_label_vector$label[mm]
    }
    
    tmp <- unique(cluster_vector)
    if (length(tmp) != 1) {
      cpy <- rainbow(length(tmp))
      names(cpy) <- tmp
    } else {
      cpy <- "Set2"
    }
    
  }
  if (network_modules == "louvain"){
    # community detection
    tmp <- igraph::graph_from_adjacency_matrix(A, mode="undirected")
    com <- igraph::cluster_louvain(tmp, weights = NULL)
    cluster_vector <- as.character(igraph::membership(com))
    degree_vector <- igraph::degree(tmp)
  
    tmp <- unique(cluster_vector)
    if (length(tmp) != 1) {
      cpy <- rainbow(length(tmp))
      names(cpy) <- tmp
    } else {
      cpy <- "Set2"
    }
  }
  
  
  #' Generate network
  net <- network::network(A, directed = FALSE)
  net %v% "Label" <- cluster_vector
  
  #' Remove communities of size 1
  if (network_modules == "louvain"){
    net <- delete.vertices(net, which(degree_vector == 0));
  }
  
  net_p <- GGally::ggnet2(net,
                           mode = "fruchtermanreingold",
                           color = "Label",
                           palette = cpy,
                           edge.alpha = 0.5, size = 2, 
                           color.legend = "Modules",
                           legend.size = 10, 
                           legend.position = "right"
  )
  
  #' network edge list 
  net <- as.matrix(net, matrix.type = "edgelist")[,1:2];
  net[,] <- cbind(as.character(data$Line[net[,1]]), as.character(data$Line[net[,2]]))

  #' Impact and betweenness
  btw <- sna::betweenness(A) # or use 'net' instead of 'A'
  impact <- apply(mat, 1, norm, type = "2") # L2 norm

  df.res <- data.frame(
    Line = data$Line,
    impact = round(impact, 3),
    betweenness = round(btw, 3),
    log.betweenness = round(log(btw + 1), 3),
    pos = factor(ifelse((impact < quantile(impact, .75)) & (log(btw + 1) < quantile(log(btw + 1), .75)), 1,
      ifelse((impact < quantile(impact, .75)) & (log(btw + 1) > quantile(log(btw + 1), .75)), 2,
        ifelse((impact > quantile(impact, .75)) & (log(btw + 1) < quantile(log(btw + 1), .75)), 3, 4)
      )
    )),
    pos.label = factor(ifelse((impact < quantile(impact, .75)) & (log(btw + 1) < quantile(log(btw + 1), .75)), "Low impact, low betweenness",
      ifelse((impact < quantile(impact, .75)) & (log(btw + 1) > quantile(log(btw + 1), .75)), "Low impact, high betweenness",
        ifelse((impact > quantile(impact, .75)) & (log(btw + 1) < quantile(log(btw + 1), .75)), "High impact, low betweenness", "High impact, high betweenness")
      )
    ))
  )
  rownames(df.res) <- data$Line[]

  q1 <- row.names(subset(df.res, (impact < quantile(impact, .75)) & (log.betweenness < quantile(log.betweenness, .75))))
  q2 <- row.names(subset(df.res, (impact < quantile(impact, .75)) & (log.betweenness > quantile(log.betweenness, .75))))
  q3 <- row.names(subset(df.res, (impact > quantile(impact, .75)) & (log.betweenness < quantile(log.betweenness, .75))))
  q4 <- row.names(subset(df.res, (impact > quantile(impact, .75)) & (log.betweenness > quantile(log.betweenness, .75))))

  #' labels 
  N <- n_labels
  lst <- list(q1, q2, q3, q4) 
  idx <- lapply(lst, function(x) {
    if (length(x) > N){
      tmp <- df.res[x, c("betweenness","impact")];
      union(rownames(tmp)[order(tmp$betweenness, decreasing = TRUE)[1:N]],
            rownames(tmp)[order(tmp$impact, decreasing = TRUE)[1:N]])
    }else{x}
  })
  idx <- unique(unlist(idx))

  df.idx <- df.res[idx, ]

  im_be_p <-
    ggplot(data = df.res, aes(x = impact, y = log.betweenness)) +
    geom_point(aes(col = pos.label), alpha = .3, size = 3) +
    scale_color_manual(values = c(
      "plum4", "palegreen4", "indianred",
      "cornflowerblue"
    )) +
    theme_linedraw() +
    theme_light() +
    theme(legend.position = "bottom") +
    guides(colour = guide_legend(nrow = 2)) +
    theme(legend.title = element_blank()) +
    geom_text_repel(data = df.idx, aes(label = Line), size = 3.5) +
    geom_vline(xintercept = quantile(df.res$impact, .75), linetype = "dashed") +
    geom_hline(
      yintercept = quantile(df.res$log.betweenness, .75),
      linetype = "dashed"
    ) +
    xlab("Impact") +
    ylab("Log(betweenness+1)")

  
  res <- list()
  res$network <- net 
  res$network.modules <- cluster_vector
  res$plot.network <- net_p # plot gene network with symbolic pheno
  res$plot.impact_betweenness <- im_be_p # plot impact betweenees
  res$stats.impact_betweenness <- df.res # impact and betweenees data
  return(res)
}





cosM <- function(x, mode = c("normal", "hybrid")) {

  #' --------------------------------------------------------------------
  #' library("pracma")
  #' pinv: Pseudoinverse (Moore-Penrose Generalized Inverse)
  pinv <- function (A, tol = .Machine$double.eps^(2/3)) {
    stopifnot(is.numeric(A), length(dim(A)) == 2, is.matrix(A))
    s <- svd(A)
    p <- (s$d > max(tol * s$d[1], 0))
    if (all(p)) {
      mp <- s$v %*% (1/s$d * t(s$u))
    } else if (any(p)) {
      mp <- s$v[, p, drop=FALSE] %*% (1/s$d[p] * t(s$u[, p, drop=FALSE]))
    } else {
      mp <- matrix(0, nrow=ncol(A), ncol=nrow(A))
    }
    return(mp)
  }

  mldivide <- function(A, B, pinv = TRUE) {
    stopifnot( is.numeric(A) || is.complex(A), is.numeric(B) || is.complex(B))
    if (is.vector(A)) A <- as.matrix(A)
    if (is.vector(B)) B <- as.matrix(B)
    if (nrow(A) != nrow(B)) {
      stop("Matrices 'A' and 'B' must have the same number of rows.")
    }
    if (pinv) {
      pinv(t(A) %*% A) %*% t(A) %*% B
    } else {
      qr.solve(A, B)
    }
  }

  mrdivide <- function(A, B, pinv = TRUE) {
    stopifnot( is.numeric(A) || is.complex(A), is.numeric(B) || is.complex(B))
    if (is.vector(A)) A <- t(A)
    if (is.vector(B)) B <- t(B)
    if (ncol(A) != ncol(B)) {
      stop("Matrices 'A' and 'B' must have the same number of columns.")
    }
    t(mldivide(t(B), t(A), pinv = pinv))
  }

  #' --------------------------------------------------------------------

  #' compute covariance
  Cov <- cov(x, use = "pairwise.complete.obs")
  n <- dim(x)[1]
  m <- dim(x)[2]

  #' compute eigenvalues
  Eig <- eigen(Cov)
  score <- mrdivide(as.matrix(x), t(Eig$vectors))

  #' wl-25-10-2020, Sun: n must be larger than m othwise this script fails.
  if (any(Eig$values < 0)) {
    stop("Cov Not Positive SemiDefinite !")
  }

  #' compute pairwise cosine similarity
  C <- rep(0, n * (n - 1) / 2)
  d <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      d <- d + 1
      #' inner product and norms
      inner <- 0
      normx1 <- 0
      normy1 <- 0
      for (l in 1:m) {
        if (mode == "normal") {
          inner <- inner + score[i, l] * score[j, l] / Eig$values[l]
        } else if (mode == "hybrid") {
          inner <- inner + score[i, l] * score[j, l]
        }
        normx1 <- normx1 + (score[i, l] * score[i, l] / Eig$values[l])
        normy1 <- normy1 + (score[j, l] * score[j, l] / Eig$values[l])
      }
      C[d] <- inner / (sqrt(normx1) * sqrt(normy1))
    }
  }

  #' wl-25-10-2020, Sun: convert to symmetric matrix
  if (T) {
    #' wl's implementation
    #' mat <- matrix(1, n, n)
    #' mat[lower.tri(mat, diag = F)] <- C
    #' mat[which(lower.tri(t(mat)), arr.ind = T)[, c(2,1)]] <- C

    #' ji's implemetation
    mat <- matrix(0, n, n)
    mat[lower.tri(mat, diag = F)] <- C
    mat <- mat + t(mat)
    diag(mat) <- 1

    dimnames(mat) <- list(rownames(x), rownames(x))
    return(mat)
  } else {
    return(C)
  }

}
#' 
#' =======================================================================
#' From R package "lsa"
#' x <- iris[, -5]
#' cosine(as.matrix(x))
#' cosine(as.matrix(t(x)))
#'
cosine <- function(x, y = NULL) {
  if (is.matrix(x) && is.null(y)) {
    co <- array(0, c(ncol(x), ncol(x)))
    f <- colnames(x)
    dimnames(co) <- list(f, f)
    for (i in 2:ncol(x)) {
      for (j in 1:(i - 1)) {
        co[i, j] <- cosine(x[, i], x[, j])
      }
    }
    co <- co + t(co)
    diag(co) <- 1
    return(as.matrix(co))
  } else if (is.vector(x) && is.vector(y)) {
    return(crossprod(x, y) / sqrt(crossprod(x) * crossprod(y)))
  } else if (is.vector(x) && is.matrix(y)) {
    co <- vector(mode = "numeric", length = ncol(y))
    names(co) <- colnames(y)
    for (i in 1:ncol(y)) {
      co[i] <- cosine(x, y[, i])
    }
    return(co)
  } else {
    stop("argument mismatch. Either one matrix or two vectors needed as input.")
  }
}
