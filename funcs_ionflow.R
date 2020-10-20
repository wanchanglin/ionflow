
#' =======================================================================
#' wl-12-10-2020, Mon: p_symd is never used. Should remove.
#' wl-20-10-2020, Tue: fix several bugs
PreProcessing <- function(data = NULL, var_id = 1, batch_id = 3, data_id = 5,
                          method_norm = c("median", "median+std", "none"),
                          control_lines = NULL,
                          control_use = c("control", "all", "control out"),
                          method_outliers = c("mad", "IQR", "log.FC.dist", "none"),
                          n_thrs = 4,
                          stand_method = c("std", "mad", "custom"),
                          stdev = NULL, symb_thr = 4) {


  ## -------------------> Import data
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

  ## -------------------> Median batch correction
  data_long$control <- rep(1, length(data_long$Concentration))
  data_long$log <- log(data_long$Concentration)

  #' ji: control use in batch correction
  #' wl-20-10-2020, Tue: select control variables: all, some or some-not.
  #' dirt-trick: use control_lines to decide using control or not.
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
        y
      })
    })
  }
  if (method_norm == "none") {
    data_long$log_corr <- data_long$log
  }

  #' get correction stats
  res <- plyr::ddply(data_long, "Ion", function(x) {
    c(round(summary(x$log_corr), 3), round(var(x$log_corr), 3))
  })
  names(res)[ncol(res)] <- "Variance"
  df_bat <- res

  ## -------------------> Outlier detection
  
  #' ji: outlier detection methods
  if (method_outliers == "IQR") {
    data_long <- plyr::ddply(data_long, "Ion", function(x) {
      res <- plyr::ddply(x, "Batch_ID", function(y) {
        lowerq <- quantile(y$log_corr, na.rm = T)[2]
        upperq <- quantile(y$log_corr, na.rm = T)[4]
        iqr <- upperq - lowerq
        extreme.t.upper <- (iqr * n_thrs) + upperq
        extreme.t.lower <- lowerq - (iqr * n_thrs)
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
        extreme.t.upper <- (med_dev * n_thrs)
        extreme.t.lower <- -(med_dev * n_thrs)
        y$Outlier <- ifelse((y$log_corr > extreme.t.upper) |
          (y$log_corr < extreme.t.lower), 1, 0)
        return(y)
      })
    })
  }
  if (method_outliers == "log.FC.dist") {
    data_long <- plyr::ddply(data_long, "Ion", function(x) {
      res <- plyr::ddply(x, "Batch_ID", function(y) {
        extreme.t.upper <- n_thrs
        extreme.t.lower <- -n_thrs
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

  ## -------------------> Standardisation

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

  #' ji: aggregate measurements at gene level (median)
  data_wide_gene_log_norm <- reshape2::dcast(data_long, Line ~ Ion,
    fun.aggregate = median,
    value.var = "log_corr"
  )

  #' ji: normalise by stds
  data_wide_gene_z_score <- data_wide_gene_log_norm
  data_wide_gene_z_score[, 2:ncol(data_wide_gene_z_score)] <-
    data_wide_gene_z_score[, 2:ncol(data_wide_gene_z_score)] / sds

  ## -------------------> Symbolisation
  symb_profiles <- data_wide_gene_z_score[, 2:ncol(data_wide_gene_z_score)]
  symb_profiles[(symb_profiles > -symb_thr) & (symb_profiles < symb_thr)] <- 0
  symb_profiles[symb_profiles >= symb_thr] <- 1
  symb_profiles[symb_profiles <= -symb_thr] <- -1

  #' wl-20-10-2020, Tue: fix a bug
  data_wide_gene_symb <- cbind(Line = data_wide_gene_z_score$Line,
                               symb_profiles)

  #' wl-20-10-2020, Tue: fix a bug
  p1 <-
    ggplot(
      data = data_long,
      aes(x = factor(Batch_ID), y = log_corr, col = factor(Batch_ID))
    ) +
    geom_point(shape = 1) +
    facet_wrap(~Ion) +
    xlab("Batch.ID") +
    ylab("log(MedianFC)") +
    theme(legend.position = "none") +
    theme(axis.text.x = element_blank())

  #' wl-20-10-2020, Tue: fix a bug
  tmp <- data_wide_gene_z_score 
  dat <- reshape2::melt(tmp, id = "Line")
  p2 <-
    ggplot(data = dat, aes(x = value)) +
    geom_histogram(binwidth = .1) +
    facet_wrap(~variable) +
    xlab("Concentration (z-score)") +
    ylab("Frequency") +
    geom_vline(xintercept = c(-3, 3), col = "red")

  ## -------------------> Output
  res                    <- list()
  res$stats.raw_data     <- df_raw                    # raw data
  res$stats.outliers     <- df_outlier                # outliers
  res$stats.batch_data   <- df_bat                    # batch corrected data
  res$data.long          <- data_long                 # with Batch_ID
  res$data.gene.logFC    <- data_wide_gene_log_norm
  res$data.gene.zscores  <- data_wide_gene_z_score
  res$data.gene.symb     <- data_wide_gene_symb
  res$plot.dot           <- p1
  res$plot.hist          <- p2
  return(res)
}

#' =======================================================================
#'
ExploratoryAnalysis <- function(data = NULL) {

  ## -------------------> Correlation
  col3 <- colorRampPalette(c("steelblue4", "white", "firebrick"))

  corrplot.mixed(cor(data[, -1], use = "complete.obs"),
    number.cex = .7,
    lower.col = "black", upper.col = col3(100)
  )
  p_corr <- recordPlot()

  ## -------------------> PCA
  #' wl-14-07-2020, Tue: Original (trust) pca computation if there is no NAs.
  dat <- t(data[, -1])
  pca <- prcomp(dat, center = T, scale. = F)

  #' variance explained
  vars <- pca$sdev^2
  vars <- vars / sum(vars) #' Proportion of Variance
  names(vars) <- colnames(pca$rotation)
  vars <- round(vars * 100, 2)
  dfn <- paste(names(vars), " (", vars[names(vars)], "%)", sep = "")

  #' PCA scores
  pca_scores <- data.frame(pca$x)
  #' names(pca_scores) <- dfn

  #' PCA loading
  PCA_loadings <- data.frame(pca$rotation)
  rownames(PCA_loadings) <- data$Line
  PCA_loadings <- PCA_loadings[, 1:2]

  #' PCA plot using ggplot2
  pca_p <-
    ggplot(data = pca_scores[, 1:2], aes(x = PC1, y = PC2)) +
    geom_point(color = "steelblue", size = 3, alpha = 0.4) +
    geom_text_repel(aes(label = row.names(pca_scores)), size = 4) +
    theme_bw() +
    theme(legend.position = "none") +
    xlab(dfn[1]) +
    ylab(dfn[2]) +
    labs(title = "PCA")

  ## pca_p

  ## -------------------> HEATMAP
  ## wl-14-08-2020, Fri:  use ggplots
  heatmap.2(as.matrix(data[, -1]),
    scale = "row", col = bluered(100),
    trace = "none", density.info = "none",
    hclustfun = function(x) hclust(x, method = "ward.D")
  )
  ## library(pheatmap)
  ## pheatmap(data[, -1], show_rownames = F, cluster_cols = T, cluster_rows = T,
  ##          legend = T, fontsize = 15, clustering_method = "ward.D",
  ##          scale = "row")
  pheat <- recordPlot()

  ## -------------------> PAIRWISE CORRELATION MAP
  col <- colorRampPalette(c("skyblue4", "white", "plum4"))(20)
  corr <- cor(na.omit(data[, -1]))
  heatmap(
    x = corr, col = col, symm = TRUE, cexRow = 1.4, cexCol = 1.4,
    main = ""
  )
  pcm <- recordPlot()

  #' -------------------> Regularized partial correlation network MAP
  #' wl-13-08-2020, Thu: there is no 'qgraph' in conda forge and bio conda.
  #' Have to plot the correlation network instead.
  if (T) {
    ## wl-14-08-2020, Fri: debug code only
    ## library(glasso)
    ## corr_reg <- glasso(corr, rho = 0.01)
    ## net <- network::network(corr_reg$w, directed = FALSE)
    net <- network::network(corr, directed = FALSE)

    #' set edge attributes
    net %e% "weight" <- corr
    net %e% "weight_abs" <- abs(corr) * 6
    net %e% "color" <- ifelse(net %e% "weight" > 0, "lightgreen", "coral")
    ## set.edge.value(net, "weight", corr)
    ## list.network.attributes(net)
    ## list.edge.attributes(net)
    ## list.vertex.attributes(net)
    net_p <-
      ggnet2(net,
        label = TRUE, mode = "spring",
        node.size = 10, edge.size = "weight_abs", edge.color = "color"
      )
    ## net_p
  } else {
    ## wl-06-07-2020, Mon: 'cor_auto' is from package qgraph(lavaan)
    ## wl-28-07-2020, Tue: cad and corr are the same
    cad <- cor_auto(data[, -1])
    suppressWarnings(qgraph(cad,
      graph = "glasso", layout = "spring",
      tuning = 0.25, sampleSize = nrow(data[, -1])
    ))
    Graph_lasso <- recordPlot()
  }

  ## -------------------> Output
  res <- list()
  res$plot.Pearson_correlation <- p_corr
  res$plot.PCA_Individual <- pca_p
  res$data.PCA_loadings <- PCA_loadings
  res$plot.heatmap <- pheat
  res$plot.pairwise_correlation_map <- pcm
  ## res$plot.regularized_partial_correlation_network <- Graph_lasso
  res$plot.correlation_network <- net_p
  return(res)
}

#' =======================================================================
#'
GeneClustering <- function(data = NULL, data_symb = NULL,
                           min_clust_size = 10, thres_anno = 5) {

  ## -------------------> Define clusters
  # group together strings of symbols at zero Hamming distance
  res.dist <- dist(data_symb[, -1], method = "manhattan")
  res.hc <- hclust(d = res.dist, method = "single")
  #' ji-25-09-2020, Fri: !
  clus <- cutree(res.hc, h = 0)

  data_symb$cluster <- clus

  ## -------------------> Subset cluster with more than 10 genes
  df <- as.data.frame(table(clus), stringsAsFactors = F)
  names(df) <- c("cluster", "nGenes")
  df_sub <- df[df$nGenes > min_clust_size, ]
  rownames(df_sub) <- c()

  #' wl-24-07-2020, Fri: cluster index satisfing threshold of cluster number
  idx <- clus %in% df_sub$cluster
  #' sum(idx)

  mat <- data[idx, ]
  mat$cluster <- clus[idx]

  mat_long <-
    reshape2::melt(mat,
      id = c("Line", "cluster"), variable.name = "Ion",
      value.name = "log_corr_norm"
    )

  res <- sapply(mat_long$cluster, function(x) {
    tmp <- df_sub[df_sub$cluster == x, ]
    tmp <- paste0("Cluster ", tmp[1], " (", tmp[2], " genes)")
  })
  mat_long$cluster <- res #' update cluster with gene numbers

  clus_p <-
    ggplot(
      data = mat_long,
      aes(x = Ion, y = log_corr_norm)
    ) +
    facet_wrap(~cluster) +
    geom_line(aes(group = Line)) +
    stat_summary(fun.data = "mean_se", color = "red") +
    labs(x = "", y = "") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text = element_text(size = 10)
    )

  #' =======================================================================
  #'
  ## -------------------> KEGG AND GO SLIM ANNOTATION
  mat <- data_symb[idx, ]
  data_GOslim$Ontology <- as.character(data_GOslim$Ontology)

  kego <- plyr::dlply(mat, "cluster", function(x) {
    ## x <- subset(mat, cluster == "15")
    inputGeneSet <- as.character(x$Line)
    N <- length(inputGeneSet)

    res <- data_GOslim %>%
      dplyr::mutate(Ontology = setNames(
        c(
          "Biological process",
          "Cellular component",
          "Molecular function"
        ),
        c("P", "C", "F")
      )[Ontology]) %>%
      dplyr::filter(ORFs %in% inputGeneSet) %>%
      dplyr::group_by(GOslim, Ontology) %>%
      dplyr::filter(GOslim != "other") %>%
      dplyr::rename(Term = GOslim) %>%
      dplyr::summarise(Count = n()) %>%
      dplyr::mutate(Percent = Count / N * 100) %>%
      dplyr::bind_rows(data_ORF2KEGG %>%
        dplyr::filter(ORF %in% inputGeneSet) %>%
        dplyr::group_by(KEGGID, Pathway) %>%
        dplyr::summarise(Count = n()) %>%
        dplyr::mutate(Ontology = "KEGG") %>%
        dplyr::rename(Term = Pathway) %>%
        dplyr::ungroup() %>%
        dplyr::filter(KEGGID != "01100") %>%
        dplyr::select(-KEGGID) %>%
        dplyr::mutate(Percent = Count / N * 100)) %>%
      dplyr::filter(!Term %in% c(
        "molecular_function", "biological_process",
        "cellular_component"
      ))
  })
  names(kego) <- paste0("Cluster ", df_sub[[1]], " (", df_sub[[2]], " genes)")

  #' wl-25-07-2020, Sat: filter annotation results. Should set threshold for
  #' Percent?
  kego <- lapply(kego, function(x) {
    ## x[(x$Percent > 5) &
    ##   (x$Ontology %in% c("Biological process", "Cellular component",
    ##                      "Molecular function")), ]
    x[x$Percent > thres_anno, ]
  })

  #' wl-04-08-2020, Tue: bind together
  kego <- dplyr::bind_rows(kego, .id = "Cluster")

  ## -------------------> GO TERMS ENRICHMENT

  #' wl-04-08-2020, Tue: re-write
  universeGenes <- as.character(data_symb$Line)
  mat <- data_symb[idx, ]

  goen <- plyr::dlply(mat, "cluster", function(x) {
    ## x <- subset(mat, cluster == "10")
    inputGeneSet <- as.character(x$Line)
    ont <- c("BP", "MF", "CC") ## wl-04-08-2020, Tue: why three?
    res <- lapply(ont, function(y) {
      params <- new("GOHyperGParams",
        geneIds = inputGeneSet,
        universeGeneIds = universeGenes,
        annotation = "org.Sc.sgd.db",
        categoryName = "GO",
        ontology = y,
        pvalueCutoff = 0.05,
        conditional = T,
        testDirection = "over"
      )
      hyperGTest(params)
    })
    names(res) <- ont

    #' extract some results and move out filtering
    res_1 <- lapply(ont, function(y) {
      hgOver <- res[[y]]
      tmp <- cbind(setNames(
        tibble(
          ID = names(pvalues(hgOver)),
          Term = Term(ID),
          pvalues = pvalues(hgOver),
          oddsRatios = oddsRatios(hgOver),
          expectedCounts = expectedCounts(hgOver),
          geneCounts = geneCounts(hgOver),
          universeCounts = universeCounts(hgOver)
        ),
        c(
          "GO_ID", "Description", "Pvalue", "OddsRatio",
          "ExpCount", "Count", "CountUniverse"
        )
      ),
      Ontology = y
      ) ## %>% dplyr::filter(Pvalue <= 0.05 & Count > 1)
    })
    res_2 <- do.call("rbind", res_1)
  })

  names(goen) <- paste0("Cluster ", df_sub[[1]], " (", df_sub[[2]], " genes)")

  #' binding and filtering
  goen <- lapply(goen, "[", -c(4, 5, 8)) %>%
    dplyr::bind_rows(.id = "Cluster") %>%
    dplyr::filter(Pvalue <= 0.05 & Count > 1)

  ## -------------------> Output
  res <- list()
  res$stats.clusters <- df_sub # selected clusters
  res$plot.profiles <- clus_p # plot cluster profiles
  res$stats.Kegg_Goslim_annotation <- kego # KEGG AND GO SLIM ANNOTATION
  res$stats.Goterms_enrichment <- goen # GO TERMS ENRICHMENT
  return(res)
}

#' =======================================================================
#' wl-19-10-2020, Mon: fix some bugs.
GeneNetwork <- function(data = NULL, data_symb = NULL,
                        min_clust_size = 10, thres_corr = 0.6,
                        method_corr = c(
                          "pearson", "spearman", "kendall",
                          "cosine", "mahal_cosine", "hybrid_mahal_cosine"
                        )) {

  #' wl-19-10-2020, Mon: add this line to make default
  method_corr <- match.arg(method_corr)

  ## Cluster of gene with same profile
  #' ji-22-09-2020, Tue: clustring already performed ?
  res.dist <- dist(data_symb[, -1], method = "manhattan")
  res.hc <- hclust(d = res.dist, method = "single")
  #' ji-25-09-2020, Fri: equivalent to zero Hamming distance !
  symb.cluster <- cutree(res.hc, h = 0) # distance 0

  data_symb$cluster <- symb.cluster

  df <- as.data.frame(table(symb.cluster), stringsAsFactors = F)
  names(df) <- c("cluster", "nGenes")
  ## filter clusters with threshold
  df_sub <- df[df$nGenes > min_clust_size, ]

  #' wl-26-07-2020, Sun: remove the largest clusters?
  ## Cluster 2 (largest cluster) contains genes with no phenotype hence not
  ## considered (input to 0)
  if (T) df_sub <- df_sub[-which.max(df_sub[, 2]), ]

  rownames(df_sub) <- c()

  #' wl-24-07-2020, Fri: cluster index satisfing threshold of cluster number
  index <- symb.cluster %in% df_sub$cluster

  ## cluster labels with info of accumulation/decumulation of Ions
  ## (high/lower abundance)
  ## Assign label
  df.symb <- data_symb[index, ]
  lab <- plyr::ddply(df.symb, "cluster", function(x) {
    mat <- x[, !names(df.symb) %in% c("Line", "cluster")]
    res <- NULL
    if (length(names(which(colSums(mat == 1) > 0))) > 0) {
      res <- paste0(names(which(colSums(mat == 1) > 0)), "(+)")
    }
    if (length(names(which(colSums(mat == -1) > 0))) > 0) {
      res <- c(res, paste0(names(which(colSums(mat == -1) > 0)), "(-)"))
    }
    res <- paste(res, collapse = ", ")
  })
  names(lab)[2] <- "label"

  #' update df_sub with labels
  df_sub <- cbind(df_sub, label = lab[, 2])

  #' get cluster+gene+label names
  label <- sapply(df.symb$cluster, function(x) {
    tmp <- df_sub[df_sub$cluster == x, ]
    res <- paste0("Cluster ", tmp[1], " (", tmp[2], " genes)")
    res <- paste(res, tmp[3], sep = ": ")
  })
  df.symb$Label <- label

  ## Compute empirical correlation matrix
  if (method_corr == "pearson" ||
      method_corr == "spearman" ||
      method_corr == "kendall") {
    corrGenes <- cor(t(as.matrix(data[, -1])), method = method_corr,
                     use = "pairwise.complete.obs")
  } else if (method_corr == "cosine") {
    corrGenes <- cosine(t(as.matrix(data[, -1])))
  } else if (method_corr == "mahal_cosine") {
    corrGenes <- cosM(t(as.matrix(data[, -1])), node = "normal")
  } else if (method_corr == "hybrid_mahal_cosine") {
    corrGenes <- cosM(t(as.matrix(data[, -1])), node = "hybrid")
  }

  ## Subset correlation matrix based on the cluster filtering
  A <- corrGenes[index, index]
  ## Diagonal value (1's) put to 0 to avoid showing edges from/to the same gene
  diag(A) <- 0

  ## Subset correlation matrix based on threshold=0.6
  ## wl-27-07-2020, Mon: need another threshold?
  A <- (A > thres_corr)
  A <- ifelse(A == TRUE, 1, 0)

  ## Generate network
  net <- network::network(A, directed = FALSE)

  #' wl-28-07-2020, Tue: add an vertex attribute and use 'Set2' in
  #'  RColorBrewer but the max. number of colors is 8 in 'Set2'
  #' wl-29-07-2020, Wed: Some layouts: circle, fruchtermanreingold,
  #'  kamadakawai, spring
  net %v% "Label" <- df.symb$Label
  tmp <- unique(df.symb$Label)
  ## fix a bug (may from ggnet2)
  if (length(tmp) != 1) {
    cpy <- rainbow(length(tmp))
    names(cpy) <- tmp
  } else {
    cpy <- "Set2"
  }
  net_p <- GGally::ggnet2(net,
    mode = "fruchtermanreingold",
    color = "Label",
    palette = cpy,
    edge.alpha = 0.5, size = 2, color.legend = "Label",
    legend.size = 10, legend.position = "right"
  )
  #' net_p

  ## Impact and betweenness
  btw <- sna::betweenness(A) # or use 'net' instead of 'A'
  impact <- apply(data[index, -1], 1, norm, type = "2") # L2 norm

  df.res <- data.frame(
    Line = data$Line[index],
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
  rownames(df.res) <- data$Line[index]

  q1 <- row.names(subset(df.res, (impact < quantile(impact, .75)) & (log.betweenness < quantile(log.betweenness, .75))))
  q2 <- row.names(subset(df.res, (impact < quantile(impact, .75)) & (log.betweenness > quantile(log.betweenness, .75))))
  q3 <- row.names(subset(df.res, (impact > quantile(impact, .75)) & (log.betweenness < quantile(log.betweenness, .75))))
  q4 <- row.names(subset(df.res, (impact > quantile(impact, .75)) & (log.betweenness > quantile(log.betweenness, .75))))

  #' idx <- unique(c(sample(q1,6),sample(q2,6),sample(q3,6),sample(q4,6)))
  #' wl-27-07-2020, Mon: random choose at least 24 genes to show in plot
  #' wl-17-07-2020, Fri: any or all of q1~q4 may be character(0)
  #' wl-14-07-2020, Tue: potential bug in sample replacement
  N <- 6
  lst <- list(q1, q2, q3, q4) #' sapply(lst, length)
  idx <- lapply(lst, function(x) {
    if (length(x) > N) sample(x, N) else x
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

  rownames(df.res) <- c()
  df.res2 <- df.res[, -c(4, 5)]
  names(df.res2) <- c("Line", "Impact", "Betweenness", "Position")

  #' wl-19-10-2020, Mon: bug here. Need to go back to PreProcesing. The
  #' following line is temporary.
  #' names(df.symb)[1] <- "Line" 
  gene.cluster <- df.symb[, c("Line", "Label")]

  names(gene.cluster) <- c("Line", "Cluster")
  df.res3 <- merge(df.res2, gene.cluster, by = "Line", all.x = TRUE)

  #' wl-28-07-2020, Tue: better to return df.tab instead of df.tab2
  df.tab <- data.frame(table(df.res3$Cluster, df.res3$Position))
  names(df.tab) <- c("Cluster", "Position", "nGenes")
  df.tab <- dplyr::arrange(df.tab, desc(nGenes))
  df.tab2 <- df.tab %>% dplyr::group_by(Cluster) %>% top_n(1, nGenes)

  ## -------------------> Output
  res <- list()
  res$plot.pnet <- net_p # plot gene network
  res$plot.impact_betweenness <- im_be_p # plot impact betweenees
  res$stats.impact_betweenness <- df.res3 # impact betweenees data
  res$stats.impact_betweenness_by_cluster <- df.tab2 # plot position by cluster
  ## wl-28-07-2020, Tue: return this one as well
  res$stats.impact_betweenness_tab <- df.tab # contingency table
  return(res)
}

#' =======================================================================
#' Mahalanobis Cosine
#' Function to compute the mahalanobis cosine between pairs of objects in an
#' n-by-m data matrix or data frame.
cosM <- function(x, mode = c("normal", "hybrid")) {
  #' library("pracma")
  #' compute covariance
  Cov <- cov(x, use = "pairwise.complete.obs")
  n <- dim(x)[1]
  m <- dim(x)[2]

  #' compute eigenvalues
  Eig <- eigen(Cov)
  score <- mrdivide(as.matrix(x), t(Eig$vectors))

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
  return(C)
}
