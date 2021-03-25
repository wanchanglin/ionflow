
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

  #' ji: aggregate measurements at gene level (median)
  data_wide_gene_log_norm <- reshape2::dcast(data_long, Line ~ Ion,
    fun.aggregate = median,
    value.var = "log_corr"
  )

  #' ji: normalise by stds
  data_wide_gene_z_score <- data_wide_gene_log_norm
  data_wide_gene_z_score[, 2:ncol(data_wide_gene_z_score)] <-
    data_wide_gene_z_score[, 2:ncol(data_wide_gene_z_score)] / sds

  #' -------------------> Symbolisation
  symb_profiles <- data_wide_gene_z_score[, 2:ncol(data_wide_gene_z_score)]
  symb_profiles[(symb_profiles > -thres_symb) & (symb_profiles < thres_symb)] <- 0
  symb_profiles[symb_profiles >= thres_symb] <- 1
  symb_profiles[symb_profiles <= -thres_symb] <- -1

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

  #' -------------------> Output
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

  #' -------------------> PCA
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
  pca_loadings <- data.frame(pca$rotation)
  rownames(pca_loadings) <- data$Line
  pca_loadings <- pca_loadings[, 1:2]

  #' PCA plot using ggplot2
  p_pca <-
    ggplot(data = pca_scores[, 1:2], aes(x = PC1, y = PC2)) +
    geom_point(color = "steelblue", size = 3, alpha = 0.4) +
    geom_text_repel(aes(label = row.names(pca_scores)), size = 4) +
    theme_bw() +
    theme(legend.position = "none") +
    xlab(dfn[1]) +
    ylab(dfn[2]) +
    labs(title = "PCA")

  #' -------------------> Correlation
  col3 <- colorRampPalette(c("steelblue4", "white", "firebrick"))
  corrplot.mixed(cor(data[, -1], use = "complete.obs"),
    number.cex = .7,
    lower.col = "black", upper.col = col3(100)
  )
  p_corr <- recordPlot()

  #' -------------------> Correlation heatmap
  col <- colorRampPalette(c("skyblue4", "white", "plum4"))(20)
  corr <- cor(na.omit(data[, -1]))
  heatmap(
    x = corr, col = col, symm = TRUE, cexRow = 1.4, cexCol = 1.4,
    main = ""
  )
  p_corr_heat <- recordPlot()

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
    net <- network::network(corr, directed = FALSE)

    #' set edge attributes
    net %e% "weight" <- corr
    net %e% "weight_abs" <- abs(corr) * 6
    net %e% "color" <- ifelse(net %e% "weight" > 0, "lightgreen", "coral")
    p_net <-
      ggnet2(net,
        label = TRUE, mode = "spring",
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
  res$plot.pca  <- p_pca
  res$data.pca.load <- pca_loadings
  res$plot.corr <- p_corr
  res$plot.corr.heat <- p_corr_heat
  res$plot.heat <- p_heat
  res$plot.net <- p_net
  return(res)
}

#' =======================================================================
#'
GeneNetwork <- function(data = NULL, data_symb = NULL,
                        min_clust_size = 10, thres_corr = 0.6,
                        method_corr = c(
                          "pearson", "spearman", "kendall",
                          "cosine", "mahal_cosine", "hybrid_mahal_cosine"
                        )) {

  method_corr <- match.arg(method_corr)

  #' Define clusters
  clust <- gene_clus(data_symb[, -1], min_clust_size = min_clust_size)
  data_symb$cluster <- clust$clus
  index <- clust$idx
  df_sub <- clust$tab_sub

  #' cluster labels with info of accumulation/decumulation of Ions
  #' (high/lower abundance) Assign label
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

  #' Compute empirical correlation matrix
  if (method_corr == "pearson" ||
      method_corr == "spearman" ||
      method_corr == "kendall") {
    corrGenes <- cor(t(as.matrix(data[, -1])), method = method_corr,
                     use = "pairwise.complete.obs")
  } else if (method_corr == "cosine") {
    corrGenes <- cosine(t(as.matrix(data[, -1])))
  } else if (method_corr == "mahal_cosine") {
    corrGenes <- cosM(data[, -1], mode = "normal")
  } else if (method_corr == "hybrid_mahal_cosine") {
    corrGenes <- cosM(data[, -1], mode = "hybrid")
  }
  
  #' wl-15-12-2020, Tue: assign names
  gene <- as.character(data[, 1])
  dimnames(corrGenes) <- list(gene, gene)

  A <- corrGenes[index, index]
  diag(A) <- 0

  A <- (A > thres_corr)
  A <- ifelse(A == TRUE, 1, 0)

  #' Generate network
  net <- network::network(A, directed = FALSE)

  net %v% "Label" <- df.symb$Label
  tmp <- unique(df.symb$Label)
  if (length(tmp) != 1) {
    cpy <- rainbow(length(tmp))
    names(cpy) <- tmp
  } else {
    cpy <- "Set2"
  }
  net_p1 <- GGally::ggnet2(net,
    mode = "fruchtermanreingold",
    color = "Label",
    palette = cpy,
    edge.alpha = 0.5, size = 2, color.legend = "Symbolic Pheno",
    legend.size = 10, legend.position = "right"
  )
  #' net_p1

  #' wl-09-11-2020, Mon: node color controoled by community detection
  tmp <- igraph::graph_from_adjacency_matrix(A, mode="undirected")
  com <- igraph::cluster_louvain(tmp, weights = NULL)
  mem <- as.factor(igraph::membership(com))

  #' wl-16-12-2020, Wed: add more stuff
  tab <- table(mem)
  mem <- sapply(mem, function(x) {
    res <- paste0("Cluster ", x, " (", tab[x], " genes)")
  })

  tmp <- unique(mem)
  if (length(tmp) != 1) {
    cpy <- rainbow(length(tmp))
    names(cpy) <- tmp
  } else {
    cpy <- "Set2"
  }
  net_p2 <- GGally::ggnet2(net,
    mode = "fruchtermanreingold",
    color = mem,
    palette = cpy,
    edge.alpha = 0.5, size = 2, color.legend = "Network Community",
    legend.size = 10, legend.position = "right"
  )
  #' net_p2

  #' Impact and betweenness
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

  gene.cluster <- df.symb[, c("Line", "Label")]

  names(gene.cluster) <- c("Line", "Cluster")
  df.res3 <- merge(df.res2, gene.cluster, by = "Line", all.x = TRUE)

  df.tab <- data.frame(table(df.res3$Cluster, df.res3$Position))
  names(df.tab) <- c("Cluster", "Position", "nGenes")
  df.tab <- dplyr::arrange(df.tab, desc(nGenes))
  df.tab2 <- df.tab %>% dplyr::group_by(Cluster) %>% top_n(1, nGenes)

  #' wl-15-12-2020, Tue: get network vertex attributes
  node_names <- net %v% "vertex.names"
  symb_pheno <- net %v% "Label"
  comm_centre <- mem
  net_node <- data.frame(Line = node_names, symb_pheno = symb_pheno, 
                         comm_centre = comm_centre)

  res <- list()
  res$plot.pnet1 <- net_p1 # plot gene network with symbolic pheno
  res$plot.pnet2 <- net_p2 # plot gene network with community detection
  res$plot.impact_betweenness <- im_be_p # plot impact betweenees
  res$stats.impact_betweenness <- df.res3 # impact betweenees data
  res$stats.impact_betweenness_by_cluster <- df.tab2 # plot position by cluster
  res$stats.impact_betweenness_tab <- df.tab # contingency table
  res$net_node <- net_node
  return(res)
}

#' =======================================================================
#' wl-04-10-2020, Sun: Hierarchical clustering
#'
gene_clus <- function(x, min_clust_size = 10) {
  dis <- stats::dist(x, method = "manhattan")
  hc <- hclust(d = dis, method = "single")
  clus <- cutree(hc, h = 0)

  tab <- as.data.frame(table(clus), stringsAsFactors = F)
  names(tab) <- c("cluster", "nGenes")
  tab_sub <- tab[tab$nGenes > min_clust_size, ]
  tab_sub <- tab_sub[order(tab_sub$nGenes, decreasing = T), ]
  rownames(tab_sub) <- NULL

  idx <- clus %in% tab_sub$cluster

  return(list(clus = clus, idx = idx, tab = tab, tab_sub = tab_sub))
}

#' =======================================================================
#' ji: Mahalanobis Cosine
#' Function to compute the mahalanobis cosine between pairs of objects in an
#' n-by-m data matrix or data frame.
#' x <- iris[, -5]
#' res <- cosM(x, mode = "normal")
#'
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

#' =======================================================================
#' Get symbolisation profiles based on z-scores in wide format
#'
symbol_data <- function(x, thres_symb = 2) {
  symb <- x[, 2:ncol(x)]
  symb[(symb > -thres_symb) & (symb < thres_symb)] <- 0
  symb[symb >= thres_symb] <- 1
  symb[symb <= -thres_symb] <- -1
  symb <- cbind(Line = x[, 1], symb)
  return(symb)
}

#' =======================================================================
#' wl-06-11-2020, Fri: Get ENTREZID  from SYMBOL
#'
get_entrez_id <- function(symbol, annot_pkg = "org.Hs.eg.db") {
  res <- AnnotationDbi::select(get(annot_pkg), keys = symbol,
                               columns = "ENTREZID", keytype = "SYMBOL")
  res <- res[,2,drop = T]
  res <- res[!is.na(res)]
  res <- res[!duplicated(res)]

  return(res)
}

#' =======================================================================
#' wl-03-10-2020, Sat: KEGG enrichment analysis for symbolization data
#' wl-06-11-2020, Fri: The first column of data must be ORF for
#'  org.Sc.sgd.db or SYMBOL for any other annotation packages.
#' wl-15-12-2020, Tue: KEGG enrichment based on network comunity centre
#' wl-15-12-2020, Tue: KEGG enrichment for group data
#'  The group information may be the symbolic clustering or network
#'  community detection
#'
kegg_enrich <- function(mat, pval = 0.05, annot_pkg = "org.Sc.sgd.db") {

  #' get the gene ids
  gene_uni <- as.character(mat[, 1])
  grp <- names(mat[2])
  gene_ids <- plyr::dlply(mat, grp, function(x) as.character(x[, 1]))
  ind <- sapply(gene_ids, function(x) ifelse(length(x) > 1, TRUE, FALSE))
  gene_ids <- gene_ids[ind]

  #' get entrez id for GOStats
  if (annot_pkg != "org.Sc.sgd.db") {
    gene_uni <- get_entrez_id(gene_uni, annot_pkg)
    gene_ids <- lapply(gene_ids, function(x){
      get_entrez_id(x, annot_pkg)
    })
  }

  #' geneIds can be ORF or ENTREZID
  enrich <- lapply(gene_ids, function(x) { #' x = gene_ids[[1]]
    params <- new("KEGGHyperGParams",
                  geneIds = x,
                  universeGeneIds = gene_uni,
                  annotation = annot_pkg,
                  categoryName = "KEGG",
                  pvalueCutoff = 1,
                  testDirection = "over")

    over <- hyperGTest(params)
  })

  #' There is no explicit methods for getting manual summary table.
  summ <- lapply(enrich, function(x) {
    tmp <- summary(x)
    if (nrow(tmp) == 0) tmp <- NULL #' wl-04-10-2020, Sun: it happens very often.
    return(tmp)
  })
  summ <- summ[!sapply(summ,is.null)]

  #' binding and filtering
  summ <- lapply(summ, "[", -c(3, 4)) %>%
    dplyr::bind_rows(.id = grp) %>%
    dplyr::filter(Pvalue <= pval & Count > 1)

  return(summ)
}

#' =======================================================================
#' wl-03-10-2020, Sat: GO enrichment analysis for symbolization data
#' wl-06-11-2020, Fri: The first column of data must be ORF for
#'  org.Sc.sgd.db or SYMBOL for any other annotation packages.
#' wl-15-12-2020, Tue: GO enrichment based on network comunity centre
#' wl-15-12-2020, Tue: GO enrichment for group data
#'
go_enrich <- function(mat, pval = 0.05, ont = "BP",
                      annot_pkg = "org.Sc.sgd.db") {

  ont <- match.arg(ont, c("BP", "MF", "CC"))

  #' get the gene ids
  gene_uni <- as.character(mat[, 1])
  grp <- names(mat[2])
  gene_ids <- plyr::dlply(mat, grp, function(x) as.character(x[, 1]))
  ind <- sapply(gene_ids, function(x) ifelse(length(x) > 1, TRUE, FALSE))
  gene_ids <- gene_ids[ind]

  #' get entrez id for GOStats
  if (annot_pkg != "org.Sc.sgd.db") {
    gene_uni <- get_entrez_id(gene_uni, annot_pkg)
    gene_ids <- lapply(gene_ids, function(x){
      get_entrez_id(x, annot_pkg)
    })
  }

  #' geneIds can be ORF or ENTREZID
  enrich <- lapply(gene_ids, function(x) { #' x = gene_ids[[1]]
    params <- new("GOHyperGParams",
                  geneIds = x,
                  universeGeneIds = gene_uni,
                  annotation = annot_pkg,
                  categoryName = "GO",
                  ontology = ont,
                  pvalueCutoff = 1,
                  conditional = T,
                  testDirection = "over")

    over <- hyperGTest(params)
  })

  summ <- lapply(enrich, function(x) {
    Pvalue        <- round(pvalues(x), digit = 4)
    ID            <- names(Pvalue)
    Description   <- Term(ID)
    OddsRatio     <- oddsRatios(x)
    ExpCount      <- expectedCounts(x)
    Count         <- geneCounts(x)
    CountUniverse <- universeCounts(x)

    tab <- cbind(ID, Description, Pvalue, OddsRatio, ExpCount, Count,
                 CountUniverse)
    rownames(tab) <- NULL
    tab <- na.omit(tab)   #' wl-03-10-2020, Sat: this is why summary fails.
    tab <- data.frame(tab, Ontology = ont)
  })

  summ <- lapply(summ, "[", -c(4, 5)) %>%
    dplyr::bind_rows(.id = grp) %>%
    dplyr::filter(Pvalue <= pval & Count > 1)

  return(summ)
}

