#' ---
#' title: "Network and enrichment analysis for the ionomics data in a publication"
#' author: "Wanchang Lin"
#' date: "`r Sys.Date()`"
#' ---
#' 
## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      fig.width = 7, fig.height = 7, 
                      fig.align = "center",
                      fig.path = "./img/",
                      dpi = 100, dev = "png", cache = TRUE)
pkgs <- c("mtExtra", "tidyverse", "igraph", "ggraph", "ggrepel", 
          "proxy", "org.Sc.sgd.db", "GO.db", "GOstats", "KEGG.db",
          "knitr", "kableExtra")
invisible(lapply(pkgs, library, character.only = TRUE))
## =======================================================================
## wl-03-10-2020, Sat: GEGG enrichment analysis
kegg_enrich <- function(data,
                        method_simil = "correlation",
                        method_hclus = "ave",
                        thres_simil = 0.6,
                        thres_clus = 10,
                        pval = 0.05,
                        annot_pkg = "org.Sc.sgd.db",
                        is_symb = FALSE) {

  ## define clusters
  if (is_symb) {
    clust <- gene_clus(data[, -1], thres_clus = thres_clus)
  } else {
    clust <- simil_clus(data[, -1],
      method_simil = method_simil,
      method_hclus = method_hclus,
      thres_simil = thres_simil,
      thres_clus = thres_clus
    )
  }

  ## update data for enrichment analysis
  data$cluster <- clust$clus
  mat <- data[clust$idx, ]

  ## get the gene ids
  gene_uni <- as.character(data[, 1, drop = T])
  gene_ids <- mat %>% 
    group_by(cluster) %>%
    group_map(~ as.character(.x[[1]])) %>%
    setNames(unique(sort(mat[, "cluster", drop = T])))

  ## get entrez id for GOStats
  if (annot_pkg != "org.Sc.sgd.db") {
    gene_uni <- get_entrez_id(gene_uni, annot_pkg)
    gene_ids <- lapply(gene_ids, function(x) {
      get_entrez_id(x, annot_pkg)
    })
  }

  ## geneIds can be ORF or ENTREZID
  enrich <- lapply(gene_ids, function(x) {
    params <- new("KEGGHyperGParams",
      geneIds = x,
      universeGeneIds = gene_uni,
      annotation = annot_pkg,
      categoryName = "KEGG",
      pvalueCutoff = 1,
      testDirection = "over"
    )

    over <- hyperGTest(params)
  })

  ## summary
  summ <- lapply(enrich, function(x) {
    tmp <- summary(x)
    if (nrow(tmp) == 0) tmp <- NULL
    return(tmp)
  })
  summ <- summ[!sapply(summ, is.null)]

  ## change names based on updated tab_sub
  tab_sub <- clust$tab_sub
  tmp <- names(summ)
  idx <- tab_sub[, 1] %in% tmp
  tab_sub <- tab_sub[idx, ]
  names(summ) <- paste0("Cluster ", tab_sub[[1]], " (", tab_sub[[2]], " genes)")

  ## binding and filtering
  summ <- lapply(summ, "[", -c(3, 4)) %>%
    dplyr::bind_rows(.id = "Cluster") %>%
    dplyr::filter(Pvalue <= pval & Count > 1)

  return(summ)
}

## =======================================================================
## wl-03-10-2020, Sat: GO enrichment analysis
go_enrich <- function(data,
                      method_simil = "correlation",
                      method_hclus = "ave",
                      thres_simil = 0.6,
                      thres_clus = 10,
                      pval = 0.05,
                      ont = "BP",
                      annot_pkg = "org.Sc.sgd.db",
                      is_symb = FALSE) {
  ont <- match.arg(ont, c("BP", "MF", "CC"))

  ## define clusters
  if (is_symb) {
    clust <- gene_clus(data[, -1], thres_clus = thres_clus)
  } else {
    clust <- simil_clus(data[, -1],
      method_simil = method_simil,
      method_hclus = method_hclus,
      thres_simil = thres_simil,
      thres_clus = thres_clus
    )
  }

  ## update data for enrichment analysis
  data$cluster <- clust$clus
  mat <- data[clust$idx, ]

  ## get the gene ids
  gene_uni <- as.character(data[, 1, drop = T])
  gene_ids <- mat %>% 
    group_by(cluster) %>%
    group_map(~ as.character(.x[[1]])) %>%
    setNames(unique(sort(mat[, "cluster", drop = T])))

  ## get entrez id for GOStats
  if (annot_pkg != "org.Sc.sgd.db") {
    gene_uni <- get_entrez_id(gene_uni, annot_pkg)
    gene_ids <- lapply(gene_ids, function(x) {
      get_entrez_id(x, annot_pkg)
    })
  }

  ## geneIds can be ORF or ENTREZID
  enrich <- lapply(gene_ids, function(x) {
    params <- new("GOHyperGParams",
      geneIds = x,
      universeGeneIds = gene_uni,
      annotation = annot_pkg,
      categoryName = "GO",
      ontology = ont,
      pvalueCutoff = 1,
      conditional = T,
      testDirection = "over"
    )

    over <- hyperGTest(params)
  })

  summ <- lapply(enrich, function(x) {

    Pvalue <- round(pvalues(x), digit = 4)
    ID <- names(Pvalue)
    Description <- Term(ID)
    OddsRatio <- oddsRatios(x)
    ExpCount <- expectedCounts(x)
    Count <- geneCounts(x)
    CountUniverse <- universeCounts(x)

    tab <- cbind(
      ID, Description, Pvalue, OddsRatio, ExpCount, Count,
      CountUniverse
    )
    rownames(tab) <- NULL
    tab <- na.omit(tab)
    tab <- data.frame(tab, Ontology = ont)
  })

  names(summ) <- paste0(
    "Cluster ", clust$tab_sub[[1]], " (",
    clust$tab_sub[[2]], " genes)"
  )

  ## binding and filtering
  summ <- lapply(summ, "[", -c(4, 5)) %>%
    dplyr::bind_rows(.id = "Cluster") %>%
    dplyr::filter(Pvalue <= pval & Count > 1)

  return(summ)
}

## =======================================================================
## wl-06-11-2020, Fri: Get ENTREZID  from SYMBOL
get_entrez_id <- function(symbol, annot_pkg = "org.Hs.eg.db") {
  res <- AnnotationDbi::select(get(annot_pkg),
    keys = symbol,
    columns = "ENTREZID", keytype = "SYMBOL"
  )
  res <- res[, 2, drop = T]
  res <- res[!is.na(res)]
  res <- res[!duplicated(res)]

  return(res)
}

## =======================================================================
## wl-04-10-2020, Sun: Hierarchical clustering on symbolic data
gene_clus <- function(x, thres_clus = 10) {
  dis <- stats::dist(x, method = "manhattan")
  hc <- hclust(d = dis, method = "single")
  clus <- cutree(hc, h = 0) # distance 0

  tab <- as.data.frame(table(clus), stringsAsFactors = F)
  names(tab) <- c("cluster", "nGenes")
  tab_sub <- tab[tab$nGenes > thres_clus, ]
  tab_sub <- tab_sub[order(tab_sub$nGenes, decreasing = T), ]
  rownames(tab_sub) <- NULL

  idx <- clus %in% tab_sub$cluster

  return(list(clus = clus, idx = idx, tab = tab, tab_sub = tab_sub))
}

## =======================================================================
## wl-14-10-2020, Wed: Hierarchical clustering based on similarity measures
simil_clus <- function(x, 
                       method_simil = "correlation", 
                       method_hclus = "ave",
                       thres_simil = 0.6, 
                       thres_clus = 5) {

  ## get similarity measures
  sim <- df_simil(x, method = method_simil)

  ## hierachical clustering based on dissimilarity
  hc <- hclust(d = as.dist(1 - sim), method = method_hclus)

  ## cut-off in threshold
  clus <- cutree(hc, h = 1 - thres_simil)

  tab <- as.data.frame(table(clus), stringsAsFactors = F)
  names(tab) <- c("cluster", "nGenes")
  tab_sub <- tab[tab$nGenes > thres_clus, ]
  tab_sub <- tab_sub[order(tab_sub$nGenes, decreasing = T), ]
  rownames(tab_sub) <- NULL

  idx <- clus %in% tab_sub$cluster

  return(list(sim = sim, clus = clus, idx = idx, tab = tab,
              tab_sub = tab_sub))
}

## ========================================================================
## wl-18-09-2020, Fri: row-wise similarity of data matrix.
df_simil <- function(x, method = "correlation") {
  method <- match.arg(method, 
                      c("correlation", "cosine", "eJaccard", "Mahalanobis"))

  if (method == "Mahalanobis") {
    res <- proxy::dist(x, x, method = "Mahalanobis")
    res <- as.data.frame.matrix(res)

    ## convert distance to similarity
    if (F) {
      res <- 1 / (1 + res)
    } else {
      max_all <- max(res, na.rm = T)
      min_all <- min(res, na.rm = T)
      ## set the MARGIN to 1:2 to operate on each cell
      res <- apply(res, 1:2, function(x) (x - min_all) / (max_all - min_all))
      res <- 1 - res
    }
  } else {
    res <- proxy::simil(x, x, method = method)
    res <- as.data.frame.matrix(res)
  }
  return(res)
}

## =======================================================================
## wl-14-10-2020, Wed: Network analysis using ionomics data
gene_net <- function(data,
                     method_simil = "correlation",
                     method_hclus = "ave",
                     thres_simil = 0.6,
                     thres_clus = 10) {

  data <- as.data.frame(data)
  rownames(data) <- data[, 1]

  ## define clusters
  clust <- simil_clus(data[, -1],
    method_simil = method_simil,
    method_hclus = method_hclus, thres_clus = thres_clus,
    thres_simil = thres_simil
  )

  ## update data
  data$cluster <- clust$clus
  index <- clust$idx
  df_sub <- clust$tab_sub
  data <- data[index, ]

  ## get similarity of data set
  sim <- clust$sim
  sim <- sim[index, index]

  ## create igraph object
  edge <- sym2long(sim, tri = "upper")
  idx <- abs(edge$var) >= thres_simil
  edge <- edge[idx, ]
  edge <- mutate(edge, mod = ifelse(var > 0, "pos", "neg"))

  ## only positive or negative?
  if (F) {
    if (T) {
      idx <- edge$var > 0 ## positive
    } else {
      idx <- edge$var < 0 ## negative
    }
    edge <- edge[idx, ]
  }

  ## get the node table
  num <- as.vector(table(clust$clus))
  num <- num[clust$clus]
  clus <- paste0(clust$clus, "(", num, " genes)")
  node <- data.frame(name = names(clust$clus), clus = clus)

  g <- graph_from_data_frame(d = edge, directed = F, vertices = node)

  ## delete isolated vertices
  g <- delete_vertices(g, which(igraph::degree(g) < 1))

  com <- cluster_louvain(g, weights = NULL)
  mem <- as.factor(membership(com))

  lay <- create_layout(g, layout = "fr")
  p <- ggraph(lay) +
    geom_edge_link(alpha = .25, aes(color = mod)) +
    geom_node_point(size = 2) +
    scale_edge_colour_discrete(name = "Edge")

  p1 <- p + geom_node_point(size = 2, aes(color = clus)) +
    scale_colour_discrete(name = "Similarity Clustering")

  p2 <- p + geom_node_point(size = 2, aes(color = mem)) +
    scale_colour_discrete(name = "Network Community")

  res <- list()
  res$plot_net_1 <- p1 # cluster by hclust
  res$plot_net_2 <- p2 # cluster by community detection
  return(res)
}

## =======================================================================
## wl-23-09-2020, Wed: Network analysis using ionomics and symbolic data.
gene_network <- function(data,
                         data_symb, 
                         method_simil = "correlation",
                         thres_simil = 0.6,
                         thres_clus = 10) {

  data <- as.data.frame(data)
  rownames(data) <- data[, 1]

  ## clustering
  clust <- gene_clus(data_symb[, -1], thres_clus = thres_clus)

  ## update data
  data_symb$cluster <- clust$clus
  index <- clust$idx
  df_sub <- clust$tab_sub
  df_symb <- data_symb[clust$idx, ]
  
  ## cluster labels with accumulation/decumulation info
  func <- function(x, ...) {
    res <- NULL
    x <- x %>% dplyr::select(where(is.numeric))
    if (any(colSums(x == 1) > 0)) {
      res <- paste0(names(which(colSums(x == 1) > 0)), "(+)")
    }
    if (any(colSums(x == -1) > 0)) {
      res <- c(res, paste0(names(which(colSums(x == -1) > 0)), "(-)"))
    }
    res <- paste(res, collapse = ", ")
    res <- as.data.frame(res)
    return(res)
  }

  lab <- df_symb %>% 
    group_by(cluster) %>% 
    group_modify(func) %>% 
    dplyr::rename(label = res) 

  df_sub <- cbind(df_sub, label = lab[, 2])

  label <- sapply(df_symb$cluster, function(x) {
    tmp <- df_sub[df_sub$cluster == x, ]
    res <- paste0("C_", tmp[1], " (", tmp[2], " genes)")
    res <- paste(res, tmp[3], sep = ": ")
  })
  df_symb$label <- label

  ## get similarity of data set for network analysis
  sim <- df_simil(data[, -1], method = method_simil)
  diag(sim) <- 0
  sim <- sim[index, index]

  sim <- sym2long(sim, tri = "upper")
  idx <- abs(sim$var) >= thres_simil
  sim <- sim[idx, ]

  node <- data.frame(name = df_symb[, 1], label = label)

  g <- graph_from_data_frame(d = sim, directed = F, vertices = node)
  g <- delete_vertices(g, which(igraph::degree(g) < 1))

  E(g)$color <- ifelse(E(g)$var > 0, "pos", "neg")

  com <- cluster_louvain(g, weights = NULL)
  mem <- as.factor(membership(com))

  lay <- create_layout(g, layout = "fr")
  p <- ggraph(lay) + ## theme_graph() +
    geom_edge_link(alpha = .25, aes(color = color)) +
    geom_node_point() +
    scale_edge_colour_discrete(name = "Edge") ## edge legend title

  ## control node color by hclust
  p1 <- p +
    geom_node_point(size = 2, aes(color = label)) +
    scale_colour_discrete(name = "Similarity Clustering") ## node legend title

  ## control node color by community detection
  p2 <- p +
    geom_node_point(size = 2, aes(color = mem)) +
    scale_colour_discrete(name = "Network Community")

  ## betweenness
  btw <- betweenness(g)
  btw <- log(btw + 1)

  ## extract the same NAME
  nam <- names(btw)
  nam_1 <- as.character(data[, 1])
  idx <- nam_1 %in% nam
  mat <- data[idx, ]

  ## impact
  imp <- apply(mat[, -1], 1, norm, type = "2")

  ## get position
  t_imp <- quantile(imp, .75)
  t_btw <- quantile(btw, .75)
  position <-
    ifelse((imp < t_imp) & (btw < t_btw),
      "Low impact, low betweenness",
      ifelse((imp < t_imp) & (btw > t_btw),
        "Low impact, high betweenness",
        ifelse((imp > t_imp) & (btw < t_btw),
          "High impact, low betweenness",
          "High impact, high betweenness"
        )
      )
    )

  imp_btw <- data.frame(
    var_id = nam,
    impact = round(imp, 3),
    log.betweenness = round(btw, 3),
    position = position
  )

  imp_btw_p1 <-
    ggplot(data = imp_btw, aes(x = impact, y = log.betweenness)) +
    geom_point(aes(col = position), alpha = .3, size = 3) +
    theme(legend.position = "bottom") +
    guides(colour = guide_legend(nrow = 2)) +
    theme(legend.title = element_blank()) +
    geom_vline(xintercept = t_imp, linetype = "dashed") +
    geom_hline(yintercept = t_btw, linetype = "dashed") +
    xlab("Impact") +
    ylab("Log(betweenness+1)")

  imp_btw_p2 <-
    imp_btw_p1 +
    geom_text_repel(data = imp_btw, aes(label = var_id), size = 3.5)

  ## get impact-betweenness tables
  imp_btw <- cbind(imp_btw, cluster = V(g)$label)

  imp_btw_tab <- with(imp_btw, data.frame(table(cluster, position)))
  names(imp_btw_tab) <- c("Cluster", "Position", "nGenes")
  imp_btw_tab <- dplyr::arrange(imp_btw_tab, desc(nGenes))

  res <- list()
  res$plot_net_1 <- p1
  res$plot_net_2 <- p2
  res$plot_impact_betweenness_1 <- imp_btw_p1 # plot impact betweenees
  res$plot_impact_betweenness_2 <- imp_btw_p2 # plot impact betweenees
  res$impact_betweenness <- imp_btw # impact betweenees data
  return(res)
}



#' 
#' ## Data preparation
#' 
#' This document explains how to perform network and enrichment analysis on
#' the processed data set in publication:
#' 
#' Yu, D., Danku, J.M.C., Baxter, I. et al. **High-resolution genome-wide scan
#' of genes, gene-networks and cellular systems impacting the yeast ionome**.
#' BMC Genomics 13, 623 (2012). (https://doi.org/10.1186/1471-2164-13-623)
#' 
#' The authors identified 1065 strains with an altered ionome, including 584
#' haploid and 35 diploid deletion strains, and 446 over expression strains.
#' 
#' To explore the ionomics pipeline, we'll use the ionomics data set from the
#' paper:
#' 
## -----------------------------------------------------------------------------
## identified 584 haploid
dat <- read_csv("paper_ko.csv")
dim(dat)


#' 
#' This data has missing values. We use a simple missing value filling function
#' to impute the missing value:
#' 
## -----------------------------------------------------------------------------
## missing valuee filling with mean
dat <- dat %>% 
  mutate(across(where(is.numeric), function(x) {
    m <- mean(x, na.rm = TRUE)
    x[is.na(x)] <- m
    x
  }))
dat

#' 
#' To symbolise the data set, we use local FDR to estimate the symbolisation
#' threshold:
## -----------------------------------------------------------------------------
res <- locfdr_filter(t(dat[, -1]), plot = 1)
res$thres
## dat <- dat[res$idx, , drop = FALSE]
## symbolise data
dat_sym <- dat %>% 
  mutate(across(where(is.numeric), ~ dat_symb(., thres = res$thres)))

#' 
#' Some of ionomics data and symbolic data are like:
#' 
## -----------------------------------------------------------------------------
dat %>% sample_n(10) %>%
  kable(caption = 'Ionomics data', digits = 2, booktabs = T) %>% 
  kable_styling(full_width = F, font_size = 10,
                latex_options = c("striped", "scale_down"))

dat_sym %>% sample_n(10) %>%
  kable(caption = 'Symbolic data', booktabs = T) %>%
  kable_styling(full_width = F, font_size = 10,
                latex_options = c("striped", "scale_down"))

#' 
#' The data has been processed and the gene are significant. They are ready for
#' network and enrichment analysis.
#' 
#' ## Clustering ##
#' 
#' The gene network and enrichment analysis are based on gene data clustering.
#' The available similarity measures are *correlation*, *correlation*,
#' *cosine*,*eJaccard* and *Mahalanobis*. The last one is converted from
#' distance and will take a considerable time to run. The method for clustering
#' are from R base's `hclust` such as *complete* and *average*.
#' 
#' One of example: 
## -----------------------------------------------------------------------------
clus <- simil_clus(x = dat[, -1], method_simil = "correlation", 
                   method_hclus = "ave", thres_simil = 0.6, thres_clus = 5)
names(clus)

#' 
#' The cluster and its number of elements(larger than `thres_clus`) are: 
#' 
## -----------------------------------------------------------------------------
clus$tab_sub

#' 
#' The function `gene_clus` is also used for clustering but only for symbolic
#' data:
#' 
## -----------------------------------------------------------------------------
clus_1 <- gene_clus(x = dat_sym[, -1], thres_clus = 5)
names(clus_1)
clus_1$tab_sub

#' 
#' \clearpage
#' 
#' ## Gene network
#' 
#' Here we set up the parameters for network and enrichment analysis. The
#' similarity measure method is one of *correlation*, *cosine*, *eJaccard* and
#' *Mahalanobis*. The hierarchical clustering method is from `hclust`. This
#' should be (an unambiguous abbreviation of) one of *ward.D*, *ward.D2*,
#' *single*, *complete*, *average*, *mcquitty*, *median*  or *centroid*.
#' 
## -----------------------------------------------------------------------------
method_simil <- "correlation"
method_hclus <- "ave"
thres_simil  <- 0.70
thres_clus   <- 10

#' 
#' ### Network analysis with ionomics data only 
#' 
#' The gene network is based on the pre-processed data and use the similarity
#' measure results to build up the network. 
#' 
## -----------------------------------------------------------------------------
net <- gene_net(data = dat,
                method_simil = method_simil,
                method_hclus = method_hclus,
                thres_simil = thres_simil,
                thres_clus = thres_clus)

#' 
#' The node colours are indicated by the similarity measures and the network
#' community detection, i.e. clustering.
#' 
## ----gene-net, fig.cap = "Netwok analysis based on ionomics data"-------------
net$plot_net_1
net$plot_net_2

#' 
#' \clearpage
#' 
#' ### Network analysis with ionomics data and symbolic data
#' 
#' `gene_network` need both ionomics data and its symbolic data. The later is
#' used for hierarchical clustering while the former is used for correlation
#' analysis. The results of two analysis are used for network build up.
#' 
## ----gene-network, fig.cap = "Netwok analysis based on symbolic data"---------
net_1 <- gene_network(data = dat, data_symb = dat_sym,
                      method_simil = method_simil,
                      thres_simil = thres_simil,
                      thres_clus = thres_clus)

net_1$plot_net_1
net_1$plot_net_2
net_1$plot_impact_betweenness_1
net_1$plot_impact_betweenness_2
head(net_1$impact_betweenness)

#' 
#' \clearpage
#' 
#' ## Enrichment analysis
#' 
#' The results of hierarchical clustering either on ionomics data or on
#' symbolic data are for enrichment analysis.
#' 
#' ### Enrichment analysis on ionomics data
#' 
#' The KEGG enrichment analysis:
#' 
## -----------------------------------------------------------------------------
kegg <- kegg_enrich(data = dat, 
                    method_simil = method_simil, 
                    method_hclus = method_hclus, 
                    thres_simil = thres_simil,
                    thres_clus = thres_clus,
                    pval = 0.05, 
                    annot_pkg = "org.Sc.sgd.db",
                    is_symb = F)
## kegg
kegg %>% 
  kable(caption = 'KEGG enrichment analysis: ionomics data',
        digits = 3, booktabs = T) %>%
  kable_styling(full_width = F, font_size = 10,
                latex_options = c("striped", "scale_down"))

#' 
#' Note that there can be none results for KRGG  enrichment analysis. Change
#' arguments such as `thres_clus` as appropriate.
#' 
#' The GO Terms enrichment analysis:
#' 
## -----------------------------------------------------------------------------
go <- go_enrich(data = dat,
                method_simil, method_hclus, thres_simil, thres_clus,
                pval = 0.05, ont = "BP", annot_pkg =  "org.Sc.sgd.db",
                is_symb = F)
## go
go %>% head() %>% 
  kable(caption = 'GO Terms enrichment analysis: ionomics data',
        digits = 3, booktabs = T) %>%
  kable_styling(full_width = F, font_size = 10,
                latex_options = c("striped", "scale_down"))

#' 
#' \clearpage
#' 
#' ### Enrichment analysis on symbolic data
#' 
#' The KEGG enrichment analysis:
#' 
## -----------------------------------------------------------------------------
kegg_1 <- kegg_enrich(data = dat_sym, 
                      method_simil = method_simil, 
                      method_hclus = method_hclus, 
                      thres_simil = thres_simil,
                      thres_clus = thres_clus,
                      pval = 0.05, 
                      annot_pkg = "org.Sc.sgd.db",
                      is_symb = T)
## kegg_1
kegg_1 %>%
  kable(caption = 'KEGG enrichment analysis: Symbolic data', digits = 3,
        booktabs = T) %>%
  kable_styling(full_width = F, font_size = 10,
                latex_options = c("striped", "scale_down"))

#' 
#' The GO Terms enrichment analysis:
#' 
## -----------------------------------------------------------------------------
go_1 <- go_enrich(data = dat_sym,
                  method_simil, method_hclus, thres_simil, thres_clus,
                  pval = 0.05, ont = "BP", annot_pkg =  "org.Sc.sgd.db",
                  is_symb = T)
## go_1
go_1 %>% head() %>%
  kable(caption = 'GO Terms enrichment analysis: Symbolic data',
        digits = 3, booktabs = T) %>%
  kable_styling(full_width = F, font_size = 10,
                latex_options = c("striped", "scale_down"))

