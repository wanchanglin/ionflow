#' wl-02-07-2020, Thu: Put package files together
#' wl-06-07-2020, Mon: 
#'  - find '<<-'. Should be '<-'
#'  - global data sets: data_GOslim and data_ORF2KEGG are used in
#'    GeneClustering. Should change

#' ==== Pre-processing ====
#'
#'Pre-processing tools for Ionomics data
#' @param \code{data} a data frame of ion's concentrations (columns) for various
#'   gene's knockout (rows)
#' @param \code{stdev} an optional pre-defined vector of standard deviations for
#'   each ion
#'
#' @return \code{stats.raw_data} statistics for the raw data ion concentration
#'   values
#' @return \code{stats.outliers} statistics (raw count and percentage of
#'   outliers) from the raw data for each ion
#' @return \code{stats.median_batch_corrected_data} statistics for the median
#'   batch corrected data
#' @return \code{stats.standardised_data} statistics for standardised data
#' @return \code{plot.logConcentration_z_scores} plot of the symbolised ion
#'   concentration profiles of the knockouts
#' @return \code{plot.logConcentration_by_batch} plot by batch of the log-ion
#'   concentration values after outlier removal and median batch correction
#' @return \code{dataR.long} a processed data frame of ion's concentrations with
#'   batch-replicated gene's knockout (long format)
#' @return \code{data.long} a processed data frame of ion's concentrations for
#'   each gene's knockout (long format)
#' @return \code{data.wide} a processed data frame of ion's concentrations
#'   (columns) for each gene's knockout (rows) (wide format)
#' @return \code{data.wide_Symb} a processed data frame of symbolised (values
#'   -1,0, or 1) ion's concentration profiles (columns) for each gene's knockout
#'   (rows) (wide format)
#'
#' @examples \code{
#' library(IonFlow)
#'
#' ### Run PreProcessing function
#' pre_proc <- PreProcessing(data=IonData,stdev=pre_defined_sd)
#'
#' # stats
#' pre_proc$stats.raw_data
#' pre_proc$stats.outliers
#' pre_proc$stats.median_batch_corrected_data
#' pre_proc$stats.standardised_data
#'
#' # plots
#' pre_proc$plot.logConcentration_by_batch
#' pre_proc$plot.logConcentration_z_scores
#'
#' # data
#' head(pre_proc$dataR.long)
#' head(pre_proc$data.long)
#' head(pre_proc$data.wide)
#' head(pre_proc$data.wide_Symb)
#' }
#' @export
PreProcessing = function(data=NULL,stdev=NULL) {

  #### -------------------> Import data
  list.s <- list()
  data2 <- data[,-c(1,2)]
  for (i in 1:14){
    list.s[[i]] <- c(names(data2)[i],round(summary(data2[,i]),3),
                     round(var(data2[,i]),3))
  }
  df.s <- data.frame(do.call(rbind,list.s))
  names(df.s) <-  c('Ion','Min','1st Quartile','Median','Mean', '3rd Quartile', 
                    'Max','Variance' )
  
  #' wl-06-07-2020, Mon: should simplify something like this:
  #' tmp <- as.data.frame(t(sapply(data2, summary)))
  
  #### -------------------> Outlier detection
  data$id <- row.names(data)
  data_long <- tidyr::gather(data,Ion,Concentration, Ca:Zn, factor_key=TRUE)

  #' wl-06-07-2020, Mon: BUG. It's dangerous to use 'levels'. It is only for
  #' factor.  Make vectors as factors before using levels function. Don't
  #' presume they are factors. Check them using:
  #' lapply(data_long, class)
  data_long$Knockout <- as.factor(data_long$Knockout)
  data_long$Ion <- as.factor(data_long$Ion)

  for (i in 1:length(levels(data_long$Ion))){
    data_long_sub <- data_long[data_long$Ion==levels(data_long$Ion)[i],]
    lowerq = quantile(data_long_sub$Concentration,na.rm=T)[2]
    upperq = quantile(data_long_sub$Concentration,na.rm=T)[4]
    iqr = upperq - lowerq
    extreme.t.upper = (iqr * 3) + upperq
    extreme.t.lower = lowerq - (iqr * 3)
    data_long[data_long$Ion==levels(data_long$Ion)[i],'Outlier'] <- 
      ifelse((data_long_sub$Concentration > extreme.t.upper || data_long_sub$Concentration < extreme.t.lower),1,0)

    #' wl-06-07-2020, Mon: potential BUG. Should use '||'. Should change '&'
    #' as '&&' as well in the rest of this script.
  }
  df_outlier <- data.frame(cbind(levels(data_long$Ion),
                                 table(data_long$Ion,data_long$Outlier),
                                 round(table(data_long$Ion,data_long$Outlier)[,2]/dim(data_long)[1]*100,2)))
  rownames(df_outlier) <- c()
  colnames(df_outlier) <- c('Ion','no outlier','outlier','outlier(%)')
  data_long_clean <- data_long[data_long$Outlier<1,]

  #### -------------------> Median batch correction
  data_long_clean$logConcentration <- log(data_long_clean$Concentration)
  med <- matrix(0,length(table(data$Batch_ID)),length(levels(data_long$Ion)))

  datasets <- list()
  j <- 1
  while (j <= length(levels(data_long_clean$Ion))){
    data_long_clean_sub <- data_long_clean[data_long_clean$Ion==levels(data_long_clean$Ion)[j],]
    for (i in 1:max(data_long_clean_sub$Batch_ID)){
      med[i,j] <- median(data_long_clean_sub$logConcentration[data_long_clean_sub$Batch_ID==i])
      data_long_clean_sub$logConcentration_corr[data_long_clean_sub$Batch_ID==i] <- 
        data_long_clean_sub$logConcentration[data_long_clean_sub$Batch_ID==i] - med[i,j]
    }
    datasets[[j]] <- data_long_clean_sub

    #' wl-06-07-2020, Mon: move out the following line
    #' data_long_clean_scaled <- do.call(rbind, datasets)
    j <- j+1
  }

  #' wl-06-07-2020, Mon: combine 
  data_long_clean_scaled <- do.call(rbind, datasets)

  p1 <- ggplot2::ggplot(data = data_long_clean_scaled, 
                        ggplot2::aes(x = factor(Batch_ID), 
                                     y = logConcentration_corr, 
                                     col=factor(Batch_ID))) +
  geom_point(shape=1) +
  facet_wrap(~Ion) +
  xlab("Batch.ID") +
  ylab("log(Concentration) (ppm)") +
  theme(legend.position="none") +
  theme(axis.text.x=element_blank())

  list.mbc <- list()

  #' for (i in 1:max(length(levels(data_long_clean_scaled$Ion)))){
  for (i in 1:length(levels(data_long_clean_scaled$Ion))){
    list.mbc[[i]] <- c(levels(data_long_clean_scaled$Ion)[i],
                       round(summary(datasets[[i]]$logConcentration_corr),3),
                       round(var(datasets[[i]]$logConcentration_corr),3))
  }
  df.mbc <- data.frame(do.call(rbind,list.mbc))
  names(df.mbc) <- c('Ion','Min','1st Quartile','Median','Mean', '3rd Quartile', 'Max','Variance' )

  #### -------------------> Standardisation
  if (is.null(stdev)) {
    sds <- rep(NA,length(levels(data_long_clean_scaled$Ion)))
    for (i in 1:length(levels(data_long_clean_scaled$Ion))){
      sds[i] <- as.numeric(sd(data_long_clean_scaled$logConcentration_corr[data_long_clean_scaled$Ion==levels(data_long_clean_scaled$Ion)[i]])) # Ions' sd of logConcentration_corr
    }
  } else {
    sds=as.numeric(as.vector(stdev[,1]))
  }

  datasets2 <- list()
  j <- 1
  while (j <= length(levels(data_long_clean_scaled$Knockout))){
    data_long_clean_scaled_sub <- data_long_clean_scaled[data_long_clean_scaled$Knockout==levels(data_long_clean_scaled$Knockout)[j],]
    for (i in 1:length(levels(data_long_clean_scaled$Ion))){
      data_long_clean_scaled_sub$logConcentration_corr_norm[data_long_clean_scaled_sub$Ion==levels(data_long_clean_scaled_sub$Ion)[i]] <-
        data_long_clean_scaled_sub$logConcentration_corr[data_long_clean_scaled_sub$Ion==levels(data_long_clean_scaled_sub$Ion)[i]] / sds[i]
    }
    datasets2[[j]] <- data_long_clean_scaled_sub

    #' wl-06-07-2020, Mon: move out the following line
    #' data_long_clean_scaled_norm <- do.call(rbind, datasets2)
    j <- j+1
  }
  
  #' wl-06-07-2020, Mon: combine 
  data_long_clean_scaled_norm <- do.call(rbind, datasets2)

  list.mbc2 <- list()

  #' for (i in 1:max(length(levels(data_long_clean_scaled_norm$Ion)))){
  for (i in 1:length(levels(data_long_clean_scaled_norm$Ion))){
    list.mbc2[[i]] <- c(levels(data_long_clean_scaled_norm$Ion)[i],
                        round(summary(datasets2[[i]]$logConcentration_corr_norm),3),
                        round(var(datasets2[[i]]$logConcentration_corr_norm),3))
  }

  df.mbc2 <- data.frame(do.call(rbind,list.mbc2))
  names(df.mbc2) <- c('Ion','Min','1st Quartile','Median','Mean', '3rd Quartile', 'Max','Variance' )

  #### -------------------> Symbolization
  data_long_clean_scaled_norm$Symb <- 
    ifelse((data_long_clean_scaled_norm$logConcentration_corr_norm > -3) && (data_long_clean_scaled_norm$logConcentration_corr_norm< 3), 
           0, ifelse(data_long_clean_scaled_norm$logConcentration_corr_norm>=3,1,-1))

  #### -------------------> Aggregation of the batch replicas
  data_long_clean_scaled_norm_unique <- 
    data.frame(aggregate(. ~ Knockout*Ion, data_long_clean_scaled_norm[,c('Knockout','Ion','logConcentration_corr_norm','Symb')], median))

  data_long_clean_scaled_norm_unique$Symb <- 
    ifelse((data_long_clean_scaled_norm_unique$Symb<0.5) && (data_long_clean_scaled_norm_unique$Symb>-0.5), 
           0, ifelse(data_long_clean_scaled_norm_unique$Symb>=0.5,1,-1))

  data_wide_clean_scaled_norm_unique <- 
    reshape2::dcast(data_long_clean_scaled_norm_unique, Knockout~ Ion, value.var="logConcentration_corr_norm")

  data_wide_clean_scaled_norm_unique_Symb <- 
    reshape2::dcast(data_long_clean_scaled_norm_unique, Knockout~ Ion, value.var="Symb")

  p2 <- ggplot2::ggplot(data = data_long_clean_scaled_norm_unique, ggplot2::aes(x = logConcentration_corr_norm)) +
  geom_histogram(binwidth=.1) +
  facet_wrap(~Ion) +
  xlab("log(Concentration) (z-score)") +
  ylab("Frequency") +
  geom_vline(xintercept=c(-3,3),col='red')

  #### -------------------> Output
  res <- list()
  class(res) = "PreProcessing"

  res$stats.raw_data <- df.s # raw data
  res$stats.outliers <- df_outlier # outliers
  res$stats.median_batch_corrected_data <- df.mbc # median batch corrected data
  res$stats.standardised_data <- df.mbc2 # standardised data

  res$dataR.long <- data_long_clean_scaled_norm
  res$data.long <- data_long_clean_scaled_norm_unique
  res$data.wide <- data_wide_clean_scaled_norm_unique
  res$data.wide_Symb <- data_wide_clean_scaled_norm_unique_Symb

  res$plot.logConcentration_by_batch <- p1
  res$plot.logConcentration_z_scores <- p2
  return(res)
}

#' ==== Exploratory Analysis ====
#'
#'Exploratory tools for Ionomics data
#' @param \code{data} a processed data frame of ion's concentrations (columns)
#'   for each gene's knockout (rows) (wide format) e.g. \code{data.wide} from
#'   the \code{PreProcessing_fn}
#'
#' @return \code{plot.Pearson_correlation} plot of Pearson's correlations across
#'   ions
#' @return \code{plot.PCA_Individual} plot of ion's profiles in a
#'   lower-dimensional space from PCA
#' @return \code{plot.heatmap} plot of ion's concentrations for each gene's
#'   knockout with added clustering dendogram
#' @return \code{plot.pairwise_correlation_map} plot of pairwise correlation
#'   coefficients across ions
#' @return \code{plot.regularized_partial_correlation_network} plot of partial
#'   correlation coefficients for pairs of ions after conditioning for the
#'   remaining ions. Red for negative, green for positive partial correlations.
#'   No connections drawn when partial correlations equals 0
#' @return \code{data.PCA_loadings} a data frame of weights of each of the
#'   original gene knockouts also known as loading vectors associated to each PC
#'
#' @examples \code{
#' library(IonFlow)
#'
#' ### Run PreProcessing function
#' pre_proc <- PreProcessing(data=IonData,stdev=pre_defined_sd)
#'
#' ### Run ExploratoryAnalysis function
#' exp_anal <- ExploratoryAnalysis(data=pre_proc$data.wide)
#'
#' # plots
#' exp_anal$plot.Pearson_correlation
#' exp_anal$plot.PCA_Individual
#' exp_anal$plot.heatmap
#' exp_anal$plot.pairwise_correlation_map
#' exp_anal$plot.regularized_partial_correlation_network
#'
#' # data
#' head(exp_anal$data.PCA_loadings)
#' }
#' @export
ExploratoryAnalysis = function(data=NULL) {

  #### -------------------> Correlation
  col3 <- colorRampPalette(c("steelblue4","white","firebrick"))
  corrplot.mixed(cor(data[,-1], use = "complete.obs"), number.cex = .7,
                 lower.col = 'black', upper.col = col3(100))

  p_corr <- recordPlot() 

  #### -------------------> PCA
  pca.X <- mixOmics::pca(t(data[,-1]),center=T,scale=F)

  loadings_PC1 <- data.frame(mixOmics::selectVar(pca.X, comp=1)$value)
  row.names(loadings_PC1) <- data$Knockout[order(loadings_PC1)]
  loadings_PC2 <- data.frame(mixOmics::selectVar(pca.X, comp=2)$value)
  row.names(loadings_PC2) <- data$Knockout[order(loadings_PC2)]
  PCA_loadings <- data.frame(cbind(loadings_PC1,loadings_PC2))
  names(PCA_loadings) <- c('PC1','PC2')

  # Individual factor map
  dtp <- data.frame('Ion' = row.names(data.frame(pca.X$variates)), pca.X$variates)
  names(dtp) <- c("Ion","PC1","PC2")

  pca_plot <- ggplot(data = dtp, aes(x = PC1, y = PC2)) +
    geom_point(color='steelblue', size=3, alpha=0.4) +
    geom_text_repel(ggplot2::aes(label = row.names(data.frame(pca.X$variates))), size=4) +
    theme_bw() +
    theme(legend.position = "none") +
    xlab(paste("PC1: ",round(pca.X$explained_variance[1]*100,2),"% expl. variance", sep='')) +
    ylab(paste("PC2: ",round(pca.X$explained_variance[2]*100,2),"% expl. variance", sep='')) +
    labs(title = "PCA")

  #### -------------------> HEATMAP
  pheatmap(data[,-1], show_rownames=F, cluster_cols=T, cluster_rows=T,
           legend = T, fontsize = 15, clustering_method = 'ward.D',
           scale="row")

  pheat <- recordPlot()

  #### -------------------> PAIRWISE CORRELATION MAP
  col <- colorRampPalette(c("skyblue4", "white", "plum4"))(20)

  corrI <- cor(na.omit(data[,-1]))
  heatmap(x = corrI, col = col, symm = TRUE, cexRow=1.4, cexCol=1.4, main = '')

  pcm <- recordPlot()

  #### -------------------> Regularized partial correlation network MAP
  cad <- cor_auto(data[,-1]) # wl-06-07-2020, Mon: from package lavaan
  suppressWarnings(qgraph(cad, graph = "glasso", layout = "spring",
         tuning = 0.25,sampleSize = nrow(data[,-1])))

  Graph_lasso <- recordPlot()

  #### -------------------> Output
  res <- list()
  class(res) = "ExploratoryAnalysis"
  res$plot.Pearson_correlation <- p_corr
  res$plot.PCA_Individual <- pca_plot
  res$data.PCA_loadings <- PCA_loadings
  res$plot.heatmap <- pheat
  res$plot.pairwise_correlation_map <- pcm
  res$plot.regularized_partial_correlation_network <- Graph_lasso

  return(res)
}

#' ==== Clustering ====
#'
#'Clustering and enrichment tools for Ionomics data
#' @param \code{data} a processed data frame of ion's concentrations (columns)
#'   for each gene's knockout (rows) (wide format) e.g. \code{data.wide} from
#'   the \code{PreProcessing_fn}
#' @param \code{data_Symb} a processed data frame of symbolised (values -1,0, or
#'   1) ion's concentration profiles (columns) for each gene's knockout (rows)
#'   (wide format) e.g. \code{data.wide_Symb} from the \code{PreProcessing_fn}
#'
#' @return \code{stats.clusters statistics} of the gene's count by cluster
#' @return \code{stats.Kegg_Goslim_annotation} statistics from the GO Slim
#'   annotation
#' @return \code{stats.Goterms_enrichment} statistics from the GO terms
#'   enrichment
#'
#' @return \code{plot.profiles} plot of the gene's profiles for each cluster
#'
#' @examples \code{
#' library(IonFlow)
#'
#' ### Run PreProcessing function
#' pre_proc <- PreProcessing(data=IonData,stdev=pre_defined_sd)
#'
#' ### Run GeneClustering function
#' gene_clust <- GeneClustering(data=pre_proc$data.wide,
#'   data_Symb=pre_proc$data.wide_Symb)
#'
#' # stats
#' gene_clust$stats.clusters
#' gene_clust$stats.Kegg_Goslim_annotation
#' gene_clust$stats.Goterms_enrichment
#'
#' # plots
#' gene_clust$plot.profiles
#' }
#' @export
GeneClustering = function(data=NULL, data_Symb=NULL) {

  #### -------------------> Define clusters
  res.dist <- dist(data_Symb[,-1], method = "manhattan")
  res.hc <- hclust(d = res.dist, method = "single")
  data_Symb$cluster <- cutree(res.hc, h = 0) # distance 0

  #### -------------------> Subset cluster with more than 10 genes
  df <- as.data.frame(table(data_Symb$cluster)); 
  names(df) <- c('cluster', 'nGenes')
  df_sub <- df[df$nGenes>10,]
  rownames(df_sub) <- c()

  #### -------------------> Plot clustering
  df.tmp_long = list()
  for(i in 1:dim(df_sub)[1]){
    tmp_wide <- data[data_Symb$cluster==df_sub$cluster[i],]
    tmp_long <- gather(tmp_wide, Ion, logConcentration_corr_norm,  Ca:Zn, factor_key=TRUE)
    tmp_long$Cluster_number <- rep(df_sub$cluster[i],nrow(tmp_long))
    label <- paste("Cluster",df_sub$cluster[i],paste("(",df_sub$nGenes[i], " genes)", sep=''), sep=' ')
    tmp_long$Cluster <- rep(label,nrow(tmp_long))
    df.tmp_long[[i]] <- tmp_long
  }
  df.tmp_long_all <- data.frame(do.call(rbind,df.tmp_long))

  p1 <- ggplot2::ggplot(data = df.tmp_long_all, ggplot2::aes(x = Ion, y = logConcentration_corr_norm)) +
    facet_wrap(~Cluster) +
    geom_line(ggplot2::aes(group = Knockout)) +
    stat_summary(fun.data = "mean_se", color = "red")+
    labs(x = "", y = "") +
    theme(legend.position="none") +
    theme(axis.text.x=element_text(angle=90, hjust=1),axis.text=element_text(size=10))

  #### -------------------> KEGG AND GO SLIM ANNOTATION
  p.data_list = list()
  for(i in 1:dim(df_sub)[1]){
    inputGeneSet <- data_Symb$Knockout[data_Symb$cluster==df_sub$cluster[i]]

    N <- as.numeric(length(inputGeneSet))
    p.data <- data_GOslim %>%
      dplyr::mutate(Ontology = setNames(c("Biological process","Cellular component","Molecular function"),c("P","C","F"))[as.character(Ontology)]) %>%
      dplyr::filter(ORFs %in% as.character(inputGeneSet)) %>%
      dplyr::group_by(GOslim,Ontology) %>%
      dplyr::filter(GOslim != "other") %>%
      dplyr::rename(Term = GOslim) %>%
      dplyr::summarise(Count = n()) %>%
      dplyr::mutate(Percent = Count/N*100) %>%
      dplyr::bind_rows(data_ORF2KEGG %>%
                         dplyr::filter(ORF %in% as.character(inputGeneSet)) %>%
                         dplyr::group_by(KEGGID, Pathway) %>%
                         dplyr::summarise(Count = dplyr::n()) %>%
                         dplyr::mutate(Ontology = "KEGG") %>%
                         dplyr::rename(Term = Pathway) %>%
                         dplyr::ungroup() %>%
                         dplyr::filter(KEGGID != "01100") %>%
                         dplyr::select(-KEGGID)%>%
                         dplyr::mutate(Percent = Count/N*100)) %>%
      dplyr::filter(!Term %in% c("molecular_function", "biological_process", "cellular_component"))
    p.data_list[[i]] = p.data
  }

  df.data_list = label_list <- list()
  for(i in 1:dim(df_sub)[1]){
    p.data_list_sub <- as.data.frame(p.data_list[[i]])
    p.data_list_sub$Cluster_number <- rep(df_sub$cluster[i],nrow(p.data_list_sub))
    label <- paste("Cluster",df_sub$cluster[i],paste("(",df_sub$nGenes[i], " genes)", sep=''), sep=' ')
    p.data_list_sub$Cluster <- rep(label,nrow(p.data_list_sub))
    sub_data <- 
      p.data_list_sub[(p.data_list_sub$Percent>5) && 
                      (p.data_list_sub$Ontology %in% c("Biological process","Cellular component","Molecular function")),]
    df.data_list[[i]] <- sub_data
    label_list[[i]] <- label
  }

  Kegg_Goslim_annotation.list_by_clust <- lapply(df.data_list,"[", c('Term','Ontology','Count','Percent'))
  names(Kegg_Goslim_annotation.list_by_clust) <- label_list

  #### -------------------> GO TERMS ENRICHMENT
  universeGenes=as.character(data_Symb$Knockout) # universo dei dati sperimentali

  p.data_list2 = label_list2 = list()
  for(i in 1:dim(df_sub)[1]){
    label <- paste("Cluster",df_sub$cluster[i],paste("(",df_sub$nGenes[i], " genes)", sep=''), sep=' ')
    inputGeneSet <- data_Symb$Knockout[data_Symb$cluster==df_sub$cluster[i]]

    ont=c("BP","MF","CC")
    results <- c()
    for(k in 1:3){
      params = new("GOHyperGParams",
                   geneIds=as.character(inputGeneSet),
                   universeGeneIds=universeGenes,
                   annotation="org.Sc.sgd.db",
                   categoryName="GO",
                   ontology=ont[k],
                   pvalueCutoff=0.05,
                   conditional=T,
                   testDirection="over")
      hgOver <- hyperGTest(params)
    }

    results <-
      rbind(results,
            cbind(setNames(dplyr::data_frame(ID=names(pvalues(hgOver)), Term = Term(ID),
                                             pvalues = pvalues(hgOver), oddsRatios = oddsRatios(hgOver),expectedCounts=expectedCounts(hgOver),
                                             geneCounts=geneCounts(hgOver),universeCounts=universeCounts(hgOver)),
                           c("GO_ID","Description","Pvalue","OddsRatio","ExpCount","Count","CountUniverse")),"Ontology"=ont[k])
            %>% dplyr::filter(Pvalue <= 0.05 && Count > 1))

    p.data_list2[[i]] <- results
    label_list2[[i]] <- label
  }
  Goterms_enrichment.list_by_clust <- lapply(p.data_list2,"[",-c(4,5,8))
  names(Goterms_enrichment.list_by_clust) <- label_list2

  #### -------------------> Output
  res <- list()
  class(res) = "GeneClustering"

  res$stats.clusters <- df_sub # selected clusters (nGenes>10)
  res$plot.profiles <- p1 # plot cluster profiles
  res$stats.Kegg_Goslim_annotation <- Kegg_Goslim_annotation.list_by_clust # KEGG AND GO SLIM ANNOTATION
  res$stats.Goterms_enrichment <- Goterms_enrichment.list_by_clust # GO TERMS ENRICHMENT

  return(res)
}

#' ==== Network ====
#'
#'Network tools for Ionomics data
#' @param \code{data} a processed data frame of ion's concentrations (columns)
#'   for each gene's knockout (rows) (wide format) e.g. \code{data.wide} from
#'   the \code{PreProcessing_fn}
#' @param \code{data_Symb} a processed data frame of symbolised (values -1,0, or
#'   1) ion's concentration profiles (columns) for each gene's knockout (rows)
#'   (wide format) e.g. \code{data.wide_Symb} from the \code{PreProcessing_fn}
#'
#' @return \code{stats.impact_betweeness} statistics of the impact betweenees
#'   data
#' @return \code{stats.impact_betweeness_by_cluster} statistics of the impact
#'   betweenees data by cluster
#'
#' @return \code{plot.pnet plot} of the gene's network
#' @return \code{plot.impact_betweenees} plot of the impact betweenees
#'
#' @examples \code{
#' library(IonFlow)
#'
#' ### Run PreProcessing function
#' pre_proc <- PreProcessing(data=IonData,stdev=pre_defined_sd)
#'
#' ### Run GeneNetwork function
#' gene_net <- GeneNetwork(data=pre_proc$data.wide,
#'   data_Symb=pre_proc$data.wide_Symb)
#'
#' # stats
#' gene_net$stats.impact_betweeness
#' gene_net$stats.impact_betweeness_by_cluster
#'
#' # plots
#' gene_net$plot.pnet
#' gene_net$plot.impact_betweenees
#' }
#' @export
GeneNetwork = function(data=NULL, data_Symb=NULL) {

  # Cluster of gene with same profile
  res.dist <- dist(data_Symb[,-1], method = "manhattan")
  res.hc <- hclust(d = res.dist, method = "single")
  data_Symb$cluster <- cutree(res.hc, h = 0) # distance 0

  df <- as.data.frame(table(data_Symb$cluster)); names(df) <- c('cluster', 'nGenes')
  for(i in 1:dim(df)[1]){
    tmp_wide <- data[data_Symb$cluster==df$cluster[i],]
    label <- paste("Cluster",df$cluster[i],paste("(",df$nGenes[i], " genes)", sep=''), sep=' ')
    data_Symb$Cluster[data_Symb$cluster==df$cluster[i]] <- rep(label,nrow(tmp_wide))
  }

  # Freq in each cluster
  Rk <- as.data.frame(table(data_Symb$Cluster))
  Rk <- Rk[with(Rk, order(-Rk$Freq)), ]

  # Compute gene cluster id (index)
  index = data_Symb$Cluster

  # Cluster 2 (largest cluster) contains genes with no phenotype hence not considered (input to 0)
  index[index==Rk$Var1[1]]=0
  # Smaller clusters (Freq<10) input to 0 as well (we consider only cluster with at least 10 genes)
  for(i in which(Rk$Freq<10)[1]:dim(Rk)[1]){
    index[index==Rk$Var1[i]]=0
  }
  # Apply the cluster filtering
  index=index>0

  # cluster labels with info of accumulation/decumulation of Ions (high/lower abundance)
  x=data_Symb$Cluster[index]
  ux=unique(x)
  df.symb <- data_Symb[data_Symb$Cluster %in% c(unique(x)),]
  df <- as.data.frame(table(df.symb$Cluster)); names(df) <- c('Cluster', 'nGenes')

  # Assign label
  labelsC2 <-labelsC <- list()
  name.label <- c()
  for(i in 1:length(ux)){
    sub_df.symb <- df.symb[(df.symb$Cluster==ux[i]), 2:15]#; rownames(sub_df.symb) <- df.symb[df.symb$cluster==ux[i],1]
    if (length(names(which(colSums(sub_df.symb == 1) > 0)))>0) {l1 <- paste0(names(which(colSums(sub_df.symb == 1) > 0)),"(+)")}
    if (length(names(which(colSums(sub_df.symb == -1) > 0)))>0) {l2 <- paste0(names(which(colSums(sub_df.symb == -1) > 0)),"(-)")}
    if (length(names(which(colSums(sub_df.symb == 1) > 0)))>0) {labelsC[[i]] <- l1}
    if (length(names(which(colSums(sub_df.symb == -1) > 0)))>0) {labelsC[[i]] <- l2}
    if ((length(names(which(colSums(sub_df.symb == 1) > 0)))>0) && (length(names(which(colSums(sub_df.symb == -1) > 0)))>0)) {labelsC[[i]] <- c(l1,l2)}
    labelsC2[[i]] <- do.call(paste, c(as.list(labelsC[[i]]), sep = ", "))
  }
  names(labelsC2) = ux

  ny <- paste(ux,unlist(labelsC2),sep = ": ")
  for(i in 1:length(ux)){
    sub_df.symb <- df.symb[(df.symb$Cluster==ux[i]),]
    df.symb$Label[(df.symb$Cluster==ux[i])] <- rep(ny[i],nrow(sub_df.symb))
  }

  x1=df.symb$Label
  ux1=unique(x1)
  cpy= rainbow(length(ux))
  names(cpy) = ux1

  # Compute empirical correlation matrix
  corrGenes <- cor(t(as.matrix(data[,-1])), method = "pearson", 
                   use = "pairwise.complete.obs")

  # Subset correlation matrix based on the cluster filtering
  A = corrGenes[index,index]
  # Diagonal value (1's) put to 0 to avoid showing edges from/to the same gene
  diag(A) <- 0
  # Subset correlation matrix based on threshold=0.6
  A=(A>0.6)
  A <- ifelse(A==TRUE,1,0)
  # Generate network
  Net = network::network(A, directed = FALSE)
  #h <- network.copy(Net)
  #summary(h)

  # Network plot

  gNEt <- function(dat=NULL ){
    GGally::ggnet2(dat, color = x1, color.legend = 'Label', palette = cpy,
                   edge.alpha = 0.5, size = 2,
                   legend.size = 10, legend.position = "right")
  }

  # Impact and betweenness
  df1 <- data[data_Symb$Cluster %in% ux,] # df of 269 genes, 10 clusters
  btw <- sna::betweenness(A)
  impact <- apply(df1[,-1],1, function(a) dist(rbind(a,rep(0,ncol(df1[-1]))))) # L2 norm

  df.res <- data.frame(
    Knockout = data$Knockout[index],
    impact = round(impact,3),
    betweenness = round(btw,3),
    log.betweenness = round(log(btw+1),3),
    pos = factor(ifelse((impact < quantile(impact,.75)) && (log(btw+1) < quantile(log(btw+1),.75)),1,
                        ifelse((impact < quantile(impact,.75)) && (log(btw+1) > quantile(log(btw+1),.75)),2,
                               ifelse((impact > quantile(impact,.75)) && (log(btw+1) < quantile(log(btw+1),.75)),3,4)))),
    pos.label = factor(ifelse((impact < quantile(impact,.75)) && (log(btw+1) < quantile(log(btw+1),.75)), 'Low impact, low betweenness',
                              ifelse((impact < quantile(impact,.75)) && (log(btw+1) > quantile(log(btw+1),.75)),'Low impact, high betweenness',
                                     ifelse((impact > quantile(impact,.75)) && (log(btw+1) < quantile(log(btw+1),.75)),'High impact, low betweenness','High impact, high betweenness')))))

  rownames(df.res) = data$Knockout[index]
  q1 <- row.names(subset(df.res, (impact < quantile(impact,.75)) && (log.betweenness < quantile(log.betweenness,.75))))
  q2 <- row.names(subset(df.res, (impact < quantile(impact,.75)) && (log.betweenness > quantile(log.betweenness,.75))))
  q3 <- row.names(subset(df.res, (impact > quantile(impact,.75)) && (log.betweenness < quantile(log.betweenness,.75))))
  q4 <- row.names(subset(df.res, (impact > quantile(impact,.75)) && (log.betweenness > quantile(log.betweenness,.75))))
  idx <- unique(c(sample(q1,6),sample(q2,6),sample(q3,6),sample(q4,6)))
  df.idx <- df.res[idx,]

  gp1 <-  ggplot2::ggplot(data=df.res,ggplot2::aes(x=impact,y=log.betweenness)) +
    geom_point(ggplot2::aes(col=pos.label), alpha=.3, size=3) +
    scale_color_manual(values=c("plum4","palegreen4","indianred","cornflowerblue")) +
    theme_linedraw() +
    theme_light() +
    theme(legend.position="bottom") +
    guides(colour = guide_legend(nrow = 2)) +
    theme(legend.title=element_blank()) +
    geom_text_repel(data=df.idx, ggplot2::aes(label=df.idx$Knockout), size=3.5) +
    geom_vline(xintercept=quantile(df.res$impact,.75),linetype="dashed") +
    geom_hline(yintercept=quantile(df.res$log.betweenness,.75),linetype="dashed") +
    xlab("Impact") +
    ylab("Log(betweenness+1)")

  rownames(df.res) <- c()
  df.res2 <- df.res[,-c(4,5)]
  names(df.res2) <- c('Knockout','Impact','Betweenness','Position')

  gene.cluster <- df.symb[,c(1,18)]
  names(gene.cluster) <- c('Knockout','Cluster')
  df.res3 <-merge(df.res2, gene.cluster, by="Knockout",all.x=TRUE)

  df.tab <- data.frame(table(df.res3$Cluster,df.res3$Position))

  #' wl-06-07-2020, Mon: add dplyr:: for clarify
  df.tab2 <- df.tab %>% dplyr::group_by(Var1) %>% top_n(1, Freq)
  names(df.tab2) <- c('Cluster','Position','nGenes')

  #### -------------------> Output
  res <- list()
  class(res) = "GeneNetwork"
  res$plot.pnet <- gNEt(dat=Net) # plot gene network
  res$plot.impact_betweenees <- gp1 # plot impact betweenees
  res$stats.impact_betweeness <- df.res3 # impact betweenees data
  res$stats.impact_betweeness_by_cluster <- df.tab2 # plot position by cluster

  return(res)
}

