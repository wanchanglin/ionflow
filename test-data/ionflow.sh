# wl-09-08-2020, Sun: Rscript test code for Linux. 
# wl-12-11-2020, Thu: Modification 

Rscript --vanilla ../ionflow.R \
  --ion_file "../test-data/iondata.tsv" \
  --var_id "1" \
  --batch_id "2" \
  --data_id "3" \
  --method_norm "median" \
  --batch_control "no" \
  --method_outliers "IQR" \
  --thres_outl 3.0 \
  --stand_method "std" \
  --std_file "../test-data/user_std.tsv" \
  --thres_symb 2.0 \
  --min_clust_size 10.0 \
  --thres_corr "0.75" \
  --method_corr "pearson" \
  --clus_id "2" \
  --pval 0.05 \
  --ont "BP" \
  --annot_pkg "org.Sc.sgd.db" \
  --pre_proc_pdf "../test-data/res/pre_proc.pdf" \
  --df_stats_out "../test-data/res/df_stats.tsv" \
  --outl_out     "../test-data/res/outl.tsv" \
  --data_wide_out "../test-data/res/data_wide.tsv" \
  --data_wide_symb_out "../test-data/res/data_wide_symb.tsv" \
  --exp_anal_pdf "../test-data/res/exp_anal.pdf" \
  --gene_net_pdf "../test-data/res/gene_net.pdf" \
  --imbe_out "../test-data/res/impact_betweenness.tsv" \
  --imbe_tab_out "../test-data/res/impact_betweenness_tab.tsv" \
  --kegg_en_out "../test-data/res/kegg_en.tsv" \
  --go_en_out "../test-data/res/go_en.tsv" 
