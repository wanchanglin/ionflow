# wl-09-08-2020, Sun: Rscript test code for Linux. 
# wl-12-11-2020, Thu: Modification 

Rscript --vanilla ../ionflow.R \
  --ion_file "~/R_lwc/r_data/icl/test-data/ionome_oe_test.tsv" \
  --var_id "1" \
  --batch_id "4" \
  --data_id "5" \
  --method_norm "median" \
  --method_outliers "IQR" \
  --thres_outl 3 \
  --stand_method "std" \
  --thres_symb 2 \
  --min_clust_size 10.0 \
  --thres_corr "0.75" \
  --method_corr "pearson" \
  --pval 0.05 \
  --ont "BP" \
  --annot_pkg "org.Sc.sgd.db" \
  --pre_proc_pdf "../test-data/res/oe_pre_proc.pdf" \
  --df_stats_out "../test-data/res/oe_df_stats.tsv" \
  --outl_out     "../test-data/res/oe_outl.tsv" \
  --data_wide_out "../test-data/res/oe_data_wide.tsv" \
  --data_wide_symb_out "../test-data/res/oe_data_wide_symb.tsv" \
  --exp_anal_pdf "../test-data/res/oe_exp_anal.pdf" \
  --gene_net_pdf "../test-data/res/oe_gene_net.pdf" \
  --imbe_out "../test-data/res/oe_impact_betweenness.tsv" \
  --imbe_tab_out "../test-data/res/oe_impact_betweenness_tab.tsv" \
  --kegg_en_out "../test-data/res/oe_kegg_en.tsv" \
  --go_en_out "../test-data/res/oe_go_en.tsv" 
