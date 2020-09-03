# wl-09-08-2020 \ Sun: Rscript test code for Linux. 
# wl-24-08-2020, Mon: test full data with default setting

Rscript --vanilla ../ionflow.R \
  --ion_file "../test-data/iondata.tsv" \
  --std_file_sel "yes" \
  --std_file "../test-data/user_std.tsv" \
  --thres_clus "10.0" \
  --thres_anno "5.0" \
  --thres_corr "0.6" \
  --pre_proc_pdf "../test-data/res/pre_proc_1.pdf" \
  --df_stats_out "../test-data/res/df_stats_1.tsv" \
  --outl_out     "../test-data/res/outl_1.tsv" \
  --data_wide_out "../test-data/res/data_wide_1.tsv" \
  --data_wide_symb_out "../test-data/res/data_wide_symb_1.tsv" \
  --exp_anal_pdf "../test-data/res/exp_anal_1.pdf" \
  --gene_clus_pdf "../test-data/res/gene_clus_1.pdf" \
  --clus_out "../test-data/res/clus_1.tsv" \
  --anno_out "../test-data/res/kegg_go_anno_1.tsv" \
  --enri_out "../test-data/res/go_enri_1.tsv" \
  --gene_net_pdf "../test-data/res/gene_net_1.pdf" \
  --imbe_out "../test-data/res/impact_betweenness_1.tsv" \
  --imbe_tab_out "../test-data/res/impact_betweenness_tab_1.tsv"