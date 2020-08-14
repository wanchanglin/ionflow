# wl-09-08-2020 \ Sun: Rscript test code for Linux. 

Rscript --vanilla ../ionflow.R \
  --ion_file "../test-data/iondata_test.tsv" \
  --std_file_sel "yes" \
  --std_file "../test-data/user_std.tsv" \
  --thres_clus "10.0" \
  --thres_anno "5.0" \
  --thres_corr "0.6" \
  --pre_proc_pdf "../test-data/res/pre_proc.pdf" \
  --df_stats_out "../test-data/res/df_stats.tsv" \
  --outl_out     "../test-data/res/outl.tsv" \
  --data_wide_out "../test-data/res/data_wide.tsv" \
  --data_wide_symb_out "../test-data/res/data_wide_symb.tsv" \
  --exp_anal_pdf "../test-data/res/exp_anal.pdf" \
  --gene_clus_pdf "../test-data/res/gene_clus.pdf" \
  --clus_out "../test-data/res/clus.tsv" \
  --anno_out "../test-data/res/kegg_go_anno.tsv" \
  --enri_out "../test-data/res/go_enri.tsv" \
  --gene_net_pdf "../test-data/res/gene_net.pdf" \
  --imbe_out "../test-data/res/impact_betweeness.tsv" \
  --imbe_tab_out "../test-data/res/impact_betweeness_tab.tsv"