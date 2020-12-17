plot_DE <- function(
  DE_results,
  DE_name,
  repeat_genes,
  gene_type,
  table_dir,
  plot_dir,
  FDR_lim = 0.05,
  FC_lim = 0.7,
  num_label = 20,
  dot_col,
  label_col
) {

  print(paste0("Plotting ", gene_type, " DE for ", DE_name))

  # isolate genes:
  if (gene_type == "non_repeat") {
    DE <- DE_results[
      !(rownames(DE_results) %in% repeat_genes),
    ]
  } else  if (gene_type == "repeat") {
    DE <- DE_results[
      rownames(DE_results) %in% repeat_genes,
    ]
  }

  # prepare data:
  DE_labelled <- prep_plot(
    DE_data = DE,
    FDR_lim = FDR_lim,
    FC_lim = FC_lim,
    num_label = num_label
  )
  
  # save up and down genes as table:
  DE_up <- DE_labelled[DE_labelled$logFC > 0,]
  DE_up <- DE_up[order(DE_up$logFC),]
  head(DE_up)
  write.table(
    DE_up,
    paste0(table_dir, DE_name, "_upregulated_", gene_type, ".txt"),
    sep = "\t",
    col.names = T,
    row.names = T,
    quote = F
  )
  
  DE_down <- DE_labelled[DE_labelled$logFC < 0,]
  DE_down <- DE_down[order(DE_down$logFC),]
  write.table(
    DE_down,
    paste0(table_dir, DE_name, "_downregulated_", gene_type, ".txt"),
    sep = "\t",
    col.names = T,
    row.names = T,
    quote = F
  )
  
  # generate volcano plot:
  DE_plot <- gen_plot(
    plot_df = DE_labelled,
    dot_col,
    label_col
  )
  
  png(paste0(plot_dir, DE_name, "_", gene_type, "_DE_volcano.png"))
    print(DE_plot)
  dev.off()

}