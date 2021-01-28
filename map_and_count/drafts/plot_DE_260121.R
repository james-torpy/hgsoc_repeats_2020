plot_DE <- function(
  DE_results,
  DE_name,
  repeat_genes,
  gene_type,
  table_dir,
  plot_dir,
  FDR_lim = 0.05,
  FC_lim = 0.7,
  num_label = 10,
  manual_lab = "none",
  dot_col,
  label_col,
  up_ctl = "none",
  up_ctl_col = "none",
  remove_label = NULL,
  act_sup_cols = FALSE,
  act_genes = "none",
  sup_genes = "none",
  plot_CPM = FALSE,
  CPM_data = NULL

) {

  print(paste0("Plotting ", gene_type, " DE for ", DE_name))

  # isolate genes:
  if (gene_type == "non_repeat") {
    DE <- DE_results[
      !(rownames(DE_results) %in% repeat_genes),
    ]
  } else {
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

  if (up_ctl[1] != "none") {
    # label upregulated controls:
    DE_labelled$up_ctl <- FALSE
    DE_labelled$up_ctl[rownames(DE_labelled) %in% up_ctl] <- TRUE
    DE_labelled$sig[rownames(DE_labelled) %in% up_ctl] <- "up_ctl"
    DE_labelled$label[rownames(DE_labelled) %in% up_ctl] <- FALSE
  }

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

  # remove any labels needed:
  DE_labelled$label[rownames(DE_labelled) %in% remove_label] <- FALSE

  if (plot_CPM) {

    DE_sig <- DE_labelled[DE_labelled$sig == "sig",]



  }

  # colour activator or supprssor genes:
  if (act_sup_cols) {

    DE_labelled$sig[
      rownames(DE_labelled) %in% act_genes & DE_labelled$sig == "sig"
    ] <- "act"
    DE_labelled$sig[
      rownames(DE_labelled) %in% sup_genes & DE_labelled$sig == "sig"
    ] <- "sup"
  
    # make sig column a factor and adjust levels:
    DE_labelled$sig <- factor(
      DE_labelled$sig, levels = c("act", "sup", "non_sig")
    )

  } else {

    # make sig column a factor and adjust levels:
    DE_labelled$sig <- factor(
      DE_labelled$sig, levels = c("up_ctl", "sig", "non_sig")
    )

  }
  
  # generate volcano plot:
  DE_plot <- gen_plot(
    plot_df = DE_labelled,
    dot_col,
    label_col,
    up_ctl_col,
    x_limits1 = "none",
    x_limits2 = "none"
  )

  
  png(
    paste0(plot_dir, DE_name, "_", gene_type, "_DE_volcano.png"),
    width = 9,
    height = 7,
    res = 300,
    units = "in"
  )
    print(DE_plot)
  dev.off()

#  pdf(paste0(plot_dir, DE_name, "_", gene_type, "_DE_volcano.pdf"))
#    print(DE_plot)
#  dev.off()

  ggsave(
    paste0(plot_dir, DE_name, "_", gene_type, "_DE_volcano.pdf"),
    width = 9,
    height = 7,
    useDingbats = FALSE
  )
    print(DE_plot)
  dev.off()

}