DE_compare <- function(
  DE_res,
  descrip,
  ref_dir,
  table_dir,
  plot_dir
) {

  library(ggplot2)

  # fetch Nature DE results:
  nat_DE <- read.table(
    paste0(ref_dir, "nat_", descrip, ".txt"),
    sep = "\t",
    header = T,
    stringsAsFactors = F
  )
  rownames(nat_DE) <- nat_DE$Gene
  nat_DE <- subset(nat_DE, select=-Gene)
  colnames(nat_DE) <- c(
    paste0(
      "nat_", colnames(nat_DE)
    )
  )
  
  # add logFC and p-vals from my DE data:
  both_DE <- merge(nat_DE, DE_res, by = 0)
  colnames(both_DE)[5:9] <- paste0(
    "my_", colnames(both_DE)[5:9]
  )
  
  # label significant genes:
  both_DE$sig <- FALSE
  both_DE$sig[both_DE$my_FDR < 0.05] <- TRUE
  both_DE$sig <- factor(
    both_DE$sig,
    levels = c(TRUE, FALSE)
  )
  
  # report significant genes and save table:
  print(
    paste0(
      "No. significant reported DE genes is ",
      length(both_DE$my_FDR[both_DE$my_FDR < 0.05]),
      " out of ",
      nrow(both_DE)
    )
  )
  write.table(
    subset(both_DE, select = -c(my_logCPM, my_LR)),
    paste0(table_dir, descrip, "_nature_vs_my_DE.txt"),
      sep = "\t",
      quote = F,
      row.names = F,
      col.names = T
  )
  
  # calculate correlation between my results and Nature:
  bm_cor <- cor.test(
    both_DE$my_logFC,
    both_DE$nat_logFC,
    method = "pearson"
  )
  cor_val <- round(bm_cor$estimate, 2)
  cor_p <- bm_cor$p.value
  
  # set up labels/annotation positions:
  if (cor_p < 0.001) {
    p_lab <- "< 0.001"
  } else {
    p_lab <- paste0("= ", as.character(round(cor_p, 8)))
  }
  x_annot_pos <- max(both_DE$my_logFC) + 
    (max(both_DE$nat_logFC)/4)
  y_annot_pos <- max(both_DE$my_logFC) + 
    (max(both_DE$my_logFC)/4)
  
  # plot nature vs my DE on density scatter:
  p <- ggplot(both_DE, aes(x=nat_logFC, y=my_logFC))
  p <- p + geom_point()
  p <- p + ylim(c(-5, 5))
  p <- p + xlim(c(-5, 5))
  p <- p + theme_cowplot(12)
  p <- p + labs(
    x="Patch et. al. logFC", 
    y="My DE method logFC"
  )
  p <- p + theme(
    text = element_text(size = 12),
    legend.title = element_blank()
  )
  p <- p + annotate(
    geom="text", 
    x=x_annot_pos, 
    y=x_annot_pos, 
    label=paste0("R2 = ", cor_val)
  )
  p <- p + annotate(
    geom="text", 
    x=x_annot_pos, 
    y=x_annot_pos-0.3, 
    label=paste0("p ", p_lab)
  )
  
  png(
    paste0(plot_dir, descrip, "_nature_vs_my_DE.png"),
    width = 10,
    height = 7,
    unit = "in",
    res = 300
  )
    print(p)
  dev.off()

  return(both_DE)

}