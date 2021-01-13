filter_thresh_test <- function(
  count_data, 
  DE_data, 
  threshs = seq(0, 5, 0.1),
  min_no_passed_thresh = 2,
  plot_dir
) {

  library(ggplot2)
  library(cowplot)

  threshs <- as.list(threshs)

  no_DE <- lapply(threshs, function(x) {

    # filter genes below threshold:
    keep <- rowSums(cpm(count_data)>x) >= 2
    keep_genes <- rownames(count_data)[keep]

    print(
      paste0(
        "No genes postfiltering with ", 
        x, 
        " CPM thresh = ", 
        length(keep_genes)
      )
    )

    # isolate DE values for kept genes:
    keep_DE <- DE_data[rownames(DE_data) %in% keep_genes,]

    # recalculate FDR and count significant genes:
    keep_DE$FDR <- p.adjust(keep_DE$PValue, method = "fdr")

    no_sig <- nrow(keep_DE[keep_DE$FDR < 0.05,])

    print(
      paste0(
        "No DE genes postfiltering with ", x, 
        " CPM thresh = ", no_sig
      )
    )

    res_df <- data.frame(
      CPM_thresh = x,
      no_sig = no_sig
    )

  })

  # bind results:
  thresh_vs_DE <- do.call("rbind", no_DE)

  # plot results:
  p <- ggplot(thresh_vs_DE, aes(x=CPM_thresh, y=no_sig))
  p <- p + geom_point()
  p <- p + theme_cowplot(12)
  p <- p + labs(x="CPM filter threshold", y="Number sig. genes (FDR<0.05)")

  png(
    paste0(
      plot_dir, 
      "CPM_threshold_vs_no_DE_genes_min_no_passed_thresh_", 
      min_no_passed_thresh, ".png"
    )
  )
    print(p)
  dev.off()

  pdf(
    paste0(
      plot_dir, 
      "CPM_threshold_vs_no_DE_genes_min_no_passed_thresh_", 
      min_no_passed_thresh, ".pdf"
    )
  )
    print(p)
  dev.off()

}