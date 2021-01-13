plot_rle <- function(count_obj, descrip, plot_dir) {

  library(reshape)
  library(ggplot2)
  library(cowplot)
  # calculate log expression values:
  logs <- log(count_obj + 1e-6)
  # calculate median expression of each gene:
  meds <- apply(logs, 1, median)
  # calculate deviations from the median:
  devs <- logs - meds
  # prepare for boxplot:
  rle_df <- melt(devs)
  colnames(rle_df) <- c("gene", "sample", "deviation")
  # plot boxplot:
  p <- ggplot(rle_df, aes(x=sample, y=deviation))
  p <- p + geom_boxplot()
  p <- p + theme_cowplot(12)
  p <- p + labs(
    x="Sample", 
    y="Deviation"
  )
  p <- p + theme(
    axis.text.x = element_text(angle = 55, , hjust = 1)
  )

  pdf(paste0(plot_dir, descrip, "_rle.pdf"))
    print(p)
  dev.off()
  
  png(paste0(plot_dir, descrip, "_rle.png"))
    print(p)
  dev.off()
}