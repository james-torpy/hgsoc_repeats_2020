gen_plot <- function(
  plot_df,
  dot_col,
  label_col
) {

  xlimit <- ceiling(max(abs(range(plot_df$logFC))))
  xlims <- c(-xlimit, xlimit)
  
  p <- ggplot(plot_df, aes(x=logFC, y=-log10(FDR), color=sig))
  p <- p + geom_point(
    plot_df, 
    mapping = aes(x=logFC, y=-log10(FDR), color=sig)
  )
  p <- p + scale_color_manual(values = c("grey", dot_col))
  p <- p + geom_text_repel(
    plot_df[plot_df$label,], 
    mapping = aes(label=symbol),
    colour = label_col
  )
  p <- p + theme_cowplot(12)
  p <- p + theme(legend.position = "none")
  p <- p + labs(x=paste0("log2 fold change"), y="-log10 FDR")
  p <- p + xlim(xlims)

  return(p)

}