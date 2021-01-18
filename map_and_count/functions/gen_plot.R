gen_plot <- function(
  plot_df,
  dot_col,
  label_col,
  up_ctl_col,
  x_limits1 = "none",
  x_limits2 = "none"
) {

  library(ggplot2)
  library(cowplot)
  library(plyr)

  xlimit <- ceiling(max(abs(range(plot_df$logFC))))
  xlims <- c(-xlimit, xlimit)

  if (any(plot_df$up_ctl)) {
    up_ctl_df <- plot_df[plot_df$up_ctl,]
  }
  
  p <- ggplot(plot_df, aes(x=logFC, y=-log10(FDR), color=sig))
  p <- p + geom_point(
    plot_df, 
    mapping = aes(x=logFC, y=-log10(FDR), color=sig)
  )

  if (any(plot_df$up_ctl)) {

    p <- p + geom_point(
      up_ctl_df, 
      mapping = aes(x=logFC, y=-log10(FDR), color=sig)
    )
    p <- p + scale_color_manual(values = c(up_ctl_col, "grey", dot_col))

  } else if ("act" %in% plot_df$sig | "sup" %in% plot_df$sig) {

    p <- p + scale_color_manual(values = c("grey", "#BF3667", "#58B9DB"))

  } else {

    p <- p + scale_color_manual(values = c("grey", dot_col))

  }

  if ("act" %in% plot_df$sig | "sup" %in% plot_df$sig) {

    act_df <- plot_df[plot_df$sig == "act",]
    sup_df <- plot_df[plot_df$sig == "sup",]

    p <- p + geom_text_repel(
      act_df[act_df$label,], 
      mapping = aes(label=symbol),
      colour = "#BF3667",
      max.overlaps = 100
    )

     p <- p + geom_text_repel(
      sup_df[sup_df$label,], 
      mapping = aes(label=symbol),
      colour = "#58B9DB",
      max.overlaps = 100
    )

  } else {

    if (x_limits1[1] == "none") {
      p <- p + geom_text_repel(
        plot_df[plot_df$label,], 
        mapping = aes(label=symbol),
        colour = c(label_col),
        max.overlaps = 100
      )
    } else {
      p <- p + geom_text_repel(
        plot_df[plot_df$label,], 
        mapping = aes(label=symbol),
        colour = c(label_col),
        max.overlaps = 100,
        xlim = x_limits1
      )
    }
  
    if (any(plot_df$up_ctl)) {
      if (x_limits2[1] == "none") { 
        p <- p + geom_text_repel(
          plot_df[plot_df$up_ctl,], 
          mapping = aes(label=symbol),
          colour = c(up_ctl_col),
          max.overlaps = 100,
          xlim  = x_limits1
        )
      } else {
        p <- p + geom_text_repel(
          plot_df[plot_df$up_ctl,], 
          mapping = aes(label=symbol),
          colour = c(up_ctl_col),
          max.overlaps = 100,
          xlim  = x_limits2
        )
      }
    }

  }
  
  p <- p + theme_cowplot(12)
  p <- p + theme(legend.position = "none")
  p <- p + labs(x=paste0("log2 fold change"), y="-log10 FDR")
  p <- p + xlim(xlims)

  return(p)

}