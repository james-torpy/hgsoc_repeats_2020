edgeR_norm_factors <- function(
  DGE_object,
  prefix,
  div_type,
  dot_col,
  plot_dir,
  mds_only = FALSE
) {

  # normalise counts and check:
  print("Normalising counts...")

  norm_y <- calcNormFactors(DGE_object)
  norm_y$samples
  
  # prepare MDS data:
  print("Generating MDS...")

  mds <- plotMDS(
    norm_y,
    labels = norm_y$samples$group, 
    method = "bcv", 
    col = as.numeric(norm_y$samples$group)
  )
  dev.off()

  plot_df <- merge(
    norm_y$samples,
    as.data.frame(mds[[3]]),
    by=0
  )
  plot_df <- subset(plot_df, select = -c(lib.size, norm.factors))
  colnames(plot_df) <- c("sample", "group", "x", "y")

  # plot MDS:
  print("Plotting MDS...")

  p <- ggplot(plot_df, aes(x=x, y=y, colour=group))
  p <- p + geom_point()
  p <- p + labs(
    x = "BCV distance 1", 
    y = "BCV distance 2", 
    colour = div_type
  )
  p <- p + scale_colour_manual(
    labels = c("HRD", "CCNE amp.", "Unknown"), values = dot_col
  )
  p <- p + labs(colour = "GIN driver")

  pdf(
    paste0(plot_dir, prefix, "_mds.pdf"),
    width = 10
  )
    print(p)
  dev.off()

  png(
    paste0(plot_dir, prefix, "_mds.png"),
    height = 7,
    width = 10,
    units = "in",
    res = 300,
  )
    print(p)
  dev.off()

  if (!mds_only) {
  
    # estimate dispersion:
    print("Estimating dispersion...")

    disp_y <- estimateCommonDisp(norm_y)
    disp_y <- estimateTagwiseDisp(disp_y)
    
    # plot tagwise bcv:
    print("Plotting dispersion vs avg. log CPM...")

    png(paste0(plot_dir, prefix, "_bcv_vs_cpm.png"))
      plotBCV(disp_y)
    dev.off()
    
    # estimate and compare dispersion by GLM:
    print("Generating GLM...")

    design <- model.matrix(~0 + norm_y$samples$group)
    colnames(design) <- levels(norm_y$samples$group)
    glm_y <- estimateGLMCommonDisp(norm_y, design)
    glm_y <- estimateGLMTrendedDisp(glm_y, design)
    glm_y <- estimateGLMTagwiseDisp(glm_y, design)
    
    # plot tagwise bcv:
    print("Plotting GLM on dispersion vs avg. log CPM...")

    png(paste0(plot_dir, prefix, "_bcv_vs_cpm_glm.png"))
      plotBCV(glm_y)
    dev.off()

    return(list(design = design, non_GLM = disp_y, GLM = glm_y))

  } else {
    return(NULL)
  }


}