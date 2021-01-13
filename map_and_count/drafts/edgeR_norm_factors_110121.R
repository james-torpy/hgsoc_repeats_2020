edgeR_norm_factors <- function(
  DGE_object,
  prefix,
  div_type,
  dot_col,
  plot_dir,
  func_dir,
  mds_only = FALSE
) {

  library(ggplot2)
  library(cowplot)
  library(RUVSeq)

  plot_rle <- dget(paste0(func_dir, "plot_rle.R"))

  # create and check RLE plot pre-normalisation:
  print("Generating pre-norm RLE...")

  plot_rle(
    count_obj = DGE_object$counts,
    descrip = "pre_norm",
    plot_dir
  )

  # normalise counts and check:
  print("Normalising counts...")

### EdgeR ###

#  norm_y <- calcNormFactors(DGE_object)
#  norm_y$samples

  # multiply samples by norm factors:
  #norm_counts <- sweep(norm_y$counts, 2, norm_y$samples$norm.factors, "*")

#  print("Generating post-norm RLE...")

#  # plot RLE:
#  plot_rle(
#    count_obj = norm_counts,
#    descrip = "post_norm",
#    plot_dir
#  )

###

### RUV_Seq ###

  set <- newSeqExpressionSet(
    DGE_object$counts, 
    phenoData = data.frame(
      sample_annot$GIN_driver, 
      row.names=colnames(DGE_object$counts)
    )
  )

  nSet <- betweenLaneNormalization(set, which="full")

  print("Generating post-norm RLE...")

  # plot RLE:
  plot_rle(
    count_obj = normCounts(nSet),
    descrip = "post_norm",
    plot_dir
  )

  # add normalised counts to DGE_object
  norm_y <- DGE_object
  norm_y$counts <- normCounts(nSet)

###

  if (!file.exists(paste0(plot_dir, prefix, "_mds.png"))) {

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
    p <- p + theme_cowplot(12)
    p <- p + labs(
      x = "BCV distance 1", 
      y = "BCV distance 2", 
      colour = div_type
    )
    p <- p + scale_colour_manual(
      labels = levels(plot_df$group), values = dot_col
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

  }


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