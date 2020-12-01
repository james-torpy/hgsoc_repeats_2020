edgeR_norm_factors <- function(
  DGE_object,
  prefix,
  plot_dir
) {

  # normalise counts and check:
  norm_y <- calcNormFactors(DGE_object)
  norm_y$samples
  
  # visualise on MDS:
  png(paste0(plot_dir, prefix, "_mds.png"))
    plotMDS(
      norm_y,
      labels = norm_y$samples$group, 
      method = "bcv", 
      col = as.numeric(norm_y$samples$group)
    )
  dev.off()
  
  # estimate dispersion:
  disp_y <- estimateCommonDisp(norm_y)
  disp_y <- estimateTagwiseDisp(disp_y)
  
  # plot tagwise bcv:
  png(paste0(plot_dir, prefix, "_bcv_vs_cpm.png"))
    plotBCV(disp_y)
  dev.off()
  
  # estimate and compare dispersion by GLM:
  design <- model.matrix(~0 + norm_y$samples$group)
  colnames(design) <- levels(norm_y$samples$group)
  glm_y <- estimateGLMCommonDisp(norm_y, design)
  glm_y <- estimateGLMTrendedDisp(glm_y, design)
  glm_y <- estimateGLMTagwiseDisp(glm_y, design)
  
  # plot tagwise bcv:
  png(paste0(plot_dir, prefix, "_bcv_vs_cpm_glm.png"))
    plotBCV(glm_y)
  dev.off()

  return(list(design = design, non_GLM = disp_y, GLM = glm_y))

}