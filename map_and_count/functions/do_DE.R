do_DE <- function(
  fit_obj,
  edgeR_obj,
  con,
  descrip,
  plot_dir
) {

  # perform DE on unfiltered counts, between CCNE and HRD:
  DE <- glmLRT(
    fit_obj, 
    contrast=con
  )
  
  # check DE numbers:
  sig <- decideTestsDGE(
    DE, adjust.method = "BH", p.value = 0.05
  )
  print(summary(sig))
  
  # generate smearplot as sanity check :
  tags <- rownames(edgeR_obj$GLM)[as.logical(sig)]
  png(paste0(plot_dir, descrip, "_smearplot.png"))
    plotSmear(DE, de.tags=tags)
    abline(h = c(-1, 1), col = "blue")
  dev.off()
  
  # fetch all DE and determine best CPM threshold (one with most sig genes)
  return(as.data.frame(topTags(DE, n=Inf)))

}