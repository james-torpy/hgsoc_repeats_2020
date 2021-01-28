fit_glm <- function(
  count_df,
  sample_annot,
  repeat_symbols,
  cols,
  div_type,
  Robject_dir,
  plot_dir,
  func_dir
) {

  edgeR_norm_factors <- dget(paste0(func_dir, "edgeR_norm_factors.R"))

  # create DGEList objects for edgeR:
  DGE_list <- list(
    all = DGEList(
      counts=count_df, 
      group = eval(parse(text = paste0("sample_annot$", div_type)))
    ),
    gc = DGEList(
      counts = count_df[!(rownames(count_df) %in% repeat_symbols), ], 
      group = eval(parse(text = paste0("sample_annot$", div_type)))
    ),
    repeats = DGEList(
      counts = count_df[rownames(count_df) %in% repeat_symbols, ], 
      group = eval(parse(text = paste0("sample_annot$", div_type)))
    )
  )
  
  # check library size:
  DGE_list$all$samples$lib.size
  
  if (!file.exists(paste0(Robject_dir, "edgeR.Rdata"))) {
  
    # calculate normalisation factors and check QC plots of different methods:
    for (i in 1:length(DGE_list)) {
  
      print(
        paste0(
          "Calculating normalisation factors for and MDS plotting ", 
          names(DGE_list)[i]
        )
      )
  
  	  if (i==1) {
  	  	edgeR_y <- edgeR_norm_factors(
          DGE_object = DGE_list[[i]], 
          prefix = names(DGE_list)[i],
          div_type = div_type,
          dot_col = col_pal,
          plot_dir,
          func_dir
        )
  	  } else {
  	  	edgeR_y <- edgeR_norm_factors(
          DGE_object = DGE_list[[i]], 
          prefix = names(DGE_list)[i],
          div_type = div_type,
          dot_col = col_pal,
          plot_dir,
          func_dir,
          mds_only = TRUE
        )
  	  }
  
      if (i==1) {
        all_edgeR <- edgeR_y
      }
  
    }
  
    saveRDS(all_edgeR, paste0(Robject_dir, "all_edgeR.Rdata"))
  
  } else {
    all_edgeR <- readRDS(paste0(Robject_dir, "all_edgeR.Rdata"))
  }
  
  # fit data to model:
  fit <- glmFit(all_edgeR$GLM, all_edgeR$design)

  return(
  	list(
  	  fit = fit,
  	  all_edgeR = all_edgeR$GLM
  	)
  )

}