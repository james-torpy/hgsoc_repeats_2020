label_top <- function(DE_res, num, thresh = 0.05) {
 
  DE_res$label <- FALSE

  DE_res_sig <- DE_res[DE_res$sig == "sig",]
  
  up_genes <- head(
    DE_res_sig$symbol[with(DE_res_sig, order(-logFC, FDR))],
    num,
  )
  
  down_genes <- head(
    DE_res_sig$symbol[with(DE_res_sig, order(logFC, FDR))],
    num,
  )
  
  # label top genes and check:
  if (num == "sig") {

    DE_res$label[DE_res$FDR < thresh] <- TRUE

  } else {

    DE_res$label[
      DE_res$symbol %in% up_genes | 
      DE_res$symbol %in% down_genes
    ] <- TRUE

  }

  
  return(DE_res)

}