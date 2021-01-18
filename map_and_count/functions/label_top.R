label_top <- function(DE_res, num, thresh = 0.05) {
 
  DE_res$label <- FALSE
  
  up_genes <- head(
    DE_res$symbol[with(DE_res, order(-logFC, FDR))],
    num,
  )
  
  down_genes <- head(
    DE_res$symbol[with(DE_res, order(logFC, FDR))],
    num,
  )
  
  # label top genes and check:
  if (num == "sig") {

    DE_res$label[DE_res$FDR < thresh] <- TRUE

  } else {

    DE_res$label[
      DE_res$symbol %in% up_genes & DE_res$sig == "sig" | 
      DE_res$symbol %in% down_genes & DE_res$sig == "sig"
    ] <- TRUE

  }

  
  return(DE_res)

}