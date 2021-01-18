return_top_DE <- function(
  in_dir,
  up_filename,
  down_filename
) {

  # load DE retrotransposons:
  up_DE <- read.table(
    paste0(in_dir, up_filename),
    header = T,
    sep = "\t",
    stringsAsFactors = F
  )
  sig_up <- up_DE[up_DE$FDR < 0.05,]
  sig_up <- sig_up[order(sig_up$logFC, decreasing = T),]
  
  down_DE <- read.table(
    paste0(in_dir, down_filename),
    header = T,
    sep = "\t",
    stringsAsFactors = F
  )
  sig_down <- down_DE[down_DE$FDR < 0.05,]
  sig_down <- sig_down[order(sig_down$logFC),]
  
  # take the top 10 up and 10 down retrotransposons:
  top_up <- sig_up[1:10,]
  top_down <- sig_down[1:10,]

  return(
    c(rownames(top_up), rownames(top_down))
  )

}