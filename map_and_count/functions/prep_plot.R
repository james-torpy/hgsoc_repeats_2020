prep_plot <- function(
  DE_data,
  FDR_lim = 0.05,
  FC_lim = 0.7,
  num_label
) {

  # extract results and convert gene symbols to names:
  DE_data$symbol <- rownames(DE_data)
  
  # colour genes with FDR < FDR_lim and -FC_lim =< logFC => FC_lim and check:
  DE_data$sig <- FALSE
  DE_data$sig[
    DE_data$FDR < FDR_lim & (
      DE_data$logFC >= FC_lim | DE_data$logFC <= -FC_lim
    )
  ] <- TRUE
  
  # make factor for colour scheme:
  DE_data$sig <- factor(DE_data$sig)
  
  # sort and identify top up and down significant genes:
  DE_data <- label_top(DE_data, num = num_label)
  
  return(DE_data)

}