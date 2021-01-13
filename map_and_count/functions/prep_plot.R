prep_plot <- function(
  DE_data,
  FDR_lim = 0.05,
  FC_lim = 0.7,
  num_label
) {

  # extract results and convert gene symbols to names:
  DE_data$symbol <- rownames(DE_data)

  # remove novel AC and AL genes:
  DE_data <- DE_data[
    grep("A[C,L,P][0-9][0-9][0-9][0-9]", DE_data$symbol, invert = T),
  ]
  # remove immune genes:
  DE_data <- DE_data[
    grep("IG[H,K,L]", DE_data$symbol, invert = T),
  ]
  # remove LINCs:
  DE_data <- DE_data[
    grep("LINC", DE_data$symbol, invert = T),
  ]
  # remove MIRs:
  DE_data <- DE_data[
    grep("MIR", DE_data$symbol, invert = T),
  ]
  
  # colour genes with FDR < FDR_lim and -FC_lim =< logFC => FC_lim and check:
  DE_data$sig <- "non_sig"
  DE_data$sig[
    DE_data$FDR < FDR_lim & (
      DE_data$logFC >= FC_lim | DE_data$logFC <= -FC_lim
    )
  ] <- "sig"
  
  # make factor for colour scheme:
  DE_data$sig <- DE_data$sig
  
  # sort and identify top up and down significant genes:
  DE_data <- label_top(DE_data, num = num_label)
  
  return(DE_data)

}