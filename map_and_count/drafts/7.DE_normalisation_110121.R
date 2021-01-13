
project_name <- "hgsoc_repeats/RNA-seq-final"

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/")
results_dir <- paste0(project_dir, "results/")
in_dir <- paste0(results_dir, "htseq/Rdata/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/map_and_count/functions/")
genome_dir <- paste0(project_dir, "genome/")
col_dir <- paste0(home_dir, "R/colour_palettes/")

sub <- FALSE
incl_ascites <- TRUE

if (sub) {

  out_dir <- paste0(
    results_dir, "DE/DE_normalisation_sub/"
  )

} else {

  out_dir <- paste0(
    results_dir, "DE/DE_normalisation/"
  )

}

if (!incl_ascites) {
  out_dir <- paste0(out_dir, "without_ascites/")
}

table_dir <- paste0(out_dir, "tables/")
system(paste0("mkdir -p ", table_dir))
plot_dir <- paste0(out_dir, "plots/")
system(paste0("mkdir -p ", plot_dir))
Robject_dir <- paste0(out_dir, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))


########################################################################
### 0.  Load packages and functions ###
########################################################################

library(edgeR)
library(DESeq2)
library(tibble)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggpointdensity)
library(ggrepel)

# function to calculate normalisation factors by edgeR:
edgeR_norm_factors <- dget(paste0(func_dir, "edgeR_norm_factors.R"))

# function to filter out genes using different CPM cutoffs and plot
# sig genes vs number of DE genes:
filter_thresh_test <- dget(paste0(func_dir, "filter_thresh_test.R"))

# function to label top 10 up and down genes:
label_top <- dget(paste0(func_dir, "label_top.R"))

# function to prepare data for plotting:
prep_plot <- dget(paste0(func_dir, "prep_plot.R"))

# volcano plot functions:
gen_plot <- dget(paste0(func_dir, "gen_plot.R"))
plot_DE <- dget(paste0(func_dir, "plot_DE.R"))

# load colours:
col_pal <- read.table(
  paste0(col_dir, "colour_palette_1.txt"),
  comment.char = "",
  stringsAsFactors = F
)[,1]


########################################################################
### 1.  Load counts and annotations ###
########################################################################

# fetch aggregated counts:
all_counts <- readRDS(
  paste0(in_dir, "all_combined_counts_ensembl_ids.Rdata")
)

if (sub) {
  all_counts <- all_counts[,1:10]
}

if (!incl_ascites) {
  all_counts <- all_counts[
    ,grep("-4|-8|-10", colnames(all_counts), invert = T)
  ]
}

# fetch sample annot:
sample_annot <- read.table(
  paste0(ref_dir, "sample_annot.txt"),
  sep = "\t",
  header = T,
  stringsAsFactors = F
)

# keep counts samples only:
sample_annot <- sample_annot[sample_annot$ID %in% colnames(all_counts),]

# only keep annotated samples:
all_counts <- all_counts[,colnames(all_counts) %in% sample_annot$ID]

# order according to all_counts:
all_counts <- all_counts[, match(sample_annot$ID, colnames(all_counts))]

# replace ensembl gene names with gene symbols:
# load exon info:
gencode_info <- read.table(
  paste0(genome_dir, "gencode.v35.basic.exon.info.txt"),
  sep = "\t",
  header = T,
  stringsAsFactors = F
)
gencode_info <- subset(gencode_info, select = -transcript_type)

# load repeat names:
repeat_symbols <- read.table(
  paste0(genome_dir, "repeats.hg38.symbols.txt"),
  header = F,
  stringsAsFactors = F
)[,1]

# keep only counts of transcripts within gencode info df, or repeats:
all_counts$transcript_id <- rownames(all_counts)
GC_counts <- all_counts[
  all_counts$transcript_id %in% gencode_info$transcript_id,
]

repeat_counts <- subset(
  all_counts[
    all_counts$transcript_id %in% repeat_symbols,
  ],
  select = -transcript_id
)

# append symbols onto counts:
GC_counts <- merge(GC_counts, gencode_info, by="transcript_id")
GC_counts <- subset(GC_counts, select = -transcript_id)

if (!file.exists(paste0(Robject_dir, "post_aggregate_GC_counts.Rdata"))) {

  # aggregate multiple counts of same symbol:
  print(
    paste0(
      "No transcripts before aggregation of genes with ",
      "multiple transcripts = ", nrow(GC_counts)
    )
  )

  GC_counts <- aggregate(.~symbol, GC_counts, sum)
  
  print(
    paste0(
      "No genes after aggregation of genes with ",
      "multiple transcripts = ", nrow(GC_counts)
    )
  )

  GC_counts <- GC_counts %>%
    column_to_rownames("symbol")

  saveRDS(
    GC_counts, 
    paste0(Robject_dir, "post_aggregate_GC_counts.Rdata")
  )

} else {

  GC_counts <- readRDS(
    paste0(Robject_dir, "post_aggregate_GC_counts.Rdata")
  )

}

# append GC and repeat counts:
formatted_counts <- rbind(GC_counts, repeat_counts)


#############################################################################
### 2. Use edgeR on raw counts to choose optimum filtering  ###
#############################################################################

fit_model <- function(
  count_df,
  sample_annot,
  repeat_symbols,
  cols,
  div_type,
  Robject_dir,
  plot_dir
) {

  # create DGEList objects for edgeR:
  DGE_list <- list(
    all = DGEList(
      counts=count_df, 
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
        edgeR_list <- list(edgeR_y)
      } else {
        edgeR_list[[i]] <- edgeR_y
      }
  
    }
    names(edgeR_list) <- names(DGE_list)
  
    all_edgeR <- edgeR_list$all
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

fit_data <- fit_model(
  count_df = formatted_counts,
  sample_annot = sample_annot,
  repeat_symbols = repeat_symbols,
  cols = col_pal,
  div_type = "GIN_driver",
  Robject_dir,
  plot_dir
)

do_DE <- function(
  fit_obj,
  edgeR_obj,
  descrip,
  plot_dir
)
  # perform DE on unfiltered counts, between CCNE and HRD:
  DE <- glmLRT(
    fit, 
    contrast=c(1, 0, -1, 0)
  )
  
  # check DE numbers:
  sig <- decideTestsDGE(
    DE, adjust.method = "BH", p.value = 0.05
  )
  print(summary(sig))
  
  # generate smearplot as sanity check :
  tags <- rownames(all_edgeR$GLM)[as.logical(sig)]
  png(paste0(plot_dir, descrip, "_smearplot.png"))
    plotSmear(DE, de.tags=tags)
    abline(h = c(-1, 1), col = "blue")
  dev.off()
  
  # fetch all DE and determine best CPM threshold (one with most sig genes)
  return(as.data.frame(topTags(DE, n=Inf)))

}

CCNE_vs_HRD_DE <- do_DE(
  fit_obj = CCNE_vs_HRD_fit$fit,
  edgeR_obj = CCNE_vs_HRD_fit$all_edgeR,
  descrip = "CCNE_vs_HRD",
  plot_dir
)

filter_thresh_file <- list.files(
  plot_dir, pattern = "CPM_threshold_vs_no_DE_genes"
)

if (length(filter_thresh_file) < 1) {
  filter_thresh_test(
    formatted_counts, 
    CCNE_vs_HRD, 
    threshs = seq(0, 5, 0.1),
    min_no_passed_thresh = 2
  )
}


# chose not to filter by CPM as this gives the most significant results


#############################################################################
### 2. Compare CCNE1 vs HRD DE to Nature paper results ###
#############################################################################

DE_compare <- function(
  DE_res,
  descrip,
  ref_dir,
  table_dir,
  plot_dir
) {

  # fetch Nature DE results:
  nat_DE <- read.table(
    paste0(ref_dir, "nat_", descrip, ".txt"),
    header = T,
    stringsAsFactors = F
  )
  rownames(nat_DE) <- nat_DE$Gene
  nat_DE <- subset(nat_DE, select=-Gene)
  colnames(nat_DE) <- c(
    paste0(
      "nat_", colnames(nat_DE)
    )
  )
  
  # add logFC and p-vals from my DE data:
  both_DE <- merge(nat_DE, DE_res, by = 0)
  colnames(both_DE)[5:9] <- paste0(
    "my_", colnames(both_DE)[5:9]
  )
  
  # label significant genes:
  both_DE$sig <- FALSE
  both_DE$sig[both_DE$my_FDR < 0.05] <- TRUE
  both_DE$sig <- factor(
    both_DE$sig,
    levels = c(TRUE, FALSE)
  )
  
  # report significant genes and save table:
  print(
    paste0(
      "No. significant reported DE genes is ",
      length(both_DE$my_FDR[both_DE$my_FDR < 0.05]),
      " out of ",
      nrow(both_DE)
    )
  )
  write.table(
    subset(both_DE, select = -c(my_logCPM, my_LR)),
    paste0(table_dir, descrip, "_nature_vs_my_DE.txt"),
      sep = "\t",
      quote = F,
      row.names = F,
      col.names = T
  )
  
  # calculate correlation between my results and Nature:
  bm_cor <- cor.test(
    both_DE$my_logFC,
    both_DE$nat_logFC,
    method = "pearson"
  )
  cor_val <- round(bm_cor$estimate, 2)
  cor_p <- bm_cor$p.value
  
  # set up labels/annotation positions:
  if (cor_p < 0.001) {
    p_lab <- "< 0.001"
  } else {
    p_lab <- paste0("= ", as.character(round(cor_p, 8)))
  }
  x_annot_pos <- max(both_DE$my_logFC) + 
    (max(both_DE$nat_logFC)/4)
  y_annot_pos <- max(both_DE$my_logFC) + 
    (max(both_DE$my_logFC)/4)
  
  # plot nature vs my DE on density scatter:
  p <- ggplot(both_DE, aes(x=nat_logFC, y=my_logFC, colour = sig))
  p <- p + geom_point()
  p <- p + ylim(c(-5, 5))
  p <- p + xlim(c(-5, 5))
  p <- p + scale_color_manual(
    values = c("#1B9E77", "grey"),
    labels = c("sig.", "non-sig.")
  )
  p <- p + theme_cowplot(12)
  p <- p + labs(
    x="Patch et. al. logFC", 
    y="My DE method logFC"
  )
  p <- p + theme(
    text = element_text(size = 12),
    legend.title = element_blank()
  )
  p <- p + annotate(
    geom="text", 
    x=x_annot_pos, 
    y=x_annot_pos, 
    label=paste0("R2 = ", cor_val)
  )
  p <- p + annotate(
    geom="text", 
    x=x_annot_pos, 
    y=x_annot_pos-0.3, 
    label=paste0("p ", p_lab)
  )
  
  png(
    paste0(plot_dir, descrip, "_nature_vs_my_DE.png"),
    width = 10,
    height = 7,
    unit = "in",
    res = 300
  )
    print(p)
  dev.off()

  return(both_DE)

}

CCNE_vs_HRD_nat_vs_me <- DE_compare(
  DE_res = CCNE_vs_HRD,
  descrip = "CCNE_vs_HRD",
  ref_dir,
  table_dir,
  plot_dir
)


#############################################################################
### 3. Compare CCNE1 vs HRD DE to Nature paper results ###
#############################################################################

HRD_vs_Unknown <- do_DE(
  count_df = formatted_counts,
  sample_annot = sample_annot,
  repeat_symbols = repeat_symbols,
  descrip = "HRD_vs_Unknown",
  cols = col_pal,
  div_type = "GIN_driver",
  Robject_dir,
  plot_dir
)

HRD_vs_Unknown_nat_vs_me <- DE_compare(
  DE_res = HRD_vs_Unknown,
  descrip = "HRD_vs_Unknown",
  ref_dir,
  table_dir,
  plot_dir
)


## find genes with most discrepency:
#both_CCNE_vs_HRD$diff <- both_CCNE_vs_HRD$nat_logFC - 
#  both_CCNE_vs_HRD$my_logFC
#disc_df <- both_CCNE_vs_HRD[order(both_CCNE_vs_HRD$diff, decreasing=F),]
#
#top_diff <- head(disc_df$Row.names, 10)
#
## fetch counts for these genes:
#diff_counts <- formatted_counts[top_diff,]
#
## split by driver:
#diff_CCNE <- diff_counts[
#  ,colnames(diff_counts) %in% 
#    sample_annot$ID[sample_annot$GIN_driver == "CCNE"]
#]
#
#diff_HRD <- diff_counts[
#  ,colnames(diff_counts) %in% 
#    sample_annot$ID[sample_annot$GIN_driver == "HRD"]
#]

