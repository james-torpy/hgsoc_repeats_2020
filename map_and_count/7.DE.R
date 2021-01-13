
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

  out_path <- paste0(
    results_dir, "DE/sub/"
  )

} else {

  out_path <- paste0(
    results_dir, "DE/"
  )

}

if (!incl_ascites) {
  out_path <- paste0(out_dir, "without_ascites/")
}

Robject_dir1 <- paste0(out_path, "Rdata/")
system(paste0("mkdir -p ", Robject_dir1))


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

# function to normalise and generate GLM using edgeR:
fit_glm <- dget(paste0(func_dir, "fit_glm.R"))

# function to perform DE using EdgeR:
do_DE <- dget(paste0(func_dir, "do_DE.R"))

# function to filter out genes using different CPM cutoffs and plot
# sig genes vs number of DE genes:
filter_thresh_test <- dget(paste0(func_dir, "filter_thresh_test.R"))

# function to compare DEs:
DE_compare <- dget(paste0(func_dir, "DE_compare.R"))

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
  all_counts <- all_counts[
    ,c(1:2, 41:42, 64:65, 80:81, 83:84, 85:86, 100:101, 105:106)
  ]
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

if (!file.exists(paste0(Robject_dir1, "post_aggregate_GC_counts.Rdata"))) {

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
    paste0(Robject_dir1, "post_aggregate_GC_counts.Rdata")
  )

} else {

  GC_counts <- readRDS(
    paste0(Robject_dir1, "post_aggregate_GC_counts.Rdata")
  )

}

# append GC and repeat counts:
formatted_counts <- rbind(GC_counts, repeat_counts)


#############################################################################
### 2. Use edgeR on raw counts to choose optimum filtering for 
# primary tumour vs FT ###
#############################################################################

Robject_dir2 <- paste0(out_path, "site/Rdata/")
system(paste0("mkdir -p ", Robject_dir2))
table_dir2 <- paste0(out_path, "site/tables/")
system(paste0("mkdir -p ", table_dir2))
plot_dir2 <- paste0(out_path, "site/plots/")
system(paste0("mkdir -p ", plot_dir2))

# normalise and fit data to GLM:
site_fit <- fit_glm(
  count_df = formatted_counts,
  sample_annot = sample_annot,
  repeat_symbols = repeat_symbols,
  cols = col_pal,
  div_type = "site",
  Robject_dir2,
  plot_dir2,
  func_dir
)

Robject_dir3 <- paste0(out_path, "site/primary_vs_FT/Rdata/")
system(paste0("mkdir -p ", Robject_dir3))
table_dir3 <- paste0(out_path, "site/primary_vs_FT/tables/")
system(paste0("mkdir -p ", table_dir3))
plot_dir3 <- paste0(out_path, "site/primary_vs_FT/plots/")
system(paste0("mkdir -p ", plot_dir3))

# perform primary vs FT DE:
primary_vs_FT_DE <- do_DE(
  fit_obj = site_fit$fit,
  edgeR_obj = site_fit$all_edgeR,
  con = c(0, -1, 0, 1),
  descrip = "primary_vs_FT",
  plot_dir = plot_dir3
)

filter_thresh_file <- list.files(
  plot_dir3, pattern = "CPM_threshold_vs_no_DE_genes"
)

if (length(filter_thresh_file) < 1) {
  filter_thresh_test(
    formatted_counts, 
    primary_vs_FT_DE, 
    threshs = seq(0, 5, 0.1),
    min_no_passed_thresh = 2,
    plot_dir = plot_dir3
  )
}

# chose not to filter by CPM as this gives the most significant results


#############################################################################
### 3. Plot primary tumour vs FT DE ###
#############################################################################

# load repeat annotation:
if (!file.exists(paste0(Robject_dir, "all_repeat_genes.Rdata"))) {

  repeat_gtf <- read.table(
    paste0(genome_dir, "custom3.repeats.hg38.gtf"),
    sep = "\t",
    header = F,
    fill = T
  )
  saveRDS(repeat_gtf, paste0(Robject_dir1, "all_repeat_genes.Rdata"))

  repeat_genes <- gsub("ID ", "", as.character(unique(repeat_gtf$V9)))

  repeat_info <- read.table(
  	paste0(genome_dir, "repeats.hg38.info.txt"),
  	sep = "\t",
  	header = T
  )
  saveRDS(repeat_info, paste0(Robject_dir1, "repeat_info.Rdata"))

} else {

  repeat_gtf <- readRDS(paste0(Robject_dir3, "all_repeat_genes.Rdata"))
  repeat_genes <- gsub("ID ", "", as.character(unique(repeat_gtf$V9)))
  repeat_info <- readRDS(paste0(Robject_dir1, "repeat_info.Rdata"))
  
}

# plot non-repeat DEs:
plot_DE(
  DE_results = primary_vs_FT_DE,
  DE_name = "primary_vs_FT",
  repeat_genes = repeat_genes,
  gene_type = "non_repeat",
  table_dir = table_dir3,
  plot_dir = plot_dir3,
  FDR_lim = 0.05,
  FC_lim = 0.7,
  num_label = 20,
  manual_lab = "none",
  dot_col = "#E6AB02",
  label_col = "#430F82",
  up_ctl = c(
  	"CDKN2A", "PTEN", "RAD51C", "PARP1", "E2F1",
  	"SBK1", "IGF1"
  ),
  up_ctl_col = "#E7298A"
)

# plot non-repeat DEs with ctls only labelled:
plot_DE(
  DE_results = primary_vs_FT_DE,
  DE_name = "primary_vs_FT_ctl_only_labelled",
  repeat_genes = repeat_genes,
  gene_type = "non_repeat",
  table_dir = table_dir3,
  plot_dir = plot_dir3,
  FDR_lim = 0.05,
  FC_lim = 0.7,
  num_label = 0,
  manual_lab = "none",
  dot_col = "#F2EF81",
  label_col = "#430F82",
  up_ctl = c(
  	"CDKN2A", "PTEN", "RAD51C", "PARP1", "E2F1",
  	"SBK1", "IGF1"
  ),
  up_ctl_col = "#E7298A"
)

# plot repeat DEs:
plot_DE(
  DE_results = primary_vs_FT_DE,
  DE_name = "primary_vs_FT",
  repeat_genes = repeat_genes,
  gene_type = "repeat",
  table_dir3,
  plot_dir3,
  FDR_lim = 0.05,
  FC_lim = 0.7,
  num_label = 20,
  manual_lab = "none",
  dot_col = "#F2EF81",
  label_col = "#430F82"
) 

# plot retrotransposons:
RT <- repeat_info$symbol[grep("LINE|SINE|LTR|SVA", repeat_info$type)]

plot_DE(
  DE_results = primary_vs_FT_DE,
  DE_name = "primary_vs_FT",
  repeat_genes = RT,
  gene_type = "retrotransposon",
  table_dir3,
  plot_dir3,
  FDR_lim = 0.05,
  FC_lim = 0.7,
  num_label = 10,
  manual_lab = "none",
  dot_col = "#F2EF81",
  label_col = "#BC8E1B",
  up_ctl = c(
  	"L1HS", "L1PA2", "AluYk2", "AluYd8", 
  	"AluYh3", "AluYb8"
  ),
  up_ctl_col = "#E7298A"
)


#############################################################################
### 4. Plot primary tumour vs FT FPKM of top 20 up and downreg 
# retrotransposons ###
#############################################################################

# calculate CPM:
all_counts <- subset(all_counts, select = -transcript_id)
reads_per_sample <- apply(all_counts, 2, sum)

repeat_CPM <- round(t(t(repeat_counts)/reads_per_sample)*1e6, 2)

# isolate top sig DE retrotransposons:
RT_DE <- primary_vs_FT_DE[RT,]



#############################################################################
### 5. Compare CCNE1 vs HRD DE to Nature paper results ###
#############################################################################

Robject_dir4 <- paste0(out_path, "GIN_driver/Rdata/")
system(paste0("mkdir -p ", Robject_dir4))
table_dir4 <- paste0(out_path, "GIN_driver/tables/")
system(paste0("mkdir -p ", table_dir4))
plot_dir4 <- paste0(out_path, "GIN_driver/plots/")
system(paste0("mkdir -p ", plot_dir4))

# normalise and fit data to GLM:
GIN_fit <- fit_glm(
  count_df = formatted_counts,
  sample_annot = sample_annot,
  repeat_symbols = repeat_symbols,
  cols = col_pal,
  div_type = "GIN_driver",
  Robject_dir4,
  plot_dir4,
  func_dir
)

Robject_dir5 <- paste0(out_path, "GIN_driver/CCNE_vs_HRD/Rdata/")
system(paste0("mkdir -p ", Robject_dir5))
table_dir5 <- paste0(out_path, "GIN_driver/CCNE_vs_HRD/tables/")
system(paste0("mkdir -p ", table_dir5))
plot_dir5 <- paste0(out_path, "GIN_driver/CCNE_vs_HRD/plots/")
system(paste0("mkdir -p ", plot_dir5))

# perform primary vs FT DE:
CCNE_vs_HRD_DE <- do_DE(
  fit_obj = GIN_fit$fit,
  edgeR_obj = GIN_fit$all_edgeR,
  con = c(1, 0, -1, 0, 0, 0),
  descrip = "CCNE_vs_HRD",
  plot_dir = plot_dir5
)

CCNE_vs_HRD_nat_vs_me <- DE_compare(
  DE_res = CCNE_vs_HRD_DE,
  descrip = "CCNE_vs_HRD",
  ref_dir,
  table_dir5,
  plot_dir5
)


#############################################################################
### 5. Compare CCNE1 vs HRD DE to Nature paper results ###
#############################################################################

Robject_dir6 <- paste0(out_path, "GIN_driver/HRD_vs_unknown/Rdata/")
system(paste0("mkdir -p ", Robject_dir6))
table_dir6 <- paste0(out_path, "GIN_driver/HRD_vs_unknown/tables/")
system(paste0("mkdir -p ", table_dir6))
plot_dir6 <- paste0(out_path, "GIN_driver/HRD_vs_unknown/plots/")
system(paste0("mkdir -p ", plot_dir6))

HRD_vs_unknown_DE <- do_DE(
  fit_obj = GIN_fit$fit,
  edgeR_obj = GIN_fit$all_edgeR,
  con = c(0, 0, 1, 0, -1, 0),
  descrip = "HRD_vs_unknown",
  plot_dir = plot_dir6
)

HRD_vs_unknown_nat_vs_me <- DE_compare(
  DE_res = HRD_vs_unknown_DE,
  descrip = "HRD_vs_unknown",
  ref_dir,
  table_dir6,
  plot_dir6
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


#############################################################################
### 6. Plot relapse tumour vs FT DE ###
#############################################################################

Robject_dir7 <- paste0(out_path, "site/primary_vs_ascites/Rdata/")
system(paste0("mkdir -p ", Robject_dir7))
table_dir7 <- paste0(out_path, "site/primary_vs_ascites/tables/")
system(paste0("mkdir -p ", table_dir7))
plot_dir7 <- paste0(out_path, "site/primary_vs_ascites/plots/")
system(paste0("mkdir -p ", plot_dir7))

primary_vs_ascites_DE <- do_DE(
  fit_obj = site_fit$fit,
  edgeR_obj = site_fit$all_edgeR,
  con = c(1, 0, 0, -1),
  descrip = "primary_vs_ascites",
  plot_dir = plot_dir7
)

# plot repeat DEs:
plot_DE(
  DE_results = primary_vs_ascites_DE,
  DE_name = "primary_vs_ascites",
  repeat_genes = repeat_genes,
  gene_type = "repeat",
  table_dir3,
  plot_dir3,
  FDR_lim = 0.05,
  FC_lim = 0.7,
  num_label = 10,
  manual_lab = "none",
  dot_col = "#F2EF81",
  label_col = "#430F82"
) 

# plot retrotransposons:
plot_DE(
  DE_results = primary_vs_ascites_DE,
  DE_name = "primary_vs_ascites",
  repeat_genes = RT,
  gene_type = "retrotransposon",
  table_dir3,
  plot_dir3,
  FDR_lim = 0.05,
  FC_lim = 0.7,
  num_label = 10,
  manual_lab = "none",
  dot_col = "#F2EF81",
  label_col = "#BC8E1B",
  up_ctl = "none",
  up_ctl_col = "#E7298A"
) 
