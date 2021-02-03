
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
incl_primary_ascites <- FALSE
dot_col <- "#D3BD39"
label_col <- "#D95F02"
up_ctl_col <- "#E7298A"

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
  out_path <- paste0(out_path, "without_ascites/")
} else if (!incl_primary_ascites) {
  out_path <- paste0(out_path, "without_primary_ascites/")
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

# fetch sample annot:
sample_annot <- read.table(
  paste0(ref_dir, "sample_annot.txt"),
  sep = "\t",
  header = T,
  stringsAsFactors = F
)

# remove ascites samples if needed:
if (!incl_ascites) {
  all_counts <- all_counts[
    ,grep("-4|-8|-10", colnames(all_counts), invert = T)
  ]
} else if (!incl_primary_ascites) {
  # identify and remove primary ascites:
  primary_ascites <- sample_annot$ID[sample_annot$site == "primary_ascites"]
  all_counts <- all_counts[
    ,!(colnames(all_counts) %in% primary_ascites)
  ]
}

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

# load repeat annotation:
if (!file.exists(paste0(Robject_dir1, "all_repeat_genes.Rdata"))) {

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

  repeat_gtf <- readRDS(paste0(Robject_dir1, "all_repeat_genes.Rdata"))
  repeat_genes <- gsub("ID ", "", as.character(unique(repeat_gtf$V9)))
  repeat_info <- readRDS(paste0(Robject_dir1, "repeat_info.Rdata"))
  
}

# isolate retrotransposons:
RT <- repeat_info$symbol[grep("LINE|SINE|LTR|SVA", repeat_info$type)]

# isolate centromere RNA
sat <- repeat_info$symbol[grep("Satellite|centromere", repeat_info$type)]
sat <- sat[!(sat %in% c("MSR1", "centromere", "SATR2"))]


# check order of annotation and counts dfs:
print(
  paste0(
    "Sample annotation in the same order as counts file? ",
    identical(sample_annot$ID, colnames(formatted_counts))
  )
)


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
if (!file.exists(paste0(Robject_dir2, "site_fit.Rdata"))) {

  site_fit <- fit_glm(
    count_df = formatted_counts,
    sample_annot = sample_annot,
    repeat_symbols = repeat_symbols,
    cols = col_pal,
    div_type = "site",
    Robject_dir = Robject_dir2,
    plot_dir = plot_dir2,
    func_dir
  )

  saveRDS(site_fit, paste0(Robject_dir2, "site_fit.Rdata"))

} else {
  site_fit <- readRDS(paste0(Robject_dir2, "site_fit.Rdata"))
}

Robject_dir3 <- paste0(out_path, "site/primary_vs_FT/Rdata/")
system(paste0("mkdir -p ", Robject_dir3))
table_dir3 <- paste0(out_path, "site/primary_vs_FT/tables/")
system(paste0("mkdir -p ", table_dir3))
plot_dir3 <- paste0(out_path, "site/primary_vs_FT/plots/")
system(paste0("mkdir -p ", plot_dir3))

if (!incl_ascites) {
  primary_vs_FT_con <- c(-1, 1)
} else if (!incl_primary_ascites) {
  primary_vs_FT_con <- c(-1, 0, 1)
} else {
  primary_vs_FT_con <- c(-1, 0, 0, 1)
}

# perform primary vs FT DE:
primary_vs_FT_DE <- do_DE(
  fit_obj = site_fit$fit,
  edgeR_obj = site_fit$all_edgeR,
  con = primary_vs_FT_con,
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
  dot_col = dot_col,
  label_col = label_col,
  up_ctl = c(
  	"PARP1", "E2F1", "SBK1", "IGF1"
  ),
  up_ctl_col = up_ctl_col
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
  dot_col = dot_col,
  label_col = label_col,
  up_ctl = c(
  	"PARP1", "E2F1",
  	"SBK1", "IGF1"
  ),
  up_ctl_col = up_ctl_col
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
  dot_col = dot_col,
  label_col = label_col
) 

# calculate CPMs:
if (!exists("CPM")) {

  if (!file.exists(paste0(Robject_dir1, "CPM.Rdata"))) {

    # calculate CPM:
    all_counts <- rbind(GC_counts, repeat_counts)
    reads_per_sample <- apply(all_counts, 2, sum)
    
    CPM <- round(t(t(all_counts)/reads_per_sample)*1e6, 2)

    saveRDS(CPM, paste0(Robject_dir1, "CPM.Rdata"))

  } else {
    CPM <- readRDS(paste0(Robject_dir1, "CPM.Rdata"))
  }

}

# plot retrotransposons:
#RT <- RT[!(RT %in% c("MER57F", "HERVL74-int"))]
plot_DE(
  DE_results = primary_vs_FT_DE,
  DE_name = "primary_vs_FT",
  repeat_genes = RT,
  gene_type = "retrotransposon",
  table_dir = table_dir3,
  plot_dir = plot_dir3,
  FDR_lim = 0.05,
  FC_lim = 0.7,
  num_label = 10,
  manual_lab = "none",
  dot_col = dot_col,
  label_col = label_col,
  up_ctl = c(
  	"L1HS", "L1PA2", "AluYk2", "AluYd8", 
  	"AluYh3", "AluYb8"
  ),
  up_ctl_col = up_ctl_col
)

# plot satellite RNA:
plot_DE(
  DE_results = primary_vs_FT_DE,
  DE_name = "primary_vs_FT",
  repeat_genes = sat,
  gene_type = "satellite_fdr_0.1",
  table_dir3,
  plot_dir3,
  FDR_lim = 0.1,
  FC_lim = 0.7,
  num_label = 10,
  manual_lab = "none",
  dot_col = dot_col,
  label_col = label_col,
  up_ctl = c(
    "HSATII"
  ),
  up_ctl_col = up_ctl_col
)


#############################################################################
### 4. Plot relapse tumour vs primary DE ###
#############################################################################

Robject_dir4 <- paste0(out_path, "site/recurrent_vs_primary/Rdata/")
system(paste0("mkdir -p ", Robject_dir4))
table_dir4 <- paste0(out_path, "site/recurrent_vs_primary/tables/")
system(paste0("mkdir -p ", table_dir4))
plot_dir4 <- paste0(out_path, "site/recurrent_vs_primary/plots/")
system(paste0("mkdir -p ", plot_dir4))

if (!incl_primary_ascites) {
  recurrent_vs_primary_con <- c(0, 1, -1)
} else {
  recurrent_vs_primary_con <- c(0, 0, 1, -1)
}

recurrent_vs_primary_DE <- do_DE(
  fit_obj = site_fit$fit,
  edgeR_obj = site_fit$all_edgeR,
  con = recurrent_vs_primary_con,
  descrip = "recurrent_vs_primary",
  plot_dir = plot_dir4
)

# plot retrotransposon DEs:
plot_DE(
  DE_results = recurrent_vs_primary_DE,
  DE_name = "recurrent_vs_primary",
  repeat_genes = RT,
  gene_type = "retrotransposon",
  table_dir4,
  plot_dir4,
  FDR_lim = 0.05,
  FC_lim = 0.7,
  num_label = 10,
  manual_lab = "none",
  dot_col = dot_col,
  label_col = label_col,
  up_ctl = c("AluYh9", "AluYk2"),
  up_ctl_col = up_ctl_col
)

# plot satellite DEs:
plot_DE(
  DE_results = recurrent_vs_primary_DE,
  DE_name = "recurrent_vs_primary",
  repeat_genes = sat,
  gene_type = "satellite_fdr_0.1",
  table_dir4,
  plot_dir4,
  FDR_lim = 0.1,
  FC_lim = 0.7,
  num_label = 10,
  manual_lab = "none",
  dot_col = dot_col,
  label_col = label_col,
  up_ctl = c(
    "HSATII"
  ),
  up_ctl_col = up_ctl_col
)


#############################################################################
### 5. Plot relapse tumour vs FT DE ###
#############################################################################

Robject_dir5 <- paste0(out_path, "site/recurrent_vs_FT/Rdata/")
system(paste0("mkdir -p ", Robject_dir5))
table_dir5 <- paste0(out_path, "site/recurrent_vs_FT/tables/")
system(paste0("mkdir -p ", table_dir5))
plot_dir5 <- paste0(out_path, "site/recurrent_vs_FT/plots/")
system(paste0("mkdir -p ", plot_dir5))

if (!incl_primary_ascites) {
  recurrent_vs_FT_con <- c(-1, 1, 0)
} else {
  recurrent_vs_FT_con <- c(-1, 0, 1, 0)
}

recurrent_vs_FT_DE <- do_DE(
  fit_obj = site_fit$fit,
  edgeR_obj = site_fit$all_edgeR,
  con = recurrent_vs_FT_con,
  descrip = "recurrent_vs_FT",
  plot_dir = plot_dir5
)

# plot retrotransposon DEs:
plot_DE(
  DE_results = recurrent_vs_FT_DE,
  DE_name = "recurrent_vs_FT",
  repeat_genes = RT,
  gene_type = "retrotransposon",
  table_dir5,
  plot_dir5,
  FDR_lim = 0.05,
  FC_lim = 0.7,
  num_label = 10,
  manual_lab = "none",
  dot_col = dot_col,
  label_col = label_col,
  up_ctl = c("AluYh9", "AluYk2", "L1PA2", "L1HS"),
  up_ctl_col = up_ctl_col
)

# plot satellite DEs:
plot_DE(
  DE_results = recurrent_vs_FT_DE,
  DE_name = "recurrent_vs_FT",
  repeat_genes = sat,
  gene_type = "satellite_fdr_0.1",
  table_dir5,
  plot_dir5,
  FDR_lim = 0.1,
  FC_lim = 0.7,
  num_label = 10,
  manual_lab = "none",
  dot_col = dot_col,
  label_col = label_col,
  up_ctl = c(
    "HSATII"
  ),
  up_ctl_col = up_ctl_col
)


#############################################################################
### 5. Compare CCNE1 vs HRD DE to Nature paper results ###
#############################################################################

Robject_dir6 <- paste0(out_path, "GIN_driver/Rdata/")
system(paste0("mkdir -p ", Robject_dir6))
table_dir6 <- paste0(out_path, "GIN_driver/tables/")
system(paste0("mkdir -p ", table_dir6))
plot_dir6 <- paste0(out_path, "GIN_driver/plots/")
system(paste0("mkdir -p ", plot_dir6))

GIN_annot <- sample_annot

# remove ascites and fallopian tissue:
GIN_annot <- GIN_annot[grep("ascites|fallopian_tissue", GIN_annot$site, invert = T), ]
GIN_counts <- formatted_counts[,colnames(formatted_counts) %in% GIN_annot$ID]

# check order of annotation and counts dfs:
print(
  paste0(
    "Sample annotation in the same order as counts file? ",
    identical(GIN_annot$ID, colnames(GIN_counts))
  )
)

# normalise and fit data to GLM:
GIN_fit <- fit_glm(
  count_df = GIN_counts,
  sample_annot = GIN_annot,
  repeat_symbols = repeat_symbols,
  cols = col_pal,
  div_type = "GIN_driver",
  Robject_dir6,
  plot_dir6,
  func_dir
)

Robject_dir7 <- paste0(out_path, "GIN_driver/CCNE_vs_HRD/Rdata/")
system(paste0("mkdir -p ", Robject_dir7))
table_dir7 <- paste0(out_path, "GIN_driver/CCNE_vs_HRD/tables/")
system(paste0("mkdir -p ", table_dir7))
plot_dir7 <- paste0(out_path, "GIN_driver/CCNE_vs_HRD/plots/")
system(paste0("mkdir -p ", plot_dir7))

# perform CCNE vs HRD DE:
CCNE_vs_HRD_DE <- do_DE(
  fit_obj = GIN_fit$fit,
  edgeR_obj = GIN_fit$all_edgeR,
  con = c(1, 0, -1, 0),
  descrip = "CCNE_vs_HRD",
  plot_dir = plot_dir7
)

CCNE_vs_HRD_nat_vs_me <- DE_compare(
  DE_res = CCNE_vs_HRD_DE,
  descrip = "CCNE_vs_HRD",
  ref_dir,
  table_dir7,
  plot_dir7
)


#############################################################################
### 6. CCNE1 & HRD vs Unknown GIN driver DE ###
#############################################################################

GIN_annot <- sample_annot

# remove ascites and fallopian tissue:
GIN_annot <- GIN_annot[grep("ascites|fallopian_tissue", GIN_annot$site, invert = T), ]
GIN_counts <- formatted_counts[,colnames(formatted_counts) %in% GIN_annot$ID]

Robject_dir8 <- paste0(out_path, "GIN_driver/unknown_vs_known_GIN/Rdata/")
system(paste0("mkdir -p ", Robject_dir8))
table_dir8 <- paste0(out_path, "GIN_driver/unknown_vs_known_GIN/tables/")
system(paste0("mkdir -p ", table_dir8))
plot_dir8 <- paste0(out_path, "GIN_driver/unknown_vs_known_GIN/plots/")
system(paste0("mkdir -p ", plot_dir8))


# merge CCNE and HRD groups:
unknown_vs_known_GIN_annot <- GIN_annot
unknown_vs_known_GIN_annot$GIN_driver[
  unknown_vs_known_GIN_annot$GIN_driver != "Unknown" &
  unknown_vs_known_GIN_annot$GIN_driver != "Fallopian_tissue"
] <- "Known"

unknown_vs_known_GIN_fit <- fit_glm(
  count_df = GIN_counts,
  sample_annot = unknown_vs_known_GIN_annot,
  repeat_symbols = repeat_symbols,
  cols = col_pal,
  div_type = "GIN_driver",
  Robject_dir = Robject_dir8,
  plot_dir = plot_dir8,
  func_dir
)

unknown_vs_known_GIN_DE <- do_DE(
  fit_obj = unknown_vs_known_GIN_fit$fit,
  edgeR_obj = unknown_vs_known_GIN_fit$all_edgeR,
  con = c(0, -1, 1),
  descrip = "unknown_vs_known",
  plot_dir = plot_dir8
)

# plot retrotransposon DE:
plot_DE(
  DE_results = unknown_vs_known_GIN_DE,
  DE_name = "unknown_vs_known_GIN",
  repeat_genes = RT,
  gene_type = "retrotransposon",
  table_dir = table_dir8,
  plot_dir = plot_dir8,
  FDR_lim = 0.05,
  FC_lim = 0.7,
  num_label = 10,
  manual_lab = "none",
  dot_col = dot_col,
  label_col = label_col,
  up_ctl = "none",
  up_ctl_col = "none"
)


############################################################################
### 7. resistant vs sensitive GIN driver DE ###
#############################################################################

GIN_annot <- sample_annot

# remove ascites and fallopian tissue:
GIN_annot <- GIN_annot[grep("ascites|fallopian_tissue", GIN_annot$site, invert = T), ]
GIN_counts <- formatted_counts[,colnames(formatted_counts) %in% GIN_annot$ID]

Robject_dir9 <- paste0(out_path, "drug_response/resistant_vs_sensitive/Rdata/")
system(paste0("mkdir -p ", Robject_dir9))
table_dir9 <- paste0(out_path, "drug_response/resistant_vs_sensitive/tables/")
system(paste0("mkdir -p ", table_dir9))
plot_dir9 <- paste0(out_path, "drug_response/resistant_vs_sensitive/plots/")
system(paste0("mkdir -p ", plot_dir9))

# remove ascites and fallopian tissue:
drug_annot <- sample_annot[grep("ascites|fallopian_tissue", sample_annot$site, invert = T), ]
drug_counts <- formatted_counts[,colnames(formatted_counts) %in% drug_annot$ID]

# check order of annotation and counts dfs:
print(
  paste0(
    "Sample annotation in the same order as counts file? ",
    identical(drug_annot$ID, colnames(drug_counts))
  )
)

drug_fit <- fit_glm(
  count_df = drug_counts,
  sample_annot = drug_annot,
  repeat_symbols = repeat_symbols,
  cols = col_pal,
  div_type = "treatment_response",
  Robject_dir = Robject_dir9,
  plot_dir = plot_dir9,
  func_dir
)

resistant_vs_sensitive_DE <- do_DE(
  fit_obj = drug_fit$fit,
  edgeR_obj = drug_fit$all_edgeR,
  con = c(0, 1, -1),
  descrip = "resistant_vs_sensitive",
  plot_dir = plot_dir9
)

# plot retrotransposon DE:
plot_DE(
  DE_results = resistant_vs_sensitive_DE,
  DE_name = "resistant_vs_sensitive_GIN",
  repeat_genes = RT,
  gene_type = "retrotransposon",
  table_dir = table_dir9,
  plot_dir = plot_dir9,
  FDR_lim = 0.05,
  FC_lim = 0.7,
  num_label = 10,
  manual_lab = "none",
  dot_col = dot_col,
  label_col = label_col,
  up_ctl = "none",
  up_ctl_col = "none"
)


#############################################################################
### 7. HRD vs CCNE & Unknown GIN driver DE ###
#############################################################################

GIN_annot <- sample_annot

# remove ascites and fallopian tissue:
GIN_annot <- GIN_annot[grep("ascites|fallopian_tissue", GIN_annot$site, invert = T), ]
GIN_counts <- formatted_counts[,colnames(formatted_counts) %in% GIN_annot$ID]

Robject_dir10 <- paste0(out_path, "GIN_driver/hrd_vs_rest_GIN/Rdata/")
system(paste0("mkdir -p ", Robject_dir10))
table_dir10 <- paste0(out_path, "GIN_driver/hrd_vs_rest_GIN/tables/")
system(paste0("mkdir -p ", table_dir10))
plot_dir10 <- paste0(out_path, "GIN_driver/hrd_vs_rest_GIN/plots/")
system(paste0("mkdir -p ", plot_dir10))


# merge CCNE and HRD groups:
hrd_vs_rest_GIN_annot <- GIN_annot
hrd_vs_rest_GIN_annot$GIN_driver[
  hrd_vs_rest_GIN_annot$GIN_driver != "HRD" &
  hrd_vs_rest_GIN_annot$GIN_driver != "Fallopian_tissue"
] <- "CCNE_and_Unknown"

hrd_vs_rest_GIN_fit <- fit_glm(
  count_df = GIN_counts,
  sample_annot = hrd_vs_rest_GIN_annot,
  repeat_symbols = repeat_symbols,
  cols = col_pal,
  div_type = "GIN_driver",
  Robject_dir = Robject_dir10,
  plot_dir = plot_dir10,
  func_dir
)

hrd_vs_rest_GIN_DE <- do_DE(
  fit_obj = hrd_vs_rest_GIN_fit$fit,
  edgeR_obj = hrd_vs_rest_GIN_fit$all_edgeR,
  con = c(-1, 1),
  descrip = "hrd_vs_rest",
  plot_dir = plot_dir10
)

# plot retrotransposon DE:
plot_DE(
  DE_results = hrd_vs_rest_GIN_DE,
  DE_name = "hrd_vs_rest_GIN",
  repeat_genes = sat,
  gene_type = "satellite_fdr_0.1",
  table_dir = table_dir10,
  plot_dir = plot_dir10,
  FDR_lim = 0.1,
  FC_lim = 0.7,
  num_label = 10,
  manual_lab = "none",
  dot_col = dot_col,
  label_col = label_col,
  up_ctl = "none",
  up_ctl_col = "none"
)



