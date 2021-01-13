
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
incl_ascites <- FALSE

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

lib_loc <- paste0(home_dir, "/R/3.6.0")

library(edgeR)
library(DESeq2, lib=lib_loc)
library(ggplot2)
library(cowplot)
library(ggpointdensity, lib.loc=lib_loc)
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




#############################################################################
### 3. Perform DE between HGSOC and FT ###
#############################################################################

if (!file.exists(paste0(Robject_dir, "HGSOC_vs_FT_edgeR.Rdata"))) {

  HGSOC_vs_FT_DGE <- DGEList(
    counts=formatted_counts, 
    group = sample_annot$site
  )
  
  HGSOC_vs_FT_y <- edgeR_norm_factors(
    DGE_object = HGSOC_vs_FT_DGE, 
    prefix = "all",
    div_type = "site",
    dot_col = col_pal,
    plot_dir
  )

  saveRDS(HGSOC_vs_FT_y, paste0(Robject_dir, "HGSOC_vs_FT_edgeR.Rdata"))

} else {
  HGSOC_vs_FT_y <- readRDS(paste0(Robject_dir, "HGSOC_vs_FT_edgeR.Rdata"))
}

# fit data to model:
fit <- glmFit(HGSOC_vs_FT_y$GLM, HGSOC_vs_FT_y$design)

# perform DE:
HGSOC_vs_FT <- glmLRT(
  fit, 
  contrast=c(-1, 1)
)

# check DE numbers:
sig <- decideTestsDGE(
  HGSOC_vs_FT, adjust.method = "BH", p.value = 0.05
)
print(summary(sig))

# generate smearplot as sanity check :
tags <- rownames(HGSOC_vs_FT_y$GLM)[as.logical(sig)]
png(paste0(plot_dir, "CCNE_vs_HRD_smearplot.png"))
  plotSmear(HGSOC_vs_FT, de.tags=tags)
  abline(h = c(-1, 1), col = "blue")
dev.off()

# fetch all DE and determine best CPM threshold (one with most sig genes)
HGSOC_vs_FT_df <- as.data.frame(topTags(HGSOC_vs_FT, n=Inf))


#############################################################################
### 4. Plot DE non-repeat genes ###
#############################################################################

# load repeat annotation:
if (!file.exists(paste0(Robject_dir, "all_repeat_genes.Rdata"))) {

  repeat_gtf <- read.table(
    paste0(genome_dir, "custom3.repeats.hg38.gtf"),
    sep = "\t",
    header = F,
    fill = T
  )
  saveRDS(repeat_gtf, paste0(Robject_dir, "all_repeat_genes.Rdata"))

} else {

  repeat_gtf <- readRDS(paste0(Robject_dir, "all_repeat_genes.Rdata"))
  repeat_genes <- gsub("ID ", "", as.character(unique(repeat_gtf$V9)))
  
}

# plot non-repeat DEs:
plot_DE(
  DE_results = HGSOC_vs_FT_df,
  DE_name = "HGSOC_vs_FT",
  repeat_genes = repeat_genes,
  gene_type = "non_repeat",
  table_dir,
  plot_dir,
  FDR_lim = 0.05,
  FC_lim = 0.7,
  num_label = 20,
  manual_lab = "none",
  dot_col = "#8BCCB5",
  label_col = "#430F82"
)
# plot repeat DEs:
plot_DE(
  DE_results = HGSOC_vs_FT_df,
  DE_name = "HGSOC_vs_FT",
  repeat_genes = repeat_genes,
  gene_type = "repeat",
  table_dir,
  plot_dir,
  FDR_lim = 0.05,
  FC_lim = 0.7,
  num_label = 20,
  manual_lab = "none",
  dot_col = "#8BCCB5",
  label_col = "#430F82"
) 