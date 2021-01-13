
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
sample_annot$type[sample_annot$site == "fallopian_tissue"] <- "fallopian_tissue"
sample_annot$treatment_response[sample_annot$site == "fallopian_tissue"] <- "Fallopian_tissue"
sample_annot$GIN_driver[sample_annot$site == "fallopian_tissue"] <- "Fallopian_tissue"

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

# create DGEList objects for edgeR:
prefilter_list <- list(
  all = DGEList(
    counts=formatted_counts, 
    group = sample_annot$GIN_driver
  ),
  gc = DGEList(
    counts = formatted_counts[!(rownames(formatted_counts) %in% repeat_symbols), ], 
    group = sample_annot$GIN_driver
  ),
  repeats = DGEList(
    counts = formatted_counts[rownames(formatted_counts) %in% repeat_symbols, ], 
    group = sample_annot$GIN_driver
  )
)

# check library size:
prefilter_list$all$samples$lib.size

if (!file.exists(paste0(Robject_dir, "prefilter_edgeR.Rdata"))) {

  # calculate normalisation factors and check QC plots of different methods:
  for (i in 1:length(prefilter_list)) {

    print(
      paste0(
        "Calculating normalisation factors for and MDS plotting ", 
        names(prefilter_list)[i]
      )
    )

    prefilter_edgeR_y <- edgeR_norm_factors(
      DGE_object = prefilter_list[[i]], 
      prefix = paste0("prefilter_", names(prefilter_list)[i]),
      div_type = "GIN_driver",
      dot_col = col_pal,
      plot_dir
    )

    if (i==1) {
      prefilter_edgeR_list <- list(prefilter_edgeR_y)
    } else {
      prefilter_edgeR_list[[i]] <- prefilter_edgeR_y
    }

  }
  names(prefilter_edgeR_list) <- names(prefilter_list)

  prefilter_edgeR <- prefilter_edgeR_list$all
  saveRDS(prefilter_edgeR, paste0(Robject_dir, "prefilter_edgeR.Rdata"))

} else {
  prefilter_edgeR <- readRDS(paste0(Robject_dir, "prefilter_edgeR.Rdata"))
}

# fit data to model:
fit <- glmFit(prefilter_edgeR$GLM, prefilter_edgeR$design)

# perform DE on unfiltered counts, between any groups:
prefilter_DE <- glmLRT(
  fit, 
  contrast=c(-1, 1, rep(0, ncol(prefilter_edgeR$design)-2))
)

# check DE numbers:
sig <- decideTestsDGE(
  prefilter_DE, adjust.method = "BH", p.value = 0.05
)
print(summary(sig))

# generate smearplot as sanity check :
tags <- rownames(prefilter_edgeR$GLM)[as.logical(sig)]
png(paste0(plot_dir, "prefilter_DE_smearplot.png"))
  plotSmear(prefilter_DE, de.tags=tags)
  abline(h = c(-1, 1), col = "blue")
dev.off()

# fetch all DE and determine best CPM threshold (one with most sig genes)
prefilter_DE_df <- as.data.frame(topTags(prefilter_DE, n=Inf))

filter_thresh_test(
  formatted_counts, 
  prefilter_DE_df, 
  threshs = seq(0, 5, 0.1),
  min_no_passed_thresh = 2
)

# chose not to filter by CPM as this gives the most significant results


#############################################################################
### 2. Filtering and normalisation using edgeR ###
#############################################################################

if (!file.exists(paste0(Robject_dir, "postfilter_edgeR.Rdata"))) {

  # calculate normalisation factors and check QC plots of different methods:
  for (i in 1:length(prefilter_list)) {

    print(
      paste0(
        "Filtering, calculating normalisation factors for ", 
        names(prefilter_list)[i]
      )
    )

    # filter out genes with less than 2 counts > 1:
    print(
      paste0(
        "No. genes before filtering ", names(prefilter_list)[i], ": ", 
        nrow(prefilter_list[[i]]$counts)
      )
    )
    
    keep <- rowSums(cpm(formatted_counts)>0.5) >= 2
    keep <- keep[names(keep) %in% rownames(prefilter_list[[i]]$counts)]

    filtered_y <- prefilter_list[[i]][keep,]

    print(
      paste0(
        "No. genes after filtering ", 
        names(prefilter_list)[i], ": ", nrow(filtered_y$counts)
      )
    )

    # adjust library size:
    filtered_y$samples$lib.size <- colSums(filtered_y$counts)
    filtered_y$samples

    if (i==1) {

      postfilter_edgeR_y <- edgeR_norm_factors(
        DGE_object = filtered_y, 
        prefix = paste0("postfilter_", names(prefilter_list)[i]),
        div_type = "GIN_driver",
        dot_col = col_pal,
        plot_dir
      )
    
      postfilter_edgeR <- postfilter_edgeR_y

    } else {

      edgeR_norm_factors(
        DGE_object = filtered_y, 
        prefix = paste0("postfilter_", names(prefilter_list)[i]),
        div_type = "GIN_driver",
        dot_col = col_pal,
        plot_dir,
        mds_only = TRUE
      )

    }

  }

  saveRDS(postfilter_edgeR, paste0(Robject_dir, "postfilter_edgeR.Rdata"))

} else {
  postfilter_edgeR <- readRDS(paste0(Robject_dir, "postfilter_edgeR.Rdata"))
}

# fit data to model:
fit <- glmFit(postfilter_edgeR$GLM, postfilter_edgeR$design)


#############################################################################
### 3. Compare CCNE1 vs HRD DE to Nature paper results ###
#############################################################################

all_groups <- as.list(colnames(postfilter_edgeR$design))

con <- c(1, 0, -1, 0)

# perform DE on unfiltered counts, between any groups:
postfilter_DE <- glmLRT(
  fit, 
  contrast=con
)

# check DE numbers:
sig <- decideTestsDGE(
  postfilter_DE, adjust.method = "BH", p.value = 0.05
)
print(summary(sig))

# generate smearplot as sanity check :
tags <- rownames(postfilter_edgeR$GLM)[as.logical(sig)]
png(
  paste0(
    plot_dir, 
    "CCNE_vs_HRD_smearplot.png"
  )
)
  plotSmear(postfilter_DE, de.tags=tags)
  abline(h = c(-0.7, 0.7), col = "blue")
dev.off()

CCNE_vs_HRD <- as.data.frame(topTags(postfilter_DE, n=Inf))

# fetch Nature DE results:
nat_CCNE_vs_HRD <- read.table(
  paste0(ref_dir, "nat_CCNE_vs_HRD.txt"),
  header = T,
  stringsAsFactors = F
)
rownames(nat_CCNE_vs_HRD) <- nat_CCNE_vs_HRD$Gene
nat_CCNE_vs_HRD <- subset(nat_CCNE_vs_HRD, select=-Gene)
colnames(nat_CCNE_vs_HRD) <- c(
  paste0(
    "nat_", colnames(nat_CCNE_vs_HRD)
  )
)

# add logFC and p-vals from my DE data:
both_CCNE_vs_HRD <- merge(nat_CCNE_vs_HRD, CCNE_vs_HRD, by = 0)
colnames(both_CCNE_vs_HRD)[5:9] <- paste0(
  "my_", colnames(both_CCNE_vs_HRD)[5:9]
)

# label significant genes:
both_CCNE_vs_HRD$sig <- FALSE
both_CCNE_vs_HRD$sig[both_CCNE_vs_HRD$my_FDR < 0.05] <- TRUE
both_CCNE_vs_HRD$sig <- factor(
  both_CCNE_vs_HRD$sig,
  levels = c(TRUE, FALSE)
)

# report significant genes and save table:
print(
  paste0(
    "No. significant reported CCNE vs HRD genes is ",
    length(both_CCNE_vs_HRD$my_FDR[both_CCNE_vs_HRD$my_FDR < 0.05]),
    " out of ",
    nrow(both_CCNE_vs_HRD)
  )
)
write.table(
  subset(both_CCNE_vs_HRD, select = -c(my_logCPM, my_LR)),
  paste0(table_dir, "CCNE1_vs_HRD_nature_vs_my_DE.txt"),
    sep = "\t",
    quote = F,
    row.names = F,
    col.names = T
)

# calculate correlation between my results and Nature:
bm_cor <- cor.test(
  both_CCNE_vs_HRD$my_logFC,
  both_CCNE_vs_HRD$nat_logFC,
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
x_annot_pos <- max(both_CCNE_vs_HRD$my_logFC) + 
  (max(both_CCNE_vs_HRD$nat_logFC)/4)
y_annot_pos <- max(both_CCNE_vs_HRD$my_logFC) + 
  (max(both_CCNE_vs_HRD$my_logFC)/4)

# plot nature vs my DE on density scatter:
p <- ggplot(both_CCNE_vs_HRD, aes(x=nat_logFC, y=my_logFC, colour = sig))
p <- p + geom_point() 
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
  paste0(plot_dir, "CCNE1_vs_HRD_nature_vs_my_DE.png"),
  width = 10,
  height = 7,
  unit = "in",
  res = 300
)
  print(p)
dev.off()

# find genes with most discrepency:
both_CCNE_vs_HRD$diff <- both_CCNE_vs_HRD$nat_logFC - 
  both_CCNE_vs_HRD$my_logFC
disc_df <- both_CCNE_vs_HRD[order(both_CCNE_vs_HRD$diff, decreasing=F),]

top_diff <- head(disc_df$Row.names, 10)

# fetch counts for these genes:
diff_counts <- formatted_counts[top_diff,]

# split by driver:
diff_CCNE <- diff_counts[
  ,colnames(diff_counts) %in% 
    sample_annot$ID[sample_annot$GIN_driver == "CCNE"]
]

diff_HRD <- diff_counts[
  ,colnames(diff_counts) %in% 
    sample_annot$ID[sample_annot$GIN_driver == "HRD"]
]



#############################################################################
### 3. Perform DE between all groups ###
#############################################################################

all_groups <- as.list(colnames(postfilter_edgeR$design))

# for each group, perform DE between that and all others:
for (i in 1:length(all_groups)) {

  for (j in 1:length(all_groups)) {

    if (i != j) {

      con <- rep(0, length(all_groups))
      con[grep(all_groups[i], colnames(postfilter_edgeR$design))] <- -1
      con[grep(all_groups[j], colnames(postfilter_edgeR$design))] <- 1

      # perform DE on unfiltered counts, between any groups:
      postfilter_DE <- glmLRT(
        fit, 
        contrast=con
      )
      
      # check DE numbers:
      sig <- decideTestsDGE(
        postfilter_DE, adjust.method = "BH", p.value = 0.05
      )
      print(summary(sig))
      
      # generate smearplot as sanity check :
      tags <- rownames(postfilter_edgeR$GLM)[as.logical(sig)]
      png(
        paste0(
          plot_dir, 
          all_groups[i], 
          "_vs_", 
          all_groups[j], 
          "_smearplot.png"
        )
      )
        plotSmear(postfilter_DE, de.tags=tags)
        abline(h = c(-0.7, 0.7), col = "blue")
      dev.off()

      if (!exists("all_DEs")) {

        all_DEs <- list(as.data.frame(topTags(postfilter_DE, n=Inf)))

      } else {

        all_DEs <- append(
          all_DEs, list(as.data.frame(topTags(postfilter_DE, n=Inf)))
        )

      }
      names(all_DEs)[length(all_DEs)] <- paste0(all_groups[i], "_vs_", all_groups[j])

    }

  }

}


#############################################################################
### 4. Check DE against Nature paper values ###
#############################################################################

if (div_type == "GIN_driver") {

  # fetch Nature DE results:
  nat_CCNE_vs_HRD <- read.table(
    paste0(ref_dir, "nat_CCNE_vs_HRD.txt"),
    header = T,
    stringsAsFactors = F
  )
  rownames(nat_CCNE_vs_HRD) <- nat_CCNE_vs_HRD$Gene
  nat_CCNE_vs_HRD <- subset(nat_CCNE_vs_HRD, select=-Gene)
  colnames(nat_CCNE_vs_HRD) <- c(
    paste0(
      "nat_", colnames(nat_CCNE_vs_HRD)
    )
  )
  
  # add logFC and p-vals from my DE data:
  both_CCNE_vs_HRD <- merge(nat_CCNE_vs_HRD, all_DEs$CCNE_vs_HRD, by = 0)
  colnames(both_CCNE_vs_HRD)[5:9] <- paste0(
    "my_", colnames(both_CCNE_vs_HRD)[5:9]
  )
  
  # report significant genes and save table:
  print(
    paste0(
      "No. significant reported CCNE vs HRD genes is ",
      length(both_CCNE_vs_HRD$my_FDR[both_CCNE_vs_HRD$my_FDR < 0.05]),
      " out of ",
      nrow(both_CCNE_vs_HRD)
    )
  )
  write.table(
    subset(both_CCNE_vs_HRD, select = -c(my_logCPM, my_LR)),
    paste0(table_dir, "CCNE1_vs_HRD_nature_vs_my_DE.txt"),
      sep = "\t",
      quote = F,
      row.names = F,
      col.names = T
  )
  
  # plot nature vs my DE on density scatter:
  p <- ggplot(both_CCNE_vs_HRD, aes(x=nat_logFC, y=my_logFC))
  p <- p + geom_pointdensity() 
  p <- p + scale_color_viridis_c()
  p <- p + theme_cowplot(12)
  p <- p + labs(
    x="Patch et. al. logFC", 
    y="My DE method logFC",
    color = "No. neighbours"
  )
  p <- p + theme(
    text = element_text(size = 12)
  )
  
  png(
    paste0(plot_dir, "CCNE1_vs_HRD_nature_vs_my_DE.png"),
    width = 10,
    height = 7,
    unit = "in",
    res = 300
  )
    print(p)
  dev.off()

}


#############################################################################
### 5. Plot DE non-repeat genes ###
#############################################################################

# load repeat annotation:
if (!file.exists(paste0(Robject_dir, "all_repeat_genes.Rdata"))) {

    repeat_gtf <- read.table(
    paste0(genome_dir, "custom3.repeats.hg38.gtf"),
    sep = "\t",
    header = F,
    fill = T
  )
  repeat_genes <- gsub("ID ", "", as.character(unique(repeat_gtf$V9)))
  saveRDS(repeat_gtf, paste0(Robject_dir, "all_repeat_genes.Rdata"))

} else {

  repeat_gtf <- readRDS(paste0(Robject_dir, "all_repeat_genes.Rdata"))
  
}

for (d in 1:length(all_DEs)) {

  # plot non-repeat DEs:
  plot_DE(
    DE_results = all_DEs[[d]],
    DE_name = names(all_DEs)[d],
    repeat_genes = repeat_genes,
    gene_type = "non_repeat",
    table_dir,
    plot_dir,
    FDR_lim = 0.05,
    FC_lim = 0.7,
    num_label = 20,
    dot_col = "#8BCCB5",
    label_col = "#430F82"
  ) 

  # plot repeat DEs:
  plot_DE(
    all_DEs[[d]],
    names(all_DEs)[d],
    repeat_genes,
    "repeat",
    table_dir,
    plot_dir,
    FDR_lim = 0.05,
    FC_lim = 0.7,
    num_label = 20,
    dot_col = "#8BCCB5",
    label_col = "#430F82"
  ) 

}


