

project_name <- "hgsoc_repeats/RNA-seq"

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/")
results_dir <- paste0(project_dir, "results/")
in_dir <- paste0(results_dir, "htseq/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/map_and_count/functions/")
genome_dir <- paste0(project_dir, "genome/")

sub <- FALSE
RNA_type <- "repeat"
div_type <- "GIN_driver"

out_dir <- paste0(results_dir, "DE/", RNA_type, "_DE_by_", div_type, "/" )
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

# volcano plot function:
gen_plot <- dget(paste0(func_dir, "gen_plot.R"))


########################################################################
### 1.  Load counts and annotations ###
########################################################################

# fetch aggregated counts:
all_counts <- read.table(
  paste0(in_dir, "combined_counts.txt"),
  sep = "\t",
  header = T
)
colnames(all_counts) <- gsub("\\.", "-", colnames(all_counts))

# fetch sample annot and keep counts samples:
sample_annot <- read.table(
  paste0(ref_dir, "sample_annot.txt"),
  sep = "\t",
  header = T,
  stringsAsFactors = F
)
sample_annot$treatment_response <- gsub(
  " ", 
  "_", 
  sample_annot$treatment_response
)

sample_annot <- sample_annot[sample_annot$ID %in% colnames(all_counts),]

# if not comparing primary tumour to ascites or met, keep only tumour:
sample_annot <- sample_annot[sample_annot$site == "Tumour",]
all_counts <- all_counts[,colnames(all_counts) %in% sample_annot$ID]


#############################################################################
### 1. Use edgeR on raw counts to choose optimum filtering  ###
#############################################################################

# create DGEList object for edgeR:
prefilter_y <- DGEList(
  counts=all_counts, 
  group = eval(parse(text=paste0("sample_annot$", div_type)))
)

# check library size:
prefilter_y$samples$lib.size

if (!file.exists(paste0(Robject_dir, "prefilter_edgeR.Rdata"))) {

  # calculate normalisation factors and check QC plots of different methods:
  prefilter_edgeR <- edgeR_norm_factors(
    DGE_object = prefilter_y, 
    prefix = "prefilter",
    plot_dir
  )

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
  all_counts, 
  prefilter_DE_df, 
  threshs = seq(0, 5, 0.1)
)

# chose CPM filtering value to be 1 as returned most amount of significant 
# genes after FDR adjustment 
# (but no filtering returned even more - is this weird?)


#############################################################################
### 2. Filtering and normalisation using edgeR ###
#############################################################################

# filter out genes with less than 2 counts > 1:
print(paste0("No. genes before filtering: ", nrow(prefilter_y$counts)))

keep <- rowSums(cpm(all_counts)>1) >= 2
filtered_y <- prefilter_y[keep,]

print(paste0("No. genes after filtering: ", nrow(filtered_y$counts)))

# adjust library size:
filtered_y$samples$lib.size <- colSums(filtered_y$counts)
filtered_y$samples

# calculate normalisation factors and check QC plots of different methods:
if (!file.exists(paste0(Robject_dir, "postfilter_edgeR.Rdata"))) {

  postfilter_edgeR <- edgeR_norm_factors(
    filtered_y, 
    prefix = "postfilter",
    plot_dir
  )

  saveRDS(postfilter_edgeR, paste0(Robject_dir, "postfilter_edgeR.Rdata"))

} else {
  postfilter_edgeR <- readRDS(paste0(Robject_dir, "postfilter_edgeR.Rdata"))
}

# fit data to model:
fit <- glmFit(postfilter_edgeR$GLM, postfilter_edgeR$design)
      

#############################################################################
### 4. Perform DE between all groups ###
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
        abline(h = c(-1, 1), col = "blue")
      dev.off()

      if (!exists("all_DEs")) {

        all_DEs <- list(as.data.frame(topTags(prefilter_DE, n=Inf)))

      } else {

        all_DEs <- append(
          all_DEs, list(as.data.frame(topTags(prefilter_DE, n=Inf)))
        )

      }
      names(all_DEs)[length(all_DEs)] <- paste0(all_groups[i], "_vs_", all_groups[j])

    }

  }

}


#############################################################################
### 3. Check DE against Nature paper values ###
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
  both_CCNE_vs_HRD <- merge(nat_CCNE_vs_HRD, DE_df, by = 0)
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
### 4. Plot DE non-repeat genes ###
#############################################################################

# load repeat annotation:
if (!file.exists(paste0(Robject_dir, "all_repeat_genes.Rdata"))) {

  repeat_gtf <- readRDS(paste0(Robject_dir, "all_repeat_genes.Rdata"))

} else {

  repeat_gtf <- read.table(
    paste0(genome_dir, "custom3.repeats.hg38.gtf"),
    sep = "\t",
    header = F,
    fill = T
  )
  repeat_genes <- gsub("ID ", "", as.character(unique(repeat_gtf$V9)))
  saveRDS(repeat_gtf, paste0(Robject_dir, "all_repeat_genes.Rdata"))

}

# isolate non-repeat genes:
non_repeat_DE <- DE_df[!(rownames(DE_df) %in% repeat_genes),]

# prepare data:
DE_labelled <- prep_plot(
  DE_data = non_repeat_DE,
  FDR_lim = 0.05,
  FC_lim = 0.7,
  num_label = 10
)

# save up and down genes as table:
DE_up <- DE_labelled[DE_labelled$logFC > 0,]
DE_up <- DE_up[order(DE_up$logFC),]
head(DE_up)
write.table(
  DE_up,
  paste0(table_dir, "upregulated_non_repeat.txt"),
  sep = "\t",
  col.names = T,
  row.names = T,
  quote = F
)

DE_down <- DE_labelled[DE_labelled$logFC < 0,]
DE_down <- DE_down[order(DE_down$logFC),]
write.table(
  DE_down,
  paste0(table_dir, "downregulated_non_repeat.txt"),
  sep = "\t",
  col.names = T,
  row.names = T,
  quote = F
)

# generate volcano plot:
DE_plot <- gen_plot(
  plot_df = DE_labelled
)

png(paste0(plot_dir, "non_repeat_DE_volcano.png"))
  print(DE_plot)
dev.off()


#############################################################################
### 4. Plot DE repeat RNA ###
#############################################################################

# isolate non-repeat genes:
repeat_DE <- DE_df[rownames(DE_df) %in% repeat_genes,]

# prepare data:
repeat_DE_labelled <- prep_plot(
  DE_data = repeat_DE,
  FDR_lim = 0.05,
  FC_lim = 0.7,
  num_label = "sig"
)

# save up and down genes as table:
repeat_DE_up <- repeat_DE_labelled[repeat_DE_labelled$logFC > 0,]
repeat_DE_up <- repeat_DE_up[order(repeat_DE_up$logFC),]
head(repeat_DE_up)
write.table(
  repeat_DE_up,
  paste0(table_dir, "upregulated_repeats.txt"),
  sep = "\t",
  col.names = T,
  row.names = T,
  quote = F
)

repeat_DE_down <- repeat_DE_labelled[repeat_DE_labelled$logFC < 0,]
repeat_DE_down <- repeat_DE_down[order(repeat_DE_down$logFC),]
write.table(
  repeat_DE_down,
  paste0(table_dir, "downregulated_repeats.txt"),
  sep = "\t",
  col.names = T,
  row.names = T,
  quote = F
)

# generate volcano plot:
repeat_DE_plot <- gen_plot(
  plot_df = repeat_DE_labelled
)

png(paste0(plot_dir, "repeat_DE_volcano.png"))
  print(repeat_DE_plot)
dev.off()





