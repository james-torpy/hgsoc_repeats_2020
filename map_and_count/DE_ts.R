
project_name <- "hgsoc_repeats/RNA-seq-final"

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/")
results_dir <- paste0(project_dir, "results/")
in_dir <- paste0(results_dir, "htseq/Rdata/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/map_and_count/functions/")
genome_dir <- paste0(project_dir, "genome/")
col_dir <- paste0(home_dir, "R/colour_palettes/")

out_dir <- paste0(results_dir, "/DE_own_vs_nat/")

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
library(ShortRead, lib.loc=lib_loc)
library(EDASeq, lib.loc=lib_loc)
library(RUVSeq, lib.loc=lib_loc)

# load colours:
col_pal <- read.table(
  paste0(col_dir, "colour_palette_1.txt"),
  comment.char = "",
  stringsAsFactors = F
)[,1]


########################################################################
### 1.  Load counts and convert ens ids to symbols ###
########################################################################

# fetch aggregated counts:
all_counts <- readRDS(
  paste0(in_dir, "all_combined_counts_ensembl_ids.Rdata")
)

# replace ensembl ids with symbols:
# split repeats and GC apart:
all_counts$split <- "repeat"
all_counts$split[grep("ENST", rownames(all_counts))] <- "GC"
split_df <- split(all_counts, all_counts$split)
split_df <- lapply(split_df, function(x) subset(x, select = -split))

# keep only lncRNAs, miRNAs and protein coding transcripts:
transcript_annot <- read.table(
  paste0(genome_dir, "gencode.v35.basic.exon.info.txt"),
  fill = T,
  sep = "\t",
  header = T
)

keep_annot <- transcript_annot[
  transcript_annot$transcript_type == "protein_coding" | 
  transcript_annot$transcript_type == "lncRNA" | 
  transcript_annot$transcript_type == "miRNA",
]

print(
  paste0(
    "No GC transcripts before filtering for protein coding, lncRNA, miRNA: ",
    nrow(split_df[[1]])
  )
)
split_df[[1]] <- split_df[[1]][
  rownames(split_df[[1]]) %in% keep_annot$transcript_id,
]
print(
  paste0(
    "No transcripts after filtering: ",
    nrow(split_df[[1]])
  )
)

# replace ensembl ids with symbols:
split_df[[1]]$id <- rownames(split_df[[1]])
split_df[[1]] <- split_df[[1]][
  split_df[[1]]$id %in% transcript_annot$transcript_id,
]

m <- match(split_df[[1]]$id, transcript_annot$transcript_id)
split_df[[1]]$id <- transcript_annot$symbol[m]

# aggregate counts from same genes:
split_df[[1]] <- aggregate(.~id, split_df[[1]], sum)

# make rownames id column:
#split_df[[1]] <- split_df[[1]][!duplicated(split_df[[1]]$id),]
rownames(split_df[[1]]) <- split_df[[1]]$id
split_df[[1]] <- subset(split_df[[1]], select = -id)

# rebind GC with repeats, remove id column and gene type from rownames:
symbol_df <- do.call("rbind", split_df)

rownames(symbol_df)[rownames(symbol_df)=="repeat.MSR1"] <- "MSR1_repeat"
rownames(symbol_df) <- gsub(
  "repeat\\.", 
  "", 
  gsub("GC\\.", "", rownames(symbol_df))
)

count_df <- symbol_df


########################################################################
### 2. Isolate HRD and CCNE1 samples ###
########################################################################

# fetch sample annot:
sample_annot <- read.table(
  paste0(ref_dir, "sample_annot.txt"),
  sep = "\t",
  header = T,
  stringsAsFactors = F
)

# keep counts samples only:
sample_annot <- sample_annot[sample_annot$ID %in% colnames(count_df),]

# get rid of FT:
sample_annot <- sample_annot[sample_annot$site != "fallopian_tissue",]

# fetch CCNE samples:
CCNE_annot <- sample_annot[sample_annot$GIN_driver == "CCNE",]
CCNE_samples <- CCNE_annot$ID

# fetch HRD samples:
HRD_annot <- sample_annot[sample_annot$GIN_driver == "HRD",]
HRD_samples <- HRD_annot$ID

# bind to make sample_annot:
sample_annot <- rbind(CCNE_annot, HRD_annot)

# isolate above samples to count df:
count_df <- count_df[
  ,colnames(count_df) %in% c(CCNE_samples, HRD_samples)
]

# order count_df as in sample_annot:
count_df <- count_df[ ,match(sample_annot$ID, colnames(count_df))]


########################################################################
### 3.  Filtering  ###
########################################################################

prefilter_y <- DGEList(
  counts=count_df, 
  group = sample_annot$GIN_driver
)

# create pre-filter MDS plot:
plot_mds <- function(DGE, prefix, col_pal, plot_dir) {

  mds <- plotMDS(
    DGE,
    labels = DGE$samples$group, 
    method = "bcv", 
    col = as.numeric(DGE$samples$group)
  )
  dev.off()
  plot_df <- merge(
    DGE$samples,
    as.data.frame(mds[[3]]),
    by=0
  )
  plot_df <- subset(plot_df, select = -c(lib.size, norm.factors))
  colnames(plot_df) <- c("sample", "group", "x", "y")

  # plot MDS:
  print(paste0("Plotting ", prefix, " MDS..."))

  p <- ggplot(plot_df, aes(x=x, y=y, colour=group))
  p <- p + geom_point()
  p <- p + labs(
    x = "BCV distance 1", 
    y = "BCV distance 2", 
    colour = "GIN driver"
  )
  p <- p + scale_colour_manual(
    labels = c("CCNE amp.", "HRD"), values = col_pal
  )
  p <- p + labs(colour = "GIN driver")

  pdf(
    paste0(plot_dir, prefix, "_mds.pdf"),
    width = 10
  )
    print(p)
  dev.off()
  
  png(
    paste0(plot_dir, prefix, "_mds.png"),
    height = 7,
    width = 10,
    units = "in",
    res = 300
  )
    print(p)
  dev.off()

}

plot_mds(
  DGE = prefilter_y, 
  prefix = "prefilter", 
  col_pal, 
  plot_dir
)

print("Filtering, calculating normalisation factors... ")

# filter out genes with less than 2 counts > 1:
print(
  paste0(
    "No. genes before filtering: ", 
    nrow(prefilter_y$counts)
  )
)

keep <- rowSums(cpm(count_df)>1) >= ncol(count_df)/3
keep <- keep[names(keep) %in% rownames(prefilter_y$counts)]
filtered_y <- prefilter_y[keep,]

print(
  paste0(
    "No. genes after filtering: ", nrow(filtered_y$counts)
  )
)

# adjust library size:
filtered_y$samples$lib.size <- colSums(filtered_y$counts)
filtered_y$samples

# plot MDS
plot_mds(
  DGE = filtered_y, 
  prefix = "postfilter", 
  col_pal, 
  plot_dir
)


########################################################################
### 4.  Normalisation using RUV-seq  ###
########################################################################

p_data <- data.frame(
  GIN_driver = factor(sample_annot$GIN_driver),
  row.names = sample_annot$ID
)

set <- newSeqExpressionSet(
  filtered_y$counts, 
  phenoData = p_data
)

# create pre-norm RLE plot:
print("Creating pre-norm RLE plot...")
pdf(file = paste0(plot_dir, "pre_norm_RLE.pdf"))
  par(mar=c(1,1,1,1))
  plotRLE(set)
dev.off()

# create RUVseq pre-norm PCA:
print("Creating pre-norm PCA plot...")
pdf(file = paste0(plot_dir, "pre_norm_PCA.pdf"))
  par(mar=c(1,1,1,1))
  plotPCA(set)
dev.off()

# perform between lane full normalisation using RUVSeq:
norm_set <- betweenLaneNormalization(set, which="full")

# create post-norm RLE plot:
print("Creating post-norm RLE plot...")
pdf(file = paste0(plot_dir, "post_norm_RLE.pdf"))
  par(mar=c(1,1,1,1))
  plotRLE(norm_set)
dev.off()

# create RUVseq post-norm PCA:
print("Creating post-norm PCA plot...")
pdf(file = paste0(plot_dir, "post_norm_PCA.pdf"))
  par(mar=c(1,1,1,1))
  plotPCA(norm_set)
dev.off()

# convert to DGE object:
RUV_norm <- DGEList(
  counts = norm_set@assayData$normalizedCounts,
  group = sample_annot$GIN_driver
)

# plot MDS:
plot_mds(
  DGE = RUV_norm, 
  prefix = "post_RUVseq_normalisation", 
  col_pal, 
  plot_dir
)


########################################################################
### 5.  Normalisation using edgeR  ###
########################################################################

# determine normalisation factors:
norm_y <- calcNormFactors(filtered_y)
norm_y$samples

plot_mds(
  DGE = norm_y,
  prefix = "post_edgeR_normalisation",
  col_pal,
  plot_dir
)


########################################################################
### 5.  Dispersion estimation/GLM using edgeR  ###
########################################################################

both_DGE <- list(
  edgeR = norm_y,
  RUVSeq = RUV_norm
)

for (i in 1:length(both_DGE)) {

  # estimate dispersion:
  print(
    paste0(
      "Estimating dispersion for ", names(both_DGE)[i], " data..."
    )
  )
  disp_y <- estimateCommonDisp(both_DGE[[i]])
  disp_y <- estimateTagwiseDisp(disp_y)
  
  # plot tagwise bcv:
  print(
    paste0(
      "Plotting dispersion vs avg. log CPM for ", 
      names(both_DGE)[i], 
      " data..."
    )
  )
  png(paste0(plot_dir, names(both_DGE)[i], "_bcv_vs_cpm.png"))
    plotBCV(disp_y)
  dev.off()
  
  # estimate and compare dispersion by GLM:
  print(paste0("Generating GLM for ", names(both_DGE)[i], " data..."))
  design <- model.matrix(~0 + both_DGE[[i]]$samples$group)
  colnames(design) <- levels(both_DGE[[i]]$samples$group)
  glm_y <- estimateGLMCommonDisp(both_DGE[[i]], design)
  glm_y <- estimateGLMTrendedDisp(glm_y, design)
  glm_y <- estimateGLMTagwiseDisp(glm_y, design)
  
  # plot tagwise bcv:
  print(
    paste0("Plotting GLM on dispersion vs avg. log CPM for ", 
    names(both_DGE)[i], " data...")
  )
  png(paste0(plot_dir, names(both_DGE)[i], "_bcv_vs_cpm_glm.png"))
    plotBCV(glm_y)
  dev.off()
  
  postfilter_edgeR <- list(
    design = design,
    non_GLM = disp_y,
    GLM = glm_y
  )
  
  # fit data to model:
  fit <- glmFit(postfilter_edgeR$GLM, postfilter_edgeR$design)
  
  
  ########################################################################
  ### 6.  Perform DE between CCNE and HRD ###
  ########################################################################
  
  # perform DE:
  con <- c(1, -1)
  
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
      names(both_DGE)[i], "_CCNE1_vs_HRD_smearplot.png"
    )
  )
    plotSmear(postfilter_DE, de.tags=tags)
    abline(h = c(-0.7, 0.7), col = "blue")
  dev.off()
  
  if (i==1) {
    CCNE_vs_HRD_DEs <- list(
      as.data.frame(topTags(postfilter_DE, n=Inf))
    )
  } else {
    CCNE_vs_HRD_DEs[[i]] <- as.data.frame(
      topTags(postfilter_DE, n=Inf)
    )
  }
  
}

names(CCNE_vs_HRD_DEs) <- names(both_DGE)


########################################################################
### 5.  Compare my results with Nature results ###
########################################################################

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

for (i in 1:length(CCNE_vs_HRD_DEs)) {

  # add logFC and p-vals from my DE data:
  both_CCNE_vs_HRD <- merge(nat_CCNE_vs_HRD, CCNE_vs_HRD_DEs[[i]], by = 0)
  colnames(both_CCNE_vs_HRD)[5:9] <- paste0(
    "my_", colnames(both_CCNE_vs_HRD)[5:9]
  )
  
  # report significant genes and save table:
  print(
    paste0(
      "No. significant reported CCNE vs HRD genes is ",
      length(both_CCNE_vs_HRD$my_FDR[both_CCNE_vs_HRD$my_FDR < 0.05]),
      " out of ",
      nrow(both_CCNE_vs_HRD), " for ", names(both_DGE)[i], " data..."
    )
  )
  write.table(
    subset(both_CCNE_vs_HRD, select = -c(my_logCPM, my_LR)),
    paste0(
      table_dir, 
      names(CCNE_vs_HRD_DEs)[i], 
      "_CCNE1_vs_HRD_nature_vs_my_DE.txt"
    ),
      sep = "\t",
      quote = F,
      row.names = F,
      col.names = T
  )
  
  # label significant genes:
  both_CCNE_vs_HRD$sig <- FALSE
  both_CCNE_vs_HRD$sig[both_CCNE_vs_HRD$my_FDR < 0.05] <- TRUE
  both_CCNE_vs_HRD$sig <- factor(both_CCNE_vs_HRD$sig)
  
  # plot nature vs my DE on density scatter:
  p <- ggplot(both_CCNE_vs_HRD, aes(x=nat_logFC, y=my_logFC, colour = sig))
  p <- p + geom_point() 
  p <- p + scale_color_manual(labels = c("non-significant", "significant"), values = c("grey", "#7CBA61"))
  p <- p + theme_cowplot(12)
  p <- p + labs(
    x="Patch et. al. logFC", 
    y="My DE method logFC",
    color = ""
  )
  p <- p + theme(
    text = element_text(size = 12)
  )
  
  png(
    paste0(
      plot_dir, 
      names(CCNE_vs_HRD_DEs)[i], 
      "_CCNE1_vs_HRD_nature_vs_my_DE.png"
    ),
    width = 10,
    height = 7,
    unit = "in",
    res = 300
  )
    print(p)
  dev.off()

  # calculate correlation:
  cor.estimate(both_CCNE_vs_HRD$my_logFC, both_CCNE_vs_HRD$nat_logFC)

}





