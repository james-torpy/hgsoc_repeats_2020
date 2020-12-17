
project_name <- "hgsoc_repeats/RNA-seq-final"

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/")
results_dir <- paste0(project_dir, "results/")
in_dir <- paste0(results_dir, "htseq/")
genome_dir <- paste0(project_dir, "genome/")

sub <- FALSE

library(org.Hs.eg.db)


########################################################################
### 0. Determine sample_names and load repeat counts ###
########################################################################

# fetch file and sample names for both types of counts:
filetypes <- list("gc", "repeats")

both_sample_names <- lapply(filetypes, function(x) {

  # fetch repeat files:
  filenames <- list.files(
    in_dir, 
    pattern = paste0("counts.", x, ".htseq.txt"), 
    full.names = T,
    recursive = T
  )
  
  # keep only files with size > 0:
  file_sizes <- file.info(filenames)$size
  filenames <- filenames[file_sizes > 0]
  
  print(paste0("The following ", x, " files with size = 0b were removed:"))
  print(filenames[file_sizes == 0])
  
  return(
    gsub(
      "/.*$", "",
      gsub(
        "^.*AOCS", "AOCS", filenames
      )
    )
  )

})

# only keep sample_names with both count files:
sample_names <- both_sample_names[[1]][
  both_sample_names[[1]] %in% both_sample_names[[2]]
]

if (sub) {
  sample_names <- sample_names[grep("sub", sample_names)]
} else {
  sample_names <- sample_names[grep("sub", sample_names, invert = T)]
}

sample_names <- as.list(sample_names)


########################################################################
### 1. Load repeat and Gencode counts ###
########################################################################

counts_list <- lapply(sample_names, function(x) {

  both_counts <- lapply(filetypes, function(y) {

    print(paste0("Loading ", y, " counts for ", x, "..."))

    counts <- read.table(
      paste0(in_dir, x, "/counts.", y, ".htseq.txt"),
      sep = "\t"
    )
  
    # remove no feature, ambiguous reads counts etc:
    print("Removing no feature reads...")
    counts <- counts[grep("__", counts$V1, invert = T),]

    # define rownames:
    rownames(counts) <- counts$V1

    return(subset(counts, select = V2))

  })

  return(do.call("rbind", both_counts))

})


########################################################################
### 2. Combine counts, change to gene symbols and save ###
########################################################################

# combine into one df:
count_df <- do.call("cbind", counts_list)
colnames(count_df) <- sample_names

# remove samples with no counts:
count_df <- count_df[,colSums(count_df) > 0]

# save non-aggregated counts:
saveRDS(
  count_df, 
  paste0(Robject_dir, "non_aggregated_combined_counts_ensembl_ids.Rdata")
)

# keep only lncRNAs, miRNAs and protein coding transcripts:
exon_annot <- read.table(
  paste0(genome_dir, "gencode.v35.basic.exon.info.txt"),
  fill = T,
  sep = "\t",
  header = T
)












## load ensembl and entrez to symbol dfs:
#ensembl_ids <- toTable(org.Hs.egENSEMBL)
#symbol_ids <- toTable(org.Hs.egSYMBOL)

# load exon annotation:
exon_annot <- read.table(
  paste0(genome_dir, "gencode.v35.basic.exons.pc.lncRNA.txt"),
  fill = T,
  sep = "\t",
  header = T
)
exon_annot$exon_id <- gsub("\\:.*$", "", exon_annot$exon_id)
exon_annot$exon_id <- gsub("\\..*$", "", exon_annot$exon_id)
exon_annot <- exon_annot[!duplicated(exon_annot$exon_id),]

# replace ensembl ids with symbols:
# split repeats and GC apart:
count_df$split <- "repeat"
count_df$split[grep("ENST", rownames(count_df))] <- "GC"
split_df <- split(count_df, count_df$split)

# replace ensembl ids with symbols:
split_df[[1]]$id <- rownames(split_df[[1]])
split_df[[1]] <- split_df[[1]][
  split_df[[1]]$id %in% exon_annot$exon_id,
]

m <- match(split_df[[1]]$id, exon_annot$exon_id)
split_df[[1]]$id <- exon_annot$symbol[m]

# remove duplicates and make rownames id column:
split_df[[1]] <- split_df[[1]][!duplicated(split_df[[1]]$id),]
rownames(split_df[[1]]) <- split_df[[1]]$id
split_df[[1]] <- subset(split_df[[1]], select = -id)

# rebind GC with repeats, remove id column and gene type from rownames:
symbol_df <- do.call("rbind", split_df)
symbol_df <- subset(symbol_df, select = -split)

rownames(symbol_df)[rownames(symbol_df)=="repeat.MSR1"] <- "MSR1_repeat"
rownames(symbol_df) <- gsub(
  "repeat\\.", 
  "", 
  gsub("GC\\.", "", rownames(symbol_df))
)

# write both as table:
if (sub) {
  out_file_ensembl <- paste0(in_dir, "sub_combined_counts_ensembl_ids.txt")
  out_file_symbol <- paste0(in_dir, "sub_combined_counts_symbols.txt")
} else {
  out_file_ensembl <- paste0(in_dir, "combined_counts_ensembl_ids.txt")
  out_file_symbol <- paste0(in_dir, "combined_counts_symbols.txt")
}

write.table(
  count_df,
  out_file_ensembl,
  sep = "\t",
  quote = F,
  row.names = T,
  col.names = T
)

write.table(
  symbol_df,
  out_file_symbol,
  sep = "\t",
  quote = F,
  row.names = T,
  col.names = T
)

