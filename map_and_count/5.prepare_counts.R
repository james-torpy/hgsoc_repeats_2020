
project_name <- "hgsoc_repeats/RNA-seq"

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/")
results_dir <- paste0(project_dir, "results/")
in_dir <- paste0(results_dir, "htseq/")

sub <- FALSE


########################################################################
### 0. Determine sample_names and load counts ###
########################################################################

# fetch files and sample names:
filenames <- list.files(
  in_dir, 
  pattern = "counts.gc.custom3.htseq.txt", 
  full.names = T,
  recursive = T
)

# keep only files with size > 0:
file_sizes <- file.info(filenames)$size
filenames <- filenames[file_sizes > 0]

sample_names <- gsub(
  "/.*$", "",
  gsub(
    "^.*AOCS", "AOCS", filenames
  )
)

if (sub) {
  sample_names <- sample_names[grep("sub", sample_names)]
} else {
  sample_names <- sample_names[grep("sub", sample_names, invert = T)]
}

sample_names <- as.list(sample_names)

counts_list <- lapply(sample_names, function(x) {

  print(paste0("Loading ", x, "..."))

  counts <- read.table(
    paste0(in_dir, x, "/counts.gc.custom3.htseq.txt"),
    sep = "\t"
  )

  # remove no feature, ambiguous reads counts etc:
  print("Removing no feature reads...")
  counts <- counts[grep("__", counts$V1, invert = T),]

  # define rownames:
  rownames(counts) <- counts$V1
  counts <- subset(counts, select = V2)

  return(counts)

})


########################################################################
### 1. Combine counts and save ###
########################################################################

# combine into one df:
count_df <- do.call("cbind", counts_list)
colnames(count_df) <- sample_names

# remove samples with no counts:
count_df <- count_df[,colSums(count_df) > 0]

# write as table:
if (sub) {
  out_file <- paste0(in_dir, "sub_combined_counts.txt")
} else {
  out_file <- paste0(in_dir, "combined_counts.txt")
}

write.table(
  count_df,
  out_file,
  sep = "\t",
  quote = F,
  row.names = T,
  col.names = T
)

