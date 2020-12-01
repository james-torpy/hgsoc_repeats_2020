
home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/hgsoc_repeats/RNA-seq/")
ref_dir <- paste0(project_dir, "/refs/")
star_dir <- paste0(project_dir, "results/star/GC/")

mdss_bams <- read.table(
  paste0(ref_dir, "mdss_bams.txt"),
)

mdss_fqs <- read.table(
  paste0(ref_dir, "mdss_fqs.txt"),
)

bams <- gsub(
  ".sub",
  "",
  gsub(
    "/.*$",
    "",
    gsub(
      "^.*/AOCS",
      "AOCS",
      list.files(
        star_dir, 
        pattern = "Aligned.out.bam", 
        recursive = T,
        full.names = T
      )
    )
  )
)

bams_to_process <- as.character(mdss_bams$V1)[
  !(as.character(mdss_bams$V1) %in% bams)
]
fqs_to_process <- as.character(mdss_fqs$V1)[
  !(as.character(mdss_fqs$V1) %in% bams)
]
fqs_to_process <- fqs_to_process[!(fqs_to_process %in% bams_to_process)]

