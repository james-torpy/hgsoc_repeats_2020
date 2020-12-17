
project_name <- "hgsoc_repeats/RNA-seq-final"

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/")
results_dir <- paste0(project_dir, "results/")
in_dir <- paste0(results_dir, "htseq/Rdata/")
star_dir <- paste0(results_dir, "star/")
gc_dir <- paste0(star_dir, "GC/")
ribo_dir <- paste0(star_dir, "ribo/")
genome_dir <- paste0(project_dir, "genome/")
ref_dir <- paste0(project_dir, "refs/")

out_path <- paste0(results_dir, "QC/")
plot_dir <- paste0(out_path, "plots/")
table_dir <- paste0(out_path, "tables/")
Robject_dir <- paste0(out_path, "Rdata/")

system(paste0("mkdir -p ", plot_dir))
system(paste0("mkdir -p ", table_dir))
system(paste0("mkdir -p ", Robject_dir))

library(ggplot2)
library(cowplot)
library(reshape)
library(tibble)
library(grid)

col_pal <- c(
  repeats = "#DDCC77", 
  protein_coding = "#F8766D", 
  non_coding = "#00B0F6", 
  other = "#AA4499", 
  ribosomal = "black"
)

##########################################################################
### 0. Load counts and mapped read nos ###
##########################################################################

# load combines counts:
all_counts <- readRDS(
  paste0(in_dir, "all_combined_counts_ensembl_ids.Rdata")
)

transcript_annot <- read.table(
  paste0(genome_dir, "gencode.v35.basic.exon.info.txt"),
  fill = T,
  sep = "\t",
  header = T
)

repeat_annot <- read.table(
  paste0(genome_dir, "repeats.hg38.info.txt"),
  fill = T,
  sep = "\t",
  header = T
)

# split all_counts:
gc_counts <- all_counts[
  rownames(all_counts) %in% transcript_annot$transcript_id,
]
repeat_counts <- all_counts[
  rownames(all_counts) %in% repeat_annot$symbol,
]

# remove ribosomal transcripts:
gc_counts <- gc_counts[
  !(rownames(gc_counts) %in% transcript_annot$transcript_id[
      grep("ribo", transcript_annot$transcript_type)
  ]),
]

# load mapped read numbers:
mapped_no <- readRDS(
  paste0(in_dir, "all_combined_mapped_read_no.Rdata")
)

# load sample annotation:
sample_annot <- read.table(
  paste0(ref_dir, "sample_annot.txt"),
  sep = "\t",
  header = T,
  stringsAsFactors = F
)


##########################################################################
### 1. Calculate total counts/mapped reads of RNA types  ###
##########################################################################

# calculate total GC counts per sample:
total_gc_counts <- data.frame(apply(gc_counts, 2, sum))

# add total protein coding, non-coding, other counts per sample:
nc_terms <- c(
  "lncRNA",
  "processed_transcript",
  "retained_intron",
  grep("pseudogene", unique(transcript_annot$transcript_type), value = T)
)

pc_counts <- gc_counts[
  rownames(gc_counts) %in% 
    transcript_annot$transcript_id[
      transcript_annot$transcript_type =="protein_coding" 
    ],
]
total_pc_counts <- data.frame(apply(pc_counts, 2, sum))

nc_counts <- gc_counts[
  rownames(gc_counts) %in% 
    transcript_annot$transcript_id[
      transcript_annot$transcript_type %in% nc_terms
    ],
]
total_nc_counts <- data.frame(apply(nc_counts, 2, sum))

other_counts <- gc_counts[
  !(rownames(gc_counts) %in% c(rownames(pc_counts), rownames(nc_counts))),
]
total_other_counts <- data.frame(apply(other_counts, 2, sum))

total_counts <- cbind(
  total_gc_counts,
  total_pc_counts,
  total_nc_counts,
  total_other_counts
)

# add total ribo and gc mapped reads per sample:
total_counts <- merge(total_counts, mapped_no, by=0)
total_counts <- total_counts %>%
  column_to_rownames(var="Row.names")

# add total repeat counts per sample:
total_counts <- merge(
  total_counts, 
  data.frame(apply(repeat_counts, 2, sum)), 
  by=0
)
total_counts <- total_counts %>%
  column_to_rownames(var="Row.names")

colnames(total_counts) <- c(
  "total_gencode", 
  "protein_coding", 
  "non_coding", 
  "other",
  "gencode_mapped",
  "ribosomal",
  "repeats"
)

# keep total gencode and repeat counts vs mapped gc:
counts_vs_mapped <- subset(
  total_counts, 
  select = c("total_gencode", "repeats", "gencode_mapped")
)
counts_vs_mapped$total_count <- counts_vs_mapped$total_gencode + 
  counts_vs_mapped$repeats

# create proportion plot df and calculate totals
RNA_df <- subset(
  total_counts, 
  select = c("repeats", "protein_coding", "non_coding", "other", "ribosomal")
)
RNA_df$total <- apply(RNA_df, 1, sum)


##########################################################################
### 2. Create total mapped vs total counts scatter plot  ###
##########################################################################

p <- ggplot(counts_vs_mapped, aes(x=gencode_mapped, y=total_count))
p <- p + geom_point()
p <- p + theme_cowplot(12)
p <- p + labs(x = "Mapped reads", y="Counted reads")

pdf(paste0(plot_dir, "/counts_vs_mapped.pdf"))
  p
dev.off()

png(paste0(plot_dir, "/counts_vs_mapped.png"))
  p
dev.off()


##########################################################################
### 3. Create composition barplot ###
##########################################################################

# order RNA_df according to sample_annot:
sample_annot <- read.table(
  paste0(ref_dir, "sample_annot.txt"),
  sep = "\t",
  header = T,
  stringsAsFactors = F
)
sample_annot <- sample_annot[
  sample_annot$ID %in% colnames(gc_counts),
]

RNA_df <- RNA_df[sample_annot$ID,]

# calculate proportions of each RNA type:
RNA_prop_df <-subset(RNA_df, select = -total)
RNA_prop_df <- apply(RNA_prop_df, 2, function(x) (x/RNA_df$total)*100)


# add both dfs to list for plotting:
plot_dfs <- list(
  absolute = subset(RNA_df, select = -total),
  proportion = as.data.frame(RNA_prop_df)
)

## set up sample annotation dimensions:
#sample_rle <- rle(sample_annot$site)
#sample_dims <- data.frame(
#  type = sample_rle$values,
#  start = cumsum(
#  	c(0, sample_rle$lengths[1:length(sample_rle$lengths)-1])
#  ),
#  length = sample_rle$lengths
#)
#sample_dims$start <- sample_dims$start/nrow(sample_annot)

# melt and create plots:
for (i in 1:length(plot_dfs)) {
 
  plot_dfs[[i]]$ID <- rownames(plot_dfs[[i]])
  plot_dfs[[i]] <- melt(plot_dfs[[i]])
  
  # create plot:
  p <- ggplot(plot_dfs[[i]], aes(x=ID, y=value))
  p <- p + geom_bar(stat = "identity", aes(fill = variable))
  p <- p + scale_fill_manual(
  	labels = gsub("_", "-", names(col_pal)), 
  	values = col_pal
  )
  if (i==2) {
  	p <- p + scale_y_continuous(expand = c(0,0.1))
  }
  p <- p + labs(x = "Samples", y = "Proportion of total", fill = "RNA type")
  p <- p + theme(
  	text = element_text(size=20),
  	axis.title.x = element_blank(),
  	#axis.text.x = element_text(angle = 90),
  	axis.text.x = element_blank(),
  	axis.ticks.x = element_blank(),
  	axis.text.y = element_text(size=18)
  )

  grob_p <- ggplotGrob(p)

  pdf(
  	paste0(plot_dir, "RNA_", names(plot_dfs)[i], "_barplot.pdf"),
  	width = 15
  )

    grid.newpage()
  
      # plot one viewport per plot:
      pushViewport(viewport(x = 0.53, y = 0.5, width = 0.95, height = 0.85))
        grid.draw(grob_p)
      popViewport()

#      pushViewport(viewport(x = 0.478, y = 0.12, width = 0.706, height = 0.025))
#        
#        grid.rect()
#
#        for (j in 1:nrow(sample_dims)) {
#
#          pushViewport(viewport(
#            x = sample_dims$start[j], 
#            y = 0.5, 
#            width = 0.0101*sample_dims$length[j], 
#            height = 1, 
#            just = "left"
#          ))
#            grid.rect()
#          popViewport()
#
#        }
#
#      popViewport()

  dev.off()

}









#######
#plot_df <- plot_dfs[[2]][30:32,]
#plot_df$ID <- rownames(plot_df)
#plot_df <- melt(plot_df)
#p <- ggplot(plot_df, aes(x=ID, y=value))
#p <- p + geom_bar(stat="identity", aes(fill=variable))
#p <- p + scale_fill_manual(values = col_pal)
##p <- p + scale_y_continuous(, limits = c(0,100))
#p <- p + scale_y_continuous(expand = c(0,0))
#p <- p + theme(
#  text = element_text(size=20), 
#  axis.ticks.x = element_blank(),
#  axis.text.x = element_blank()
#)
#p
#######


### 5. Create composition barplots ###

# convert Counts to a dataframe:
countsDF <- as.data.frame(t(do.call("cbind", Counts)))
# add rownames as column and melt dataframe:
countsDF$gene_type <- rownames(countsDF)

# calculate percentages as new data frame:
perCountsDF <- as.data.frame(apply(countsDF[,1:ncol(countsDF)-1], 2, function(x) {
  return(as.integer(x)/as.integer(sum(x))*100)
}))
perCountsDF$gene_type <- countsDF$gene_type
 

# create composition barplots of CountsDF and perCountsDF:
cDFs <- list(countsDF, perCountsDF)
Plots <- list()
for (i in 1:2) {
  pCounts <- melt(cDFs[i], variable.name = "sample")
  pCounts$gene_type <- factor(pCounts$gene_type, levels = c("non-coding", "ribosomal", "other", "repeat", "protein-coding"))
  
  # plot data as barplot:
  p <- ggplot(pCounts, aes(x=sample, y=value))
  p <- p + geom_bar(stat="identity", aes(fill=gene_type))
  p <- p + scale_fill_manual(values = c("#00BF7D", "#AA4499", "#00B0F6", "#DDCC77", "#F8766D"))
  Plots[[i]] <- p + theme(text = element_text(size=20), axis.text.x = element_text(angle = 90))
  i=i+1
}



pdf(file = paste0(plotDir, "compBarplotCounts.pdf"), height=20, width=35)
Plots[[1]]
dev.off()

pdf(file = paste0(plotDir, "compBarplotPercent.pdf"), height=20, width=35)
Plots[[2]]
dev.off()
  
  
save.image(file=paste0(RobjectDir, "/QCimage.RData"))







