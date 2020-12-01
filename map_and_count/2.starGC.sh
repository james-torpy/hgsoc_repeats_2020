#!/bin/bash

date

module load gi/zlib/1.2.8
module load phuluu/samtools/1.4

# fetch input variables:
sample_name=$1
cores=$2

#sample_name="AOCS-152-2"
#cores=10

# define dirs:
home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/hgsoc_repeats/RNA-seq"
genome_dir="$project_dir/genome/star"
in_dir="$project_dir/raw_files/fastqs"
out_dir="$project_dir/results/star/GC/$sample_name/"
mkdir -p $out_dir

# fetch sample_name and create out_dir:
fq1=$in_dir/$sample_name\_1.fastq.gz
fq2=$in_dir/$sample_name\_2.fastq.gz

echo "This is the sample_name:"
echo $sample_name
echo "This is the fq1:"
echo $fq1
echo "This is the fq2:"
echo $fq2
echo "This is the genome_dir:"
echo $genome_dir
echo -e
echo "This is the out_dir:"
echo $out_dir

# align reads of input files with STAR, output into .bam files:
starline="/home/jamtor/local/bin/STAR --runMode alignReads \
      --readFilesCommand zcat \
      --genomeDir $genome_dir \
    --outFilterType BySJout \
    --outSAMattributes NH HI AS NM MD\
    --outFilterMultimapNmax 999 \
    --outMultimapperOrder Random \
    --runRNGseed 666 \
    --outSAMmultNmax 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1500000 \
    --alignMatesGapMax 1500000 \
    --alignSJoverhangMin 6 \
    --alignSJDBoverhangMin 1 \
    --readFilesIn $fq1 $fq2 \
    --outFileNamePrefix $out_dir \
    --runThreadN $cores \
    --outFilterMatchNmin 76 \
  --chimSegmentMin 25 \
    --chimJunctionOverhangMin 25 \
    --chimScoreMin 0 \
    --chimScoreDropMax 20 \
    --chimScoreSeparation 10 \
    --chimScoreJunctionNonGTAG -1 \
    --outSAMtype BAM Unsorted"

echo -e
echo "This is the starline:"
echo $starline

# run command:
/home/jamtor/local/bin/STAR --runMode alignReads \
      --readFilesCommand zcat \
      --genomeDir $genome_dir \
    --outFilterType BySJout \
    --outSAMattributes NH HI AS NM MD\
    --outFilterMultimapNmax 999 \
    --outMultimapperOrder Random \
    --runRNGseed 666 \
    --outSAMmultNmax 1 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1500000 \
    --alignMatesGapMax 1500000 \
    --alignSJoverhangMin 6 \
    --alignSJDBoverhangMin 1 \
    --readFilesIn $fq1 $fq2 \
    --outFileNamePrefix $out_dir \
    --runThreadN $cores \
    --outFilterMatchNmin 76 \
  --chimSegmentMin 25 \
    --chimJunctionOverhangMin 25 \
    --chimScoreMin 0 \
    --chimScoreDropMax 20 \
    --chimScoreSeparation 10 \
    --chimScoreJunctionNonGTAG -1 \
    --outSAMtype BAM Unsorted

rm $fq1
rm $fq2

date