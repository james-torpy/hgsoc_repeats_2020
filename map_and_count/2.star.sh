#!/bin/bash

module load gi/zlib/1.2.8
module load phuluu/samtools/1.4

# fetch input variables:
sample_name=$1
cores=$2
type=$3

#sample_name="AOCS-172-31-sub"
#cores=10
#type="GC"

# define dirs:
home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/hgsoc_repeats/RNA-seq-final"
in_dir="$project_dir/raw_files/fastqs"

out_dir="$project_dir/results/star/$type/$sample_name/"
mkdir -p $out_dir

if [ $type = "GC" ]; then
  echo "Mapping to Gencode..."
  genome_dir="$project_dir/genome/GC/star"
  comp_file="$out_dir/star_complete"
elif [ $type = "ribo" ]; then
  echo "Mapping to ribosome..."
  genome_dir="$project_dir/genome/ribo/star"
  comp_file="$out_dir/star_ribo_complete"
elif [ $type = "L1MD3" ]; then
  echo "Mapping to L1MD3..."
  genome_dir="$project_dir/genome/L1MD3/star"
  comp_file="$out_dir/star_L1MD3_complete"
fi

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

if [ ! $type = "ribo" ]; then
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

else

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
    --outSAMtype BAM Unsorted"

  echo -e
  echo "This is the starline:"
  echo $starline
  echo -e
  
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
    --outSAMtype BAM Unsorted

fi

touch $comp_file
