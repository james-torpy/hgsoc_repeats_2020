#!/bin/bash

source /home/jamtor/.bashrc

conda activate py2.7

module load gi/boost/1.53.0
module load gi/zlib/1.2.8
module load aledre/samtools/prebuilt/1.10

sample_name=$1
cores=$2
type=$3

#sample_name="AOCS-172-31-sub"
#cores=2
#type="GC"

echo "Sample name: $sample_name"
echo "Number cores: $cores"
echo "RNA type to be counted? $type"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/hgsoc_repeats/RNA-seq-final"
results_dir="$project_dir/results"
bam_dir="$results_dir/star/GC/$sample_name"
genome_dir="$project_dir/genome"
out_dir="$results_dir/htseq/$sample_name"
log_dir="$project_dir/logs/"

mkdir -p $log_dir
mkdir -p $out_dir

if [ $type = "GC" ]; then
  gff="$genome_dir/gencode.v35.basic.annotation.gff3"
  bam_file="$bam_dir/sorted.by.name.bam"
  out_file="$out_dir/counts.gc.htseq.txt"
  feature_type="exon"
else
  gff="$genome_dir/repeats.hg38.gff"
  bam_file="$bam_dir/sorted.by.name.mmappers.bam"
  out_file="$out_dir/counts.repeats.htseq.txt"
  feature_type="sequence_feature"
fi

echo This is the bam file:
echo $bam_file

echo This is the gff file:
echo $gff

echo This is the out file:
echo $out_file

htseq_line="htseq-count \
  -f bam \
  -i ID \
  -t $feature_type \
  --stranded=no \
  -a 0 \
  $bam_file  \
  $gff >> $out_file"

echo This is the htseq_line:
echo $htseq_line
echo -e

htseq-count \
  -f bam \
  -i ID \
  -t $feature_type \
  --stranded=no \
  -a 0 \
  $bam_file  \
  $gff >> $out_file

# remove in file:
echo "Removing $bam_file"
rm $bam_file


