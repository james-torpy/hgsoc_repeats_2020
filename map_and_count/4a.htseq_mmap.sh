#!/bin/bash

source /home/jamtor/.bashrc

conda activate py2.7

module load gi/boost/1.53.0
module load gi/zlib/1.2.8
module load aledre/samtools/prebuilt/1.10

sample_name=$1
cores=$2
split=$3

#sample_name="AOCS-134-4"
#cores=2
#split="no"

echo "Sample name: $sample_name"
echo "Number cores: $cores"
echo "Split bam for parallelisation? $split"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/hgsoc_repeats/RNA-seq"
results_dir="$project_dir/results"
bam_dir="$results_dir/star/GC/$sample_name"
genome_dir="$project_dir/genome"
out_dir="$results_dir/htseq/$sample_name"
log_dir="$project_dir/logs/"

mkdir -p $log_dir
mkdir -p $out_dir

gff="$genome_dir/custom3.repeats.hg38.gtf"

echo This is the out_dir:
echo $out_dir

echo This is the gffFile:
echo $gff

bam_file="$bam_dir/sorted.by.name.mmappers.bam"
echo This is the bam_file:
echo $bam_file

if [ $split = "yes" ]; then

  # split bam for parallel processing:
  bamtools split -in $bam_file -reference
  
  # fetch filenames:
  filepaths=( $(ls -d $bam_dir/* | grep chr | grep Aligned | grep -v K | \
    grep -v GL | grep -v JH | grep -v chrUn | grep -v chrM) )
  
  filenames=( $(ls $bam_dir | grep chr | grep Aligned | grep -v K | \
    grep -v GL | grep -v JH | grep -v chrUn | grep -v chrM) )

  # remove original bam:
  if [ -f $bam_file ]; then
    rm $bam_file
  fi
  
  for s in ${filenames[@]}; do 
  	
    sname=$(echo $s | sed "s/^.*REF_//g" | sed "s/.bam//g")
  
    #perform analysis on bam files using the reference genome:
    htseq_line="htseq-count \
      -f bam \
      -i ID \
      -t sequence_feature \
      --stranded=no \
      -a 0 \
      $bam_dir/$s \
      $gff >> $out_dir$sname.custom3.htseq.txt"
    echo This is the htseq_line1:
    echo $htseq_line
    echo -e
  
    qsub -N $sname.htseq -b y -wd $log_dir -j y -R y -pe smp $cores -V $htseq_line
  
  done

else 

 	htseq_line="htseq-count \
      -f bam \
      -i ID \
      -t sequence_feature \
      --stranded=no \
      -a 0 \
      $bam_dir/sorted.by.name.mmappers.bam  \
      $gff >> $out_dir/counts.gc.custom3.htseq.txt"
    echo This is the htseq_line:
    echo $htseq_line
    echo -e

    htseq-count \
      -f bam \
      -i ID \
      -t sequence_feature \
      --stranded=no \
      -a 0 \
      $bam_dir/sorted.by.name.mmappers.bam  \
      $gff >> $out_dir/counts.gc.custom3.htseq.txt

#    # remove no_feature and non_unique reads from sam file:
#    grep -v "no_feature" $out_dir/counted.reads.sam | grep -v "not_unique" | grep -v "ambiguous" > $out_dir/counted.reads.final.sam
#    rm $out_dir/counted.reads.sam

  # remove in file:
  rm $bam_file

fi

