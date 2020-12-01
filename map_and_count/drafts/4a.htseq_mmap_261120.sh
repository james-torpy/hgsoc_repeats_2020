#!/bin/bash

source /home/jamtor/.bashrc

conda activate py2.7

module load gi/boost/1.53.0
module load gi/zlib/1.2.8
module load aledre/samtools/prebuilt/1.10

sample_name=$1
cores=$2
split=$3

#sample_name="AOCS-088-2"
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

gff="$genome_dir/gc.v35.genes.custom3.repeats.gff3"

echo This is the out_dir:
echo $out_dir

echo This is the gffFile:
echo $gff

if [ -f "$bam_dir/sorted.by.name.bam" ]; then

  bam_file="$bam_dir/sorted.by.name.bam"
  echo This is the bam_file:
  echo $bam_file

  # change all NH:i:>1 to 1:
  \time samtools view -h $bam_file | \
    sed "s/NH:i:[0-9]*/NH:i:1/g" | \
    samtools view -b > "$bam_dir/sorted.by.name.mmappers.bam"

  samtools quickcheck -v $bam_dir/*.bam > $bam_dir/bad_bams.fofn && echo 'all ok' || 
    echo 'some files failed check, see bad_bams.fofn'

else

  bam_file="$bam_dir/Aligned.out.bam"
  echo This is the bam_file:
  echo $bam_file

  # sort by name and change all NH:i:>1 to 1:
  echo "Sorting by name, adjusting NH flags..."
  \time samtools view -h $bam_file | \
    samtools sort -n -T "$bam_dir/$sample_name" -O sam | \
    sed "s/NH:i:[0-9]*/NH:i:1/g" | \
    samtools view -b > "$bam_dir/sorted.by.name.mmappers.bam"

  samtools quickcheck -v $bam_dir/*.bam > $bam_dir/bad_bams.fofn && echo 'all ok' || 
    echo 'some files failed check, see bad_bams.fofn'

fi

## check NH values - should only be NH:i:1:
#samtools view "$bam_dir/sorted.by.name.mmappers.bam" | awk '{print $12}' | \
#  sort | uniq

# remove original bam, then add dummy for snakemake:
#if [ -f $bam_file ]; then
#  rm $bam_file
#  touch $bam_file
#fi

if [ $split = "yes" ]; then

  # split bam for parallel processing:
  bamtools split -in $bam_file -reference
  
  # fetch filenames:
  filepaths=( $(ls -d $bam_dir/* | grep chr | grep Aligned | grep -v K | \
    grep -v GL | grep -v JH | grep -v chrUn | grep -v chrM) )
  
  filenames=( $(ls $bam_dir | grep chr | grep Aligned | grep -v K | \
    grep -v GL | grep -v JH | grep -v chrUn | grep -v chrM) )

  # remove original bam:
  if [ -f "$bam_dir/sorted.by.name.mmappers.bam" ]; then
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
      -o $out_dir$sname.counted.reads.sam \
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
      -o $out_dir/counted.reads.sam \
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
      -o $out_dir/counted.reads.sam \
      $bam_dir/sorted.by.name.mmappers.bam  \
      $gff >> $out_dir/counts.gc.custom3.htseq.txt

    # remove no_feature and non_unique reads from sam file:
    grep -v "no_feature" $out_dir/counted.reads.sam | grep -v "not_unique" | grep -v "ambiguous" > $out_dir/counted.reads.final.sam
  	rm $out_dir/counted.reads.sam

fi

