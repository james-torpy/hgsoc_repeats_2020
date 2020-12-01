#!/bin/bash

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/hgsoc_repeats/RNA-seq"
out_dir="$project_dir/raw_files/fastqs"

gadi_dir="/g/data1a/ku3/jt3341"
in_dir="jt3341/projects/hgsoc_repeats/RNA-seq/raw_files/fastqs"
temp_dir="$gadi_dir/mdss/hgsoc_repeats/RNA-seq/raw_files/fastqs"

sample_name=$1

#sample_name="AOCS-076-2"

mkdir -p $out_dir
cd $out_dir

for i in 1 2; do 

  filename=$sample_name\_$i.fastq.gz

  echo "Fetching $filename from mdss..."
  ssh jt3341@gadi.nci.org.au mdss get $in_dir/$filename $temp_dir
  
  echo "Transferring $filename from gadi..."
  rsync -avPS jt3341@gadi.nci.org.au:$temp_dir/$filename $out_dir
  
  echo "Deleting $filename from gadi..."
  ssh jt3341@gadi.nci.org.au rm $temp_dir/$filename
  
done