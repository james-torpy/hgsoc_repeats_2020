#!/bin/bash

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/hgsoc_repeats/RNA-seq"
out_dir="$project_dir/results/star/GC"

gadi_dir="/g/data1a/ku3/jt3341"
in_dir="jt3341/projects/hgsoc_repeats/RNA-seq/results/star/GC"
temp_dir="$gadi_dir/mdss/hgsoc_repeats/RNA-seq/results/star/GC"

#filename=$1
#sample_name=$(echo $filename | sed "s/\.tar.gz//")
#sample_name=$1
#sample_name="AOCS-141-10"
filename=$sample_name.tar.gz

mkdir -p $out_dir
cd $out_dir

echo "Fetching $filename from mdss..."
ssh jt3341@gadi.nci.org.au mdss get $in_dir/$filename $temp_dir

echo "Transferring $filename from gadi..."
rsync -avPS jt3341@gadi.nci.org.au:$temp_dir/$filename $out_dir

echo "Deleting $filename from gadi..."
ssh jt3341@gadi.nci.org.au rm $temp_dir/$filename

echo "Untarring $filename..."
tar -zxvf $out_dir/$filename

echo "Deleting tar file..."
rm $out_dir/$filename

mkdir -p $out_dir/$sample_name/temp/

if [ -f "$out_dir/$sample_name/Aligned.novosortedByName.out.bam" ]; then
  echo "Removing all except name sorted bam..."
  mv "$out_dir/$sample_name/Aligned.novosortedByName.out.bam" "$out_dir/$sample_name/temp/sorted.by.name.bam"
  rm $out_dir/$sample_name/*
  mv "$out_dir/$sample_name/temp/sorted.by.name.bam" "$out_dir/$sample_name"
  rm -r "$out_dir/$sample_name/temp"
  # create dummy file for snakemake:
  touch "$out_dir/$sample_name/Aligned.out.bam"
else
  echo "Removing all except non-sorted bam..."
  mv "$out_dir/$sample_name/Aligned.out.bam" "$out_dir/$sample_name/temp/"
  rm $out_dir/$sample_name/*
  mv "$out_dir/$sample_name/temp/Aligned.out.bam" "$out_dir/$sample_name"
  rm -r "$out_dir/$sample_name/temp"
fi

