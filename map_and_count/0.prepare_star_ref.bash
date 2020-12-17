#!/bin/bash

#make modules to load here

module load gi/samtools/1.2

#number of cores
numcores=6

#genome directories
species="L1MD3"
projectname="hgsoc_repeats/RNA-seq-final2"

homeDir="/share/ScratchGeneral/jamtor"
projectDir="$homeDir/projects/$projectname"
genomeDir="$projectDir/genome/"
genomeFile="$genomeDir/$species/star/$species.fa"
annotationFile="$genomeDir/$species.gtf"
outDir="$genomeDir/$species/star"

mkdir -p $outDir

#log directory
logDir="$projectDir/scripts/map_and_count/logs"
mkdir -p $logDir

echo This is the genomeDir:
echo $genomeDir
echo -e
echo This is the genomeFile:
echo $genomeFile
echo -e
echo This is the annotationFile:
echo $annotationFile
echo -e
echo This is the outDir
echo $outDir
echo -e
echo "This is the logDir:"
echo $logDir

#generate the star reference files:
star_ref_line="STAR --runMode genomeGenerate \
	--genomeDir $genomeDir \
	--genomeFastaFiles $genomeFile --runThreadN $numcores --outFileNamePrefix \
	$outDir/star_ref_"

echo -e
echo This is the star_ref_line:
echo $star_ref_line

#submit job with name 'RSEM_count_$sample' to 15 cluster cores:
qsub -N STAR_ref_$genomeName -wd $logDir -b y -j y -R y -pe smp $numcores -V $star_ref_line
