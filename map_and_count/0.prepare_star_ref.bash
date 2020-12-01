#!/bin/bash

#make modules to load here

module load gi/samtools/1.2

#number of cores
numcores=6

#genome directories
genomeName="GRCh38.primary_assembly.genome"
annotationName="gencode.v35.annotation"
projectname="hgsoc_repeats/RNA-seq"

homeDir="/share/ScratchGeneral/jamtor"
projectDir="$homeDir/projects/$projectname"
genomeDir="$projectDir/genome/"
genomeFile="$genomeDir/star/$genomeName.fa"
annotationFile="$genomeDir/$annotationName.gtf"
outDir="$genomeDir/star"

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
	--sjdbGTFfile $annotationFile --sjdbOverhang 74 --genomeDir $genomeDir \
	--genomeFastaFiles $genomeFile --runThreadN $numcores --outFileNamePrefix \
	$outDir/star_ref_"

echo -e
echo This is the star_ref_line:
echo $star_ref_line

#submit job with name 'RSEM_count_$sample' to 15 cluster cores:
qsub -N STAR_ref_$genomeName -wd $logDir -b y -j y -R y -pe smp $numcores -V $star_ref_line
