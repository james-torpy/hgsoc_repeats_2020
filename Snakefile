# Run command:
# snakemake --reason --cores 160 --cluster-config cluster.json --cluster 'qsub -pe smp {cluster.cores} -N rep_quant.smk -wd '/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/logs' -b y -j y -V -P TumourProgression' -j 80

# DAG command:
# snakemake --dag | dot -Tsvg > dag.svg

### This script counts repeats, protein coding and ribosomal read pairs, 
# then performs a range of DE ###

import glob
import os

# define directories:
R_dir = '/share/ClusterShare/software/contrib/briglo/octR/src/R-3.6.0/builddir/bin/'
in_dir = './results/star/GC/'

# fetch sample names:
# from fastqs:
SAMPLES = list([

#    #1 - done
#    "AOCS-034-2", "AOCS-034-4", "AOCS-055-2", "AOCS-056-2", 
#    "AOCS-057-2", "AOCS-058-2", "AOCS-059-2", "AOCS-060-2"

    #2 - screen 7907.pts-15.dice01
#    "AOCS-061-2", "AOCS-063-2", "AOCS-064-2", "AOCS-064-4",
#    "AOCS-065-2", "AOCS-065-4", "AOCS-075-2", "AOCS-078-2"

    #3 - screen 8412.pts-12.dice01
#    "AOCS-081-2", "AOCS-083-2", "AOCS-085-2", "AOCS-086-2", 
#    "AOCS-086-4", "AOCS-088-2", "AOCS-088-4", "AOCS-090-2"

    #4 - screen 18942.pts-23.dice01
#    "AOCS-172-31", "AOCS-173-31", "AOCS-174-31", "AOCS-175-31", 
#    "AOCS-176-31", "AOCS-177-31", "AOCS-178-31", "AOCS-167-29"

    #5 - screen 20231.pts-23.dice01
#    "AOCS-114-2", "AOCS-115-2", "AOCS-116-2", "AOCS-117-4", 
#    "AOCS-119-4", "AOCS-120-4", "AOCS-122-2", "AOCS-123-2"

    #1 - screen 28737.pts-9.dice01
    "AOCS-105-2", "AOCS-106-2", "AOCS-107-2", "AOCS-108-2", 
    "AOCS-109-2", "AOCS-111-2", "AOCS-112-2", "AOCS-113-2"

    #3
#    "AOCS-124-2", "AOCS-126-2", "AOCS-128-2", "AOCS-131-2", 
#    "AOCS-132-2", "AOCS-133-2", "AOCS-135-4", "AOCS-135-10"

    #4
#    "AOCS-138-4", "AOCS-139-2", "AOCS-139-7", "AOCS-139-24", 
#    "AOCS-139-26", "AOCS-141-10", "AOCS-143-2", "AOCS-144-2"

    #1
#    "AOCS-145-2", "AOCS-146-2", "AOCS-147-2", "AOCS-148-2", 
#    "AOCS-149-2", "AOCS-150-4","AOCS-152-2", "AOCS-155-4"

    #2
#    "AOCS-157-2", "AOCS-159-2", "AOCS-160-2", "AOCS-163-2", 
#    "AOCS-164-2", "AOCS-167-27", "AOCS-170-2", "AOCS-170-4"

    #3
#    "AOCS-091-2", "AOCS-092-2", "AOCS-092-4", "AOCS-093-2", 
#    "AOCS-093-4", "AOCS-095-2", "AOCS-095-4", "AOCS-097-2"

    #4
#    "AOCS-125-2", "AOCS-130-2", "AOCS-137-10", "AOCS-091-4", 
#    "AOCS-104-2", "AOCS-064-4", "AOCS-093-8", "AOCS-150-10"

    #1
#    "AOCS-079-2", "AOCS-076-2", "AOCS-096-2", "AOCS-134-4", 
#    "AOCS-142-4", "AOCS-167-4", "AOCS-137-4", "AOCS-141-4"

])

# counting parameters:
split = 'no'

#rule all:
#    input:
#        expand(
#            'raw_files/fastqs/{sample}_1.fastq.gz',
#            sample=SAMPLES
#        ),
#        expand(
#            'raw_files/fastqs/{sample}_2.fastq.gz',
#            sample=SAMPLES
#        )

#rule all:
#    input:
#        expand(
#            'results/star/GC/{sample}/Aligned.out.bam',
#            sample=SAMPLES
#        )

#rule all:
#    input:
#        expand(
#            'results/star/GC/{sample}/sorted.by.name.mmappers.bam',
#            sample=SAMPLES
#        )

rule all:
    input:
        expand(
            'results/htseq/{sample}/counts.gc.custom3.htseq.txt',
            sample=SAMPLES
        )

rule transfer:
    output:
        fq1='raw_files/fastqs/{sample}_1.fastq.gz',
        fq2='raw_files/fastqs/{sample}_2.fastq.gz'
    threads: 2
    shell:
        'cd logs; ' + 
        ' ../scripts/map_and_count/1b.transfer_fq.sh' + 
        ' {wildcards.sample}' +
        ' 2> {wildcards.sample}.mdss.transfer.errors'

rule star:
    input:
        fq1='raw_files/fastqs/{sample}_1.fastq.gz',
        fq2='raw_files/fastqs/{sample}_2.fastq.gz'
    output:
        'results/star/GC/{sample}/Aligned.out.bam'
    threads: 10
    shell:
        'cd logs; ' + 
        ' ../scripts/map_and_count/2.starGC.sh' + 
        ' {wildcards.sample}' +
        ' 10' + # cores
        ' {split} 2> {wildcards.sample}.mdss.star.errors'

rule sort:
    input:
        'results/star/GC/{sample}/Aligned.out.bam'
    output:
        'results/star/GC/{sample}/sorted.by.name.mmappers.bam'
    threads: 10
    shell:
        'cd logs; ' + 
        ' ../scripts/map_and_count/3.sort_and_format.sh' + 
        ' {wildcards.sample} 2> {wildcards.sample}.mdss.sort.errors'

rule htseq:
    input:
        'results/star/GC/{sample}/sorted.by.name.mmappers.bam'
    output:
        'results/htseq/{sample}/counts.gc.custom3.htseq.txt'
    threads: 2
    shell:
        'cd logs; ' + 
        ' ../scripts/map_and_count/4a.htseq_mmap.sh' + 
        ' {wildcards.sample}' +
        ' 2' + # cores
        ' {split} 2> {wildcards.sample}.mdss.fqs.htseq.errors'


