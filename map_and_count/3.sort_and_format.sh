
module load aledre/samtools/prebuilt/1.10

sample_name=$1

#sample_name="AOCS-108-2"

echo "Sample name: $sample_name"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/hgsoc_repeats/RNA-seq"
results_dir="$project_dir/results"
bam_dir="$results_dir/star/GC/$sample_name"
log_dir="$project_dir/logs/"

mkdir -p $log_dir

bam_file="$bam_dir/Aligned.out.bam"
echo This is the bam_file:
echo $bam_file

# sort by name and change all NH:i:>1 to 1:
echo "Sorting by name, adjusting NH flags..."
samtools view -h $bam_file | \
  samtools sort -n -@ 9 -T "$bam_dir/$sample_name" -O sam | \
  sed "s/NH:i:[0-9]*/NH:i:1/g" | \
  samtools view -b > "$bam_dir/sorted.by.name.mmappers.bam"
samtools quickcheck -v $bam_dir/*.bam > $bam_dir/bad_bams.fofn && echo 'all ok' || 
  echo 'some files failed check, see bad_bams.fofn'

## check NH values - should only be NH:i:1:
#samtools view "$bam_dir/sorted.by.name.mmappers.bam" | awk '{print $12}' | \
#  sort | uniq

# remove original bam, then add dummy for snakemake:
#if [ -f $bam_file ]; then
#  rm $bam_file
#  touch $bam_file
#fi