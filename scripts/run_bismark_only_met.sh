#!/bin/bash
# get data from the adaptor-and-quality trimmed fastq files, and run bismark
# We run bismark single-end, non-directional mapping (which invovles both R1 and R2)


# Then, we deduplicate the bismark output files, merge all 2 files,
# extract the methylation information with bismark_methylation_extractor
# and finally generate the coverage report using coverage2cytosine


## This is needed to conform to the file name format in bismark. It will use the the base name in R1 to generate output file names
input_dir=$1 # the entire path to the directory containing R1 and R2
trim_fq1=$2 # just the file name for R1
trim_fq2=$3 # just the file name for R2
out_dir=$4
sample=$5 #sample name, or data_des. Can be arbitrary. Do not use to identify a file
recompute=$6
multicore=$7
genome_ref=$8


module load samtools

cov_dir=$out_dir/coverage
mkdir -p $cov_dir

bismark_out=$out_dir/bismark
mkdir -p $bismark_out

## below, we have anticipated the output file format from bismark. They seems crazy, but do not change them!!!
R1_base_name=`basename ${trim_fq1} .fq.gz` # bismark automatically use the base name from R1 to generate the output files
R2_base_name=`basename ${trim_fq2} .fq.gz` # bismark automatically use the base name from R1 to generate the output files
bismark_out_R1=$bismark_out/mapping_outcome_R1/${R1_base_name}_bismark_bt2 #single-end mapping output
bismark_out_R2=$bismark_out/mapping_outcome_R2/${R2_base_name}_bismark_bt2 #single-end mapping output


## single-end mapping for un-mapping R1
echo "--------------single-end mapping for un-mapping R1----------------------"
if [ -s  ${bismark_out_R1}.deduplicated.bam ] && [ -s  ${bismark_out_R1}.bam ] && [[ $recompute -eq 0 ]]
then
    echo "Use precomputed file: ${bismark_out_R1}.bam"
else
    bismark  --multicore $multicore   --non_directional   --unmapped                 \
        --phred33-quals                            \
        --output_dir $bismark_out/mapping_outcome_R1 --temp_dir $bismark_out/mapping_outcome_R1         \
        $genome_ref   ${input_dir}/$trim_fq1
fi


echo "---------single-end mapping for un-mapping R1: deduplication-----------"
if [ -s ${bismark_out_R2}.bam ] && [ -s ${bismark_out_R1}.deduplicated.bam ] && [[ $recompute -eq 0 ]]
then
    echo "Use precomputed file: ${bismark_out_R1}.deduplicated.bam"
else
    deduplicate_bismark -s --bam  --output_dir $bismark_out/mapping_outcome_R1 ${bismark_out_R1}.bam  # single-end
fi


## single-end mapping for un-mapping R2
echo "--------------single-end mapping for un-mapping R2----------------------"
if [ -s  ${bismark_out_R2}.deduplicated.bam ] && [ -s ${bismark_out_R2}.bam ] && [[ $recompute -eq 0 ]]
then
    echo "Use precomputed file: ${bismark_out_R2}.bam"
else
    bismark --multicore $multicore    --non_directional   --unmapped                 \
        --phred33-quals                          \
        --output_dir $bismark_out/mapping_outcome_R2 --temp_dir $bismark_out/mapping_outcome_R2         \
        $genome_ref   ${input_dir}/$trim_fq2
fi

echo "---------single-end mapping for un-mapping R2: deduplication-----------"
if [ -s $bismark_out/$sample.deduplicated.bam ] && [ -s ${bismark_out_R2}.deduplicated.bam ] && [[ $recompute -eq 0 ]]
then
    echo "Use precomputed file: ${bismark_out_R2}.deduplicated.bam"
else
    deduplicate_bismark -s --bam  --output_dir $bismark_out/mapping_outcome_R2 ${bismark_out_R2}.bam  # single-end
fi


echo "---------merge all bam files-----------"
if [ -s ${cov_dir}/${sample}.deduplicated.bismark.cov.gz  ] && [ -s $bismark_out/$sample.deduplicated.bam ] && [[ $recompute -eq 0 ]]
then
    echo "Use precomputed file: $bismark_out/$sample.deduplicated.bam"
else
    samtools merge -f $bismark_out/$sample.deduplicated.bam            \
                           ${bismark_out_R1}.deduplicated.bam                      \
                           ${bismark_out_R2}.deduplicated.bam
fi


echo "---------run bismark_methylation_extractor-----------"
if [  -f ${cov_dir}/coverage.done ] && [ -s ${cov_dir}/${sample}.deduplicated.bismark.cov.gz ]  && [[ $recompute -eq 0 ]]
then
    echo "Use precomputed file: ${cov_dir}/${sample}.deduplicated.bismark.cov.gz"
else
    bismark_methylation_extractor --multicore 4 -s --bedGraph --counts  --buffer_size 10G  --CX  --genome_folder $genome_ref -o $cov_dir $bismark_out/$sample.deduplicated.bam
fi

echo "coverage_file_generated" > ${cov_dir}/coverage.done
