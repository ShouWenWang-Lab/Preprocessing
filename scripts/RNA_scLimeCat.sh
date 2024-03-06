#!/bin/bash


## need to allocate at least 20G (since mapping would consume ~17G)
## the R1 and R2 in scLimeCat-seq is the opposite: R1 corresponds to R2 in 10X

out_dir=$1
R1_file_0=$2
R2_file_0=$3
script_dir=$4
cell_barcode=$5
sample_name=$6
recompute=$7
protocol=$8
reference=$9

clean_dir=${out_dir}/clean_data
STAR_dir=${out_dir}/STAR/$sample_name
count_dir=${out_dir}/count

mkdir -p ${clean_dir}
mkdir -p ${STAR_dir}
mkdir -p ${count_dir}






#sample_name=simple
# specify the number of cores to use
cores=6

genome=$reference/star_1
gtf=$reference/genes/genes.gtf


module load gcc perl
#module load star/2.5.2b
module load samtools

echo "Step 1: extract barcode and UMI"
# UMI_tools is supposed to be a python tool in this environment (snakemake)
# --read2-stdout : throw away simple.R2.fq results, and only keep simple.R1.fq results
if [ -f ${clean_dir}/${sample_name}.fq.gz ] && [[ $recompute -eq 0 ]]
then
    echo "Use precomputed file: ${clean_dir}/${sample_name}.fq.gz"
else
    if [ "${protocol}" == "10X" ]
    then
        echo "Use 10X protocol: 16bp BC + 12bp UMI"
        bc_pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNNN #16BC + 12UMI
        R1_file=$R1_file_0
        R2_file=$R2_file_0
        umi_tools extract --bc-pattern=$bc_pattern \
                      --stdin $R1_file \
                      --stdout ${clean_dir}/${sample_name}.fq.gz \
                      --read2-stdout \
                      --read2-in $R2_file
    else
        echo "Use scLimeCat protocol: 8bp BC + 8bp UMI; change file R1 <-> R2"
        bc_pattern=CCCCCCCCNNNNNNNN # scLimeCat
        R2_file=$R1_file_0
        R1_file=$R2_file_0
        umi_tools extract --bc-pattern=$bc_pattern \
                  --stdin $R1_file \
                  --stdout ${clean_dir}/${sample_name}.fq.gz \
                  --read2-stdout \
                  --read2-in $R2_file \
                  --filter-cell-barcode \
                  --whitelist=$cell_barcode
    fi

fi


# We trim TSO and polyA to ensure that these artifical seqeunces do not interfere with alignment to genome
echo "Step 2: Trim TSO and polyA, discard low-quality reads"
if [ -f ${clean_dir}/${sample_name}.clean.fq.gz ]  && [[ $recompute -eq 0 ]]
then
    echo "Use precomputed file: ${clean_dir}/${sample_name}.clean.fq.gz"
else
    perl ${script_dir}/fastq_preprocessing/trim_TSO_polyA.pl ${clean_dir}/${sample_name}.fq.gz  ${clean_dir}/${sample_name}_trim.fq.gz 0
    # Trim low-quality bases from both ends using the Phred algorithm (seqtk is pre-installed)
    # see tutorial here https://github.com/lh3/seqtk
    seqtk trimfq ${clean_dir}/${sample_name}_trim.fq.gz | gzip - > ${clean_dir}/${sample_name}.clean.fq.gz
fi


echo "Step 3: run STAR -- Alignment"
# to run this, the $genome should be a star-generated reference
if [ -f ${STAR_dir}/${sample_name}.Aligned.sortedByCoord.out.bam ]  && [[ $recompute -eq 0 ]]
then
    echo "Use precomputed file: ${STAR_dir}/${sample_name}.Aligned.sortedByCoord.out.bam"
else
    STAR --runThreadN $cores \
         --genomeDir $genome \
         --readFilesIn ${clean_dir}/${sample_name}.clean.fq.gz \
         --readFilesCommand zcat \
         --outFilterMultimapNmax 1 \
         --outFileNamePrefix ${STAR_dir}/${sample_name}. \
         --outSAMtype BAM SortedByCoordinate

fi

echo "Step 4: add features"
if [ -f ${STAR_dir}/${sample_name}.Aligned.sortedByCoord.out.bam.featureCounts.bam ] && [[ $recompute -eq 0 ]]
then
    echo "Use precomputed file: ${STAR_dir}/${sample_name}.Aligned.sortedByCoord.out.bam.featureCounts.bam"
else
    featureCounts -a $gtf -o ${STAR_dir}/gene_assigned -R BAM ${STAR_dir}/${sample_name}.Aligned.sortedByCoord.out.bam -T $cores -g gene_name
fi

echo "Step 5: sort and index bam file"
if [ -f ${STAR_dir}/${sample_name}.assigned_sorted.bam ] && [[ $recompute -eq 0 ]]
then
    echo "Use precomputed file: ${STAR_dir}/${sample_name}.assigned_sorted.bam"
else
    samtools sort -m 15G ${STAR_dir}/${sample_name}.Aligned.sortedByCoord.out.bam.featureCounts.bam  -O BAM -o ${STAR_dir}/${sample_name}.assigned_sorted.bam
    samtools index ${STAR_dir}/${sample_name}.assigned_sorted.bam
fi

echo "Step 6: get UMI count table"
# if we force snakemake to rerun, it will first delete this file
if [ -f ${count_dir}/${sample_name}.UMI_counts.tsv ] && [[ $recompute -eq 0 ]]
then
    echo "Use precomputed file: ${count_dir}/${sample_name}.UMI_counts.tsv"
else
    umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell --wide-format-cell-counts -I ${STAR_dir}/${sample_name}.assigned_sorted.bam -S ${count_dir}/${sample_name}.UMI_counts.tsv
fi
