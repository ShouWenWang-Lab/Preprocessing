Scripts for preprocessing both the scBS-seq data and the scRNA-seq data.
# scBS-seq preprocessing
In `scripts/run_bismark_only_met.sh`,  we get data from the adaptor-and-quality trimmed fastq files, and run bismark. We run bismark single-end, non-directional mapping (which invovles both R1 and R2). Then, we deduplicate the bismark output files, merge all 2 files,
extract the methylation information with bismark_methylation_extractor
 and finally generate the coverage report using coverage2cytosine

### variables
- input_dir=$1 # the entire path to the directory 
- containing R1 and R2
- trim_fq1=$2 # just the file name for R1
- trim_fq2=$3 # just the file name for R2
- out_dir=$4
- sample=$5 #sample name, or data_des. Can be arbitrary. Do not use to identify a file
- recompute=$6
- multicore=$7
- genome_ref=$8 # reference genome

# scRNA-seq preprocessing
In `scripts/RNA_scLimeCat.sh`, we preprocess the fastq files from scRNA-seq obtained in the scSTRT-seq protocol, and generate the cell-by-count matrix.

### Variables:
- out_dir=$1
- R1_file_0=$2
- R2_file_0=$3
- script_dir=$4
- cell_barcode=$5 # white_list cell barcode directory
- sample_name=$6
- recompute=$7
- protocol=$8
- reference=$9 # reference genome
