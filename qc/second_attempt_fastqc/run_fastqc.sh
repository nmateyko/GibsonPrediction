#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --account=st-cdeboer-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=test
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out
#SBATCH --mail-user=nick.mateyko@ubc.ca
#SBATCH --mail-type=ALL

cd $SLURM_SUBMIT_DIR

conda activate misc

# Set base directory
data_dir="/arc/project/st-cdeboer-1/nick/GSC_seq_20240313"

# Declare associative array: sample_name => sample_index
declare -A samples=(
    [linear_rep_1]="TAGAGGCT-CGTGGAGT"
    [linear_rep_2]="CGATGGTG-CACTACTT"
    [input_rep_1]="CCGGAACA-TTCCACCG"
    [input_rep_2]="TTAGTGGT-CAGACCTT"
    [assembled_rep_1]="CACAATCG-TGAGTTGA"
    [assembled_rep_2]="TAATGCAG-AACGCCTT"
)

# Output directory for FastQC results
output_dir="fastqc_results"
mkdir -p "$output_dir"

# Load fastqc module if needed
# module load fastqc

for sample_name in "${!samples[@]}"; do
    sample_index="${samples[$sample_name]}"
    
    # File paths for R1 and R2
    r1="${data_dir}/PX3270_${sample_index}/150bp/PX3270_${sample_index}_1_150bp_3_lanes.merge.fastq.gz"
    r2="${data_dir}/PX3270_${sample_index}/150bp/PX3270_${sample_index}_2_150bp_3_lanes.merge.fastq.gz"

    # Run FastQC
    echo "Running FastQC on $sample_name"
    fastqc -o "$output_dir" "$r1" "$r2"

    # Rename FastQC output files
    mv "${output_dir}/PX3270_${sample_index}_1_150bp_3_lanes.merge_fastqc.html" "${output_dir}/${sample_name}_R1_fastqc.html"
    mv "${output_dir}/PX3270_${sample_index}_1_150bp_3_lanes.merge_fastqc.zip"  "${output_dir}/${sample_name}_R1_fastqc.zip"
    mv "${output_dir}/PX3270_${sample_index}_2_150bp_3_lanes.merge_fastqc.html" "${output_dir}/${sample_name}_R2_fastqc.html"
    mv "${output_dir}/PX3270_${sample_index}_2_150bp_3_lanes.merge_fastqc.zip"  "${output_dir}/${sample_name}_R2_fastqc.zip"
done
