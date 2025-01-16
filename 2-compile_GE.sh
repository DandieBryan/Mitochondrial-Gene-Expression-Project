#!/bin/bash
#SBATCH --time=96:00:00
#SBATCH --mail-type=ALL



##########################################################################################
### MADE BY BRYAN LE, JULY 2024

# This part of the pipeline takes gene expression information from all samples and
# compiles it into a single file



###########################################################################################
# Define the directory containing the output files
output_directory="INSERT PATH TO OUTPUT DIRECTORY"
output_file="rnaseQC_1kg_compiled_GE_tpm.tab"

# Get a list of all sample files
sample_files=($(ls "$output_directory"))

# Create a temporary directory to hold intermediate files
temp_dir=$(mktemp -d)

# Declare an associative array to group files by sample ID
declare -A sample_files_map

# Group files by the first 7 characters of their name (sample ID)
for sample_file in "${sample_files[@]}"; do
    sample_id=$(basename "$sample_file" .txt | cut -c1-7)
    sample_files_map["$sample_id"]+="$output_directory/$sample_file "
done

# Function to average TPM values and sort by gene ID, rounding to six decimal places
average_tpm() {
    awk 'NR > 2 { # Skip the first two lines (GCT metadata)
        if ($3 ~ /^[0-9.]+$/) {
            sum[$1] += $3;
            count[$1] += 1;
        }
    } END {
        for (gene in sum) {
            avg = sum[gene] / count[gene];
            printf "%s\t%.6f\n", gene, avg; # Round to six decimal places
        }
    }' "$@" | sort -k1,1
}

# Process each sample group
for sample_id in "${!sample_files_map[@]}"; do
    # Get the list of files for this sample
    files=(${sample_files_map[$sample_id]})
    
    # Average the TPM values for this sample, sort by gene ID, and save to a temporary file
    average_tpm "${files[@]}" > "$temp_dir/${sample_id}_TPM.txt"
done

# Initialize the header with "SampleID"
header="GeneID"

# Sort the sample IDs such that 'HG' samples come first, followed by 'GM' samples
sorted_sample_ids=($(printf "%s\n" "${!sample_files_map[@]}" | grep '^HG' | sort -n))
sorted_sample_ids+=($(printf "%s\n" "${!sample_files_map[@]}" | grep '^GM' | sort -n))

for sample_id in "${sorted_sample_ids[@]}"; do
    header="$header\t$sample_id"
done

# Initialize the combined data with the first sample's gene IDs and averaged TPMs
cp "$temp_dir/${sorted_sample_ids[0]}_TPM.txt" "$temp_dir/combined_TPM.txt"

# Add the TPM columns from the other samples
for sample_id in "${sorted_sample_ids[@]:1}"; do
    awk '{print $2}' "$temp_dir/${sample_id}_TPM.txt" | paste "$temp_dir/combined_TPM.txt" - > "${temp_dir}/combined_TPM_new.txt"
    mv "${temp_dir}/combined_TPM_new.txt" "$temp_dir/combined_TPM.txt"
done

# Add the header to the final output file
echo -e "$header" > "$output_directory/$output_file"

# Add the combined data to the final output file
cat "$temp_dir/combined_TPM.txt" >> "$output_directory/$output_file"

# Clean up temporary files
rm -r "$temp_dir"

# Extract reads mapping to mitochondrial genes
sbatch 3-MT_GE_comparison.sh