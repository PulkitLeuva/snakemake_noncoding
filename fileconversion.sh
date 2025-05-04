#!/bin/bash

# Usage: ./process_benchmarking_data.sh input_file output_folder

# Input parameters
input_file=$1
output_folder=$2

# Ensure the output folder exists
mkdir -p "$output_folder"

###### CSCAPE ######
# Process for CSCAPE
awk '!/^#/ && length($4)==1 && length($5)==1 && $4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/ {print $1, $2, $4, $5}' OFS=',' "$input_file" > "$output_folder/cscape_input.txt"
#awk -F'\t' 'BEGIN {OFS=","} NR > 1 { $3=""; $6=""; print $0 }' "$input_file" | sed 's/,,/,/g' | sed 's/^,//' | sed 's/,$//' > "$output_folder/cscape_input.txt"
echo "File processed and saved as $output_folder/cscape_input.txt"

###### CANDRA ######
# Extract only CHROM, POS, REF, and ALT columns for CANDRA
awk -F'\t' 'BEGIN {OFS="\t"} {print $1, $2, $4, $5}' "$input_file" > "$output_folder/candra_input.txt"
echo "File processed and saved as $output_folder/candra_input.txt"

###### funseq2 ######
# Process for funseq2 with VCF format
{
  echo "##fileformat=VCFv4.0"
  echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  awk -F'\t' 'BEGIN {OFS="\t"} NR > 1 {print $1, $2, $3, $4, $5, ".", ".", "."}' "$input_file"
} > "$output_folder/funseq2_input.txt"
echo "File processed and saved as $output_folder/funseq2_input.txt"

###### annovar ######
# Process for annovar with VCF format
{
  echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  awk -F'\t' 'BEGIN {OFS="\t"} NR > 1 {print $1, $2, $3, $4, $5, ".", ".", "."}' "$input_file"
} > "$output_folder/annovar_input.vcf"
echo "File processed and saved as $output_folder/annovar_input.vcf"

###### CADD ######
# Process for CADD with VCF format
{
  echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  awk -F'\t' 'BEGIN {OFS="\t"} NR > 1 {print $1, $2, $3, $4, $5, ".", ".", "."}' "$input_file"
} > "$output_folder/cadd_input.vcf"
echo "File processed and saved as $output_folder/cadd_input.vcf"

###### fathmm ######
# Process for fathmm
awk -F'\t' '!/^#/ { print $1, $2, $4, $5 }' OFS=',' "$input_file" > "$output_folder/fathmm_input.txt"
#awk -F'\t' 'BEGIN {OFS=","} { $3=""; $6=""; print $0 }' "$input_file" | sed 's/,,/,/g' | sed 's/^,//' | sed 's/,$//' > "$output_folder/fathmm_input.txt"
echo "File processed and saved as $output_folder/fathmm_input.txt"

###### sift ######
# Process for sift with VCF format
{
  echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  awk -F'\t' 'BEGIN {OFS="\t"} NR > 1 {print $1, $2, $3, $4, $5, ".", ".", "."}' "$input_file" | sort -k1,1 -k2,2n
} > "$output_folder/sift_input.vcf"
echo "File processed and saved as $output_folder/sift_input.vcf"
