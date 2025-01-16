#!/bin/bash

#!/bin/bash

# Input file
input_file=$1

# Output file
output_file="exons.sorted.bed"

# Step 1: Extract specific columns (3, 10, 11, 2, 12, 4) using awk and save to the output file
awk 'BEGIN {FS="\t"; OFS="\t"} {print $3, $10, $11, $2, $12, $4}' "$input_file" > "$output_file"

# Step 2: Process the extracted columns to split the exonStarts and exonEnds into individual exons
awk '
BEGIN {FS="\t"; OFS="\t"}
{
    # Split exonStarts and exonEnds into arrays
    n = split($2, starts, ","); # Split exonStarts into array "starts"
    split($3, ends, ",");       # Split exonEnds into array "ends"

    # Iterate over the arrays and create a new line for each exon
    for (i = 1; i <= n; i++) {
        # Only print rows with valid exonStarts and exonEnds
        if (starts[i] != "" && ends[i] != "") {
            print $1, starts[i], ends[i];
        }
    }
}
' "$output_file" > "$output_file.tmp"

# Step 3: Remove header if it exists
# Check if the first line starts with "chrom" (or any other header)
if head -n 1 "$output_file.tmp" | grep -q "^chrom"; then
    grep -v "^chrom" "$output_file.tmp" > "$output_file.tmpi"
else
    mv "$output_file.tmp" "$output_file.tmpi"
fi

# Step 4: Sort the BED file by chromosome and start position
sort -k1,1 -k2,2n "$output_file.tmpi" > "$output_file"
output_file="${output_file%.tmp}"

# Clean up temporary files
rm -f "$output_file.tmpi"
rm "$output_file.tmp"

# Step 5: Filter chromosomes to include only standard chromosomes (1-22, X, Y)
filtered_file="filtered_$output_file"
awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ {print $0}' "$output_file" > "$filtered_file"
rm $output_file
# Step 6: Compress and index the filtered file
final_file="targets.bed"
mv "$filtered_file" "$final_file"
bgzip -k "$final_file" > "$final_file.gz"
tabix -p bed "$final_file.gz"

# Print completion message
echo "Conversion complete. Output saved to $final_file"

