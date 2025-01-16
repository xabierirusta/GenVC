#!/bin/bash

base_dir=$(pwd)
output="${base_dir}/results"

mkdir -p $output

cd $output
# Define the main directories
directories=("alignment" "trimmed" "post_alignment" "plots" "variant_call" "annotated")

# Define the subdirectories for variant_call
sub_directories=("delly" "cnvkit" "haplotypecaller" "mutect2")

# Loop through main directories
for dir in "${directories[@]}"; do
    mkdir -p $dir
    
    # Check if the directory is 'variant_call' to create subdirectories
    if [ "$dir" == "variant_call" ]; then
        cd $dir  # Change into 'variant_call' directory

        # Loop through subdirectories and create them
        for subdir in "${sub_directories[@]}"; do
            mkdir -p $subdir
        done

        cd ..  # Return to the parent directory after creating subdirectories
    # Handle the 'plots' directory to create 'cnvkit' subdirectory
    elif [ "$dir" == "plots" ]; then
        cd $dir
        mkdir -p "cnvkit"
        cd ..  # Return to the parent directory after creating 'cnvkit'
    fi
done
