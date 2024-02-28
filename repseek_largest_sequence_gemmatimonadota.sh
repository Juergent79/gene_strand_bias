#!/bin/bash

# Define the input folder containing the *.fna files
input_folder="/home/aap/mnt/Juergen/Projects/Gemmatimonas_AP64/genome_structure/genome_files/Gemmatimonadota"

for file in "$input_folder"/*/*[0-9]_genomic.fna; do
    # Get the filename without extension
    filename=$(basename "$file" .fna)
    folder=$(dirname "$file")


    # Use awk to find the largest sequence in the FASTA file
    largest_sequence=$(awk '/^>/{if (NR>1) print header ORS seq; seq=""; header=$0;} {if (NR>1) seq = seq $0} END{print header ORS seq}' "$file" | \
                      awk -v RS=">" '{if (length($0) > max) {max = length($0); seq = $0}} END {print seq}' | \
                      sed -n '2,$p')
   

    # Create a temporary FASTA file for the largest sequence
    temp_fasta="/tmp/largest_sequence_$filename.fst"
    echo ">largest_sequence" > "$temp_fasta"
    echo "$largest_sequence" >> "$temp_fasta"

    # Run repseek on the temporary FASTA file and save the result in a .rep file
    /home/aap/scripts/repseek -l 32 "$temp_fasta" > "$folder/$filename.rep"

    # Remove the temporary FASTA file
    rm "$temp_fasta"

    echo "Processed $file and saved results to $filename.rep"
done

