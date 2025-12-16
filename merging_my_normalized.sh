#!/bin/bash

#SBATCH --account=prudenlab
#SBATCH --partition=normal_q
#SBATCH --mem=20G
#SBATCH -t 8:00:00
#SBATCH -n 50
#SBATCH -N 1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=dceglio@vt.edu
# Output file name

output_file="ChezLiz_merged_args_new.txt"

cd /projects/ciwars/ChezLiz/diamond_output/Abundances

# Find and merge files

for file in *_abundance_new.csv; do
    # Skip the output file to avoid self-merging
    if [ "$file" != "$output_file" ]; then

        # Print filename as a comment in the merged file. nope i made it a column so that i can just treat it as a csv directly
        #echo "# $file" >> "$output_file"

        # Append contents to the merged file, skipping the header line
        #tail -n +2 "$file" >> "$output_file"

        #so i think i just make it 1 instead of 2 but we'll see ig
        tail -n +1 "$file" >> "$output_file"
    fi
done


echo "Merge complete. Merged file: $output_file"
