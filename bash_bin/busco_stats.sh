#!/bin/bash

# Base directory
base_dir="/scratch/alpine/beyo2625/arthur_truffles/busco_analysis"

# Print header
echo -e "Sample\tDatabase\tMeasure\tValue"

for sample_dir in "$base_dir"/*; do
    [ -d "$sample_dir" ] || continue

    sample=$(basename "$sample_dir")

    for db in fung asco; do
        summary_file="$sample_dir/$db/short_summary.specific.${db}i_odb12.${db}.txt"
        
        if [ -f "$summary_file" ]; then
            awk -v sample="$sample" -v db="$db" '
            /^[[:space:]]*[0-9]+[[:space:]]+Complete BUSCOs/ {
                print sample"\t"db"\tComplete\t"$1
            }
            /^[[:space:]]*[0-9]+[[:space:]]+Complete and single-copy BUSCOs/ {
                print sample"\t"db"\tComplete_Single_Copy\t"$1
            }
            /^[[:space:]]*[0-9]+[[:space:]]+Complete and duplicated BUSCOs/ {
                print sample"\t"db"\tComplete_Duplicated\t"$1
            }
            /^[[:space:]]*[0-9]+[[:space:]]+Fragmented BUSCOs/ {
                print sample"\t"db"\tFragmented\t"$1
            }
            /^[[:space:]]*[0-9]+[[:space:]]+Missing BUSCOs/ {
                print sample"\t"db"\tMissing\t"$1
            }
            /^[[:space:]]*[0-9]+[[:space:]]+Total BUSCO groups searched/ {
                print sample"\t"db"\tTotal\t"$1
            }
            ' "$summary_file"
        fi
    done
done
