#!/bin/bash

# Default line length
line_length=60

# Help message
usage() {
  echo "Usage: $0 -s species_file -i orthogroup_dir -o output_dir"
  exit 1
}

# Parse options
while getopts "s:i:o:" opt; do
  case $opt in
    s) species_file="$OPTARG" ;;
    i) orthogroup_dir="$OPTARG" ;;
    o) output_dir="$OPTARG" ;;
    *) usage ;;
  esac
done

# Check required arguments
if [[ -z $species_file || -z $orthogroup_dir || -z $output_dir ]]; then
  usage
fi

# Begin script logic
for orthogroup in "$orthogroup_dir"/*.out.trim; do
    echo "Processing $orthogroup"
    
    existing_species=$(grep ">" "$orthogroup" | sed 's/>//; s/_.*//')
    readarray -t all_species < "$species_file"
    missing_species=()

    for species in "${all_species[@]}"; do
        if ! printf '%s\n' "${existing_species[@]}" | grep -q "^$species$"; then
            missing_species+=("$species")
        fi
    done

    max_length=0
    current_length=0
    while read -r line; do
        if [[ "$line" == \>* ]]; then
            if (( current_length > max_length )); then
                max_length=$current_length
            fi
            current_length=0
        else
            current_length=$((current_length + ${#line}))
        fi
    done < "$orthogroup"

    if (( current_length > max_length )); then
        max_length=$current_length
    fi

    output_file="$output_dir/$(basename "$orthogroup")"
    cp "$orthogroup" "$output_file"

    for species in "${missing_species[@]}"; do
        echo "Adding missing species: $species"
        echo ">${species}_00|" >> "$output_file"
        padding=$(printf '%*s' "$max_length" | tr ' ' '-')
        while [ -n "$padding" ]; do
            echo "${padding:0:$line_length}" >> "$output_file"
            padding="${padding:$line_length}"
        done
    done
done
