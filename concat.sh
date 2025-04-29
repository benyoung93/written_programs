#!/usr/bin/perl
use strict;
use warnings;

# Check command-line arguments
my ($expected_file, $output_file, @ortholog_files) = @ARGV;
die("Usage: $0 <expected_file> <output_file> <ortholog_files...>\n") unless @ortholog_files;

# Read expected species identifiers
my %expected;
open(my $efh, '<', $expected_file) or die "Could not open '$expected_file': $!";
while (<$efh>) {
    chomp;
    s/^>//;  # Remove '>'
    $expected{$_} = "";  # Initialize with empty string
}
close($efh);

# Concatenate sequences for each species
foreach my $ortholog_file (@ortholog_files) {
    open(my $ofh, '<', $ortholog_file) or die "Could not open '$ortholog_file': $!";
    my $current_species = '';
    my $current_sequence = '';

    while (<$ofh>) {
        chomp;
        if (/^>(\S+)_(\d+)\|/) {
            # If we encounter a new header, save the current species sequence
            if ($current_species && exists $expected{$current_species}) {
                $expected{$current_species} .= $current_sequence;  # Append sequence
            }

            $current_species = $1;  # Get species name (spp)
            $current_sequence = '';  # Reset for new species
        } elsif ($current_species) {
            $current_sequence .= $_;  # Build sequence
        }
    }

    # Save the last species sequence
    if ($current_species && exists $expected{$current_species}) {
        $expected{$current_species} .= $current_sequence;
    }

    close($ofh);
}

# Write concatenated sequences and their lengths to output file
open(my $outfh, '>', $output_file) or die "Could not open '$output_file': $!";
foreach my $species (keys %expected) {
    my $sequence = $expected{$species};
    my $length = length($sequence);  # Get length of the concatenated sequence
    print $outfh ">$species\n$sequence\n";  # Write species and its concatenated sequence
    print "$species: $length bp\n";  # Print length to console
}
close($outfh);

print "Concatenated FASTA file created: $output_file\n";
