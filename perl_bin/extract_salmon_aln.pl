#!/usr/bin/perl
use strict;
use warnings;
use File::Find;

# Check for proper arguments
die "Usage: $0 <base_directory> <output_file>\n" unless @ARGV == 2;

my ($base_dir, $output_file) = @ARGV;

open my $out, '>', $output_file or die "Cannot open output file '$output_file': $!\n";
print $out "Sample\tMapping_Percentage\n";

find(\&process_log, $base_dir);

sub process_log {
    return unless /salmon_quant\.log$/;

    my $log_file = $File::Find::name;
    my $sample;

    # Extract sample name from path (e.g., SRR24633217_loca)
    if ($log_file =~ m|/([^/]+?)/logs/salmon_quant\.log$|) {
        $sample = $1;
    } else {
        warn "Could not parse sample name from $log_file\n";
        return;
    }

    open my $fh, '<', $log_file or do {
        warn "Could not open $log_file: $!";
        return;
    };

    my $mapping;
    while (<$fh>) {
        if (/Mapping rate\s*=\s*([\d.]+)%/) {
            $mapping = $1;
            last;
        }
    }
    close $fh;

    if (defined $mapping) {
        print $out "$sample\t$mapping%\n";
    } else {
        warn "No mapping rate found in $log_file\n";
    }
}

close $out;
