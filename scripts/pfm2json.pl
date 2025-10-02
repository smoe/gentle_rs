#!/usr/bin/perl -w

# Small script to transform JASPAR database to JSON.
# See https://jaspar2022.genereg.net/api/v1/matrix/MA0564.1/ for a reference
# on the JSON format.

use strict;
use warnings;
use JSON;

local $/ = "\n>";  # Sets record split

my @pfms;

while (<>) {
    s/^>//;  # Remove heading ">" (only in first entry)
    s/>$//;  # Remove tailing ">" (all entries)
    next if /^\s*$/;  # Skip first (now empty) entry
    my @lines = split /\n/;
    my $header = shift @lines;
    my ($matrix_id, $name) = split /\s+/, $header, 2;
    my %pfm = (matrix_id => $matrix_id, name => $name, pfm => {});

    foreach my $line (@lines) {
        my ($base, $values) = $line =~ /^(\w)\s+\[\s*(.*?)\s*\]$/;
        if (!defined $base || !defined $values) {
            die "Ungültige Zeile: $line\n";
        }
        my @counts = split /\s+/, $values;
        $pfm{pfm}{$base} = \@counts;
    }

    push @pfms, \%pfm;
}

print to_json(\@pfms, { pretty => 1 });
