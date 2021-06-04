#!/usr/bin/perl

#####################################################################
# Reads a pedigree from STDIN and outputs a pedigree with identical
# structure of parent-child relationships but genomes shortened to
# the command-line argument
#####################################################################

use strict;
use warnings;
use v5.16;

my $LEN = shift;
for (<STDIN>) {
    chomp;
    if (/^(.*-i.*)-g\s+(\d+)\s+([^\-]*)(.*)$/) {
        my @gen = grep {$_} split /\s+/, $3;
        print "$1-g $LEN @gen[0..$LEN-1]$4\n";
    } elsif (/^(.*)-B\s+(\d+)(.*)$/) { print "$1-B $LEN$3\n"; }
    else { print "$_\n"; }
}
