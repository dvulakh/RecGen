#!/usr/bin/perl

#####################################################################
# Reads a pedigree from STDIN and outputs a pedigree with identical
# structure of parent-child relationships but non-extant genomes set
# to all-zeroes
#####################################################################

use strict;
use warnings;
use v5.16;

my %ext;
my @log = <STDIN>;
for (@log) {
	chomp;
	$ext{$1} = undef if (/^c.*-m\s+2\s+(\d+)\s+(\d+)/ && $1 eq $2);
}
for (@log) {
	chomp;
	if (/^i\s+-i\s+(\d+)(.*)-g\s+(\d+)\s+(.*)$/) {
		print "i -i $1$2 -g $3 ".
			(exists $ext{$1} ? $4 : join ' ', (0) x $3) . "\n";
	} else { print "$_\n"; }
}
