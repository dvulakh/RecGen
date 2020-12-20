#!/usr/bin/perl

#####################################################################
# Reads an original and reconstructed poisson pedigree from STDIN
# and reads the treediff data file given as an argument. Prints a
# pedigree with topology identical to the original pedigree, in
# which the genome of each couple is the genome of its image in the
# reconstructed pedigree (or all 0s if the couple was not recovered)
#####################################################################

use strict;
use warnings;
use v5.16;

### PATH TO RECGEN ###
my $PATH_TO_RECGEN = @ARGV > 1 ? pop @ARGV : '.';

### PEDIGREE INFORMATION ###
# Bijection
my %O_TO_R;
# Lists of recovered genes for each couple (pair)
my %R_TO_GEN;
# Original pedigree string
my $PED;

### READ PEDIGREES ###
my $SW = 0;
for (<STDIN>) {
	chomp;
	# Switch to reading reconstructed pedigree on tilde
	$SW = 1 if /~/;
	# Store original pedigree
	if ($SW == 0) {
		$PED .= "$_\n";
	}
	# Read genome for recovered pedigree
	elsif (/^i.*-c\s+(\d+).*-g\s+(\d+)\s+([\d\s]+)/) {
		push @{$R_TO_GEN{$1}}, [ split /\s+/, $3 ];
	}
}

### READ BIJECTION ###
for (<>) { $O_TO_R{$1} = $2 if /bijection\s+\((\d+)o[^\d]*(\d+)r\)/; }

### SPLICE GENOMES ONTO ORIGINAL ###
for (split /\n/, $PED) {
	if (/^(i.*-c\s+(\d+).*-g\s+(\d+)\s+)([\d\s]+)$/) {
		my $l = $1; my $c = $2; my $ng = $3; my @g;
		if (exists $O_TO_R{$c} && exists $R_TO_GEN{$O_TO_R{$c}} && @{$R_TO_GEN{$O_TO_R{$c}}} > 0) {
			@g = @{pop @{$R_TO_GEN{$O_TO_R{$c}}}};
		} else {
			@g = (0) x $ng;
		}
		print $l.(join ' ', @g)."\n";
	} else { print "$_\n"; }
}
