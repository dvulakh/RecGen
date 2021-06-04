#!/usr/bin/perl

#####################################################################
# Reads an original and reconstructed poisson pedigree from STDIN
# and reads the REC-GEN data file given as an argument. Prints
# information about the sorts of errors made in recovering the
# parents of the extant population to STDOUT.
#####################################################################

use strict;
use warnings;
use v5.16;
use List::Util 'max';
use experimental 'smartmatch';

### PATH TO RECGEN ###
my $PATH_TO_RECGEN = @ARGV > 2 ? pop @ARGV : '.';

### PEDIGREE INFORMATION ###
# The couple ids associated with the original and reconstructed pedigrees
my %ID_TO_O; my %O_TO_ID;
my %ID_TO_R; my %R_TO_ID;
# Genome associated with each extant node
my %BLOCKS;
# Lists of actual and predicted siblings
my %SIB_O; my %SIB_R;
# Set of candidate pairs
my %SIB_PAIR;
# Set of hyperedge triples
my %SIB_TRIP;
# Set of missing edges
my %SIB_QUERY;

### TYPES OF MISTAKES ###
# Number of false negatives/false positives in pair detection
my $PAIR_NEG = 0; my $PAIR_POS = 0;
# Number of false negatives/false positives in triple detection
my $TRIP_NEG = 0; my $TRIP_POS = 0;
# Number of misplaced vertices (majority of reconstructed siblings are wrong)
my $CHANGELING = 0; my $ORPHAN = 0;
# Number of vertices where all children have same gene
my $LOST_GENES = 0;

### READ PEDIGREES ###
my $id_to = \%ID_TO_O;
my $to_id = \%O_TO_ID;
my $sib = \%SIB_O;
my $recon = 0;
while (<STDIN>) {
    # Switch to reading reconstructed pedigree on tilde
    if (/~/) {
        $recon = 1;
        $id_to = \%ID_TO_R;
        $to_id = \%R_TO_ID;
        $sib = \%SIB_R;
    }
    # Create new individual node entry
    elsif (/^i\s+-i\s+(\d+).*-g\s+\d+([\s\d]+)/) {
        $BLOCKS{$1} = [grep {$_} split /\s+/, $2] if $recon < 1;
    }
    # Process a couple
    elsif (/^c\s+-i\s+(\d+).*-m\s+2\s+(\d+)\s+(\d+).*-c\s+\d+([\s\d]+)/) {
        # Rebind captures
        my $cid = $1, my $m1 = $2, my $m2 = $3, my $ch_raw = $4;
        chomp $ch_raw;
        my @ch = grep {$_} split /\s+/, $ch_raw;
        # If members are the same, this is the couple of an extant node
        if ($m1 == $m2) {
            $$id_to{$m1} = $cid;
            $$to_id{$cid} = $m1;
        }
        # Otherwise, this is some ancestor -- add children as siblings
        else {
            for my $v (@ch) {
                $$sib{$v} = [ @ch ];
            }
        }
    }
}

### READ DATA FILE ###
my $finished_ext = 0;
while (<>) {
    # Only process first generation of rec-gen data file
    if ($finished_ext < 2) {
        # Process a candidate pair
        if (/Found candidate pair \((\d+),\s*(\d+)\)/) {
            my $u = $R_TO_ID{$1}, my $v = $R_TO_ID{$2};
            if ($finished_ext != 0) {
                $finished_ext = 2;
                next;
            }
            $SIB_PAIR{canon($u, $v)} = undef;
            $PAIR_POS++ if !($u ~~ $SIB_O{$v});
        }
        # Process a hyperedge triple
        elsif (/Inserting hypergraph edge \((\d+),\s*(\d+),\s*(\d+)\)/) {
            my $u = $R_TO_ID{$1}, my $v = $R_TO_ID{$2}, my $w = $R_TO_ID{$3};
            $SIB_TRIP{canon($u, $v, $w)} = undef;
            $TRIP_POS++ if !($u ~~ $SIB_O{$v} && $w ~~ $SIB_O{$v});
        }
        # Process gene summary
        elsif (/found genes \d+ and \d+ \(frequency:\s*(\d+)\s+(\d+)\)/) {
            $LOST_GENES++ if $2 == 0;
            $finished_ext = 1;
        }
    }
    # Process lines from tree-diff data file
    if (/Missing edge.*to \((\d+)o/) {
        $SIB_QUERY{$O_TO_ID{$1}} = undef if exists $O_TO_ID{$1};
    }
}

### COUNT FALSE NEGATIVES ###
for my $i (keys %ID_TO_O) {
    for my $j (@{$SIB_O{$i}}) {
        if ($i < $j) {
            $PAIR_NEG++ if !exists $SIB_PAIR{canon($i, $j)};
            for my $k (@{$SIB_O{$i}}) {
                $TRIP_NEG++ if $j < $k && !exists $SIB_TRIP{canon($i, $j, $k)};
            }
        }
    }
}

### COUNT NUMBER OF MISPLACED VERTICES ###
for my $i (keys %ID_TO_O) {
    if (!exists $SIB_R{$i}) {
        $ORPHAN++;
    } else {
        $CHANGELING++ if (grep {$_ ~~ $SIB_O{$i}} @{$SIB_R{$i}}) < @{$SIB_R{$i}} / 2 ||
            (grep {$_ ~~ $SIB_R{$i}} @{$SIB_O{$i}}) < @{$SIB_O{$i}} / 2;
    }
}

### OUTPUT SUMMARY ###
print "Pairs false negatives:   $PAIR_NEG\n";
print "Pairs false positives:   $PAIR_POS\n";
print "Triples false negatives: $TRIP_NEG\n";
print "Triples false positives: $TRIP_POS\n";
print "Changeling vertex:       $CHANGELING\n";
print "Orphaned vertex:         $ORPHAN\n";
print "Lost genes:              $LOST_GENES\n";

### PROCESS QUERIES ###
my $qout = '';
for my $q (keys %SIB_QUERY) {
    my $nsib = $#{$SIB_O{$q}};
    my $nrec = max((grep {$_ ~~ $SIB_O{$q}} @{$SIB_R{$q}}) - 1, 0);
    my $nmis = grep {!($_ ~~ $SIB_R{$q})} @{$SIB_O{$q}};
    my $nwsb = grep {!($_ ~~ $SIB_O{$q})} @{$SIB_R{$q}};
    $qout .= "Vertex $q: siblings $nsib reconstructed $nrec missing $nmis incorrect $nwsb\n";
}
print `printf '$qout' | $PATH_TO_RECGEN/outfmt 3`;

### CANONIZATION HELPER ###
sub canon {
    return join ',', sort @_;
}
