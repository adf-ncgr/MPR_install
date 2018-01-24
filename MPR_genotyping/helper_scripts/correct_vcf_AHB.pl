#!/usr/bin/env perl
use strict;
use Getopt::Long;
my $parentA;
my $parentB;
GetOptions(
    "parentA=s" => \$parentA,
    "parentB=s" => \$parentB,
);
my $correction_output = shift;
open(C, $correction_output) || die $!;
#for now, we can assume the files can be read in parallel (ie rows and columns are consistent between the two, modulo some missing headers and non-essential columns)
#skip head
while (<C>) {
    last if /^CHROM_POS/;
}

my $parentA_idx;
my $parentB_idx;
use constant GT_START_IDX => 9;
while (my $c_gt = <C>) {
        chomp $c_gt;
        my @c_gt = split /\t/, $c_gt;
        my $matched_it = 0;
        while (! $matched_it) {
            $_ = <>;
            die "could not match $c_gt something must be wrong!\n" unless $_;
            if (/^#/) {
                print;
                if (/^#CHROM/) {
                    chomp;
                    my @headers = split /\t/;
                    for (my $i = GT_START_IDX; $i < @headers; $i++) {
                        if ($headers[$i] eq $parentA) {
                            $parentA_idx = $i;
                        }
                        elsif ($headers[$i] eq $parentB) {
                            $parentB_idx = $i;
                        }
                    }
                    if (defined $parentA && !defined $parentA_idx) {
                        die "could not find $parentA in VCF headers\n";
                    }
                    if (defined $parentB && !defined $parentB_idx) {
                        die "could not find $parentB in VCF headers\n";
                    }
                }
                next;
            }
            else {
                chomp;
                my @data = split /\t/;
                my ($chr, $pos) = ($c_gt[0] =~ /(.*)_(\d+)/);
                next unless $chr eq $data[0] && $pos == $data[1];
                $matched_it = 1;
                my $parentA_gt;
                if (defined $parentA_idx) {
                    ($parentA_gt) = ($data[$parentA_idx] =~ /^([^:]+)/);
                    if ($parentA_gt eq ".") {
                        $parentA_gt = undef;
                    }
                }
                my $parentB_gt;
                if (defined $parentB_idx) {
                    ($parentB_gt) = ($data[$parentB_idx] =~ /^([^:]+)/);
                    if ($parentB_gt eq ".") {
                        $parentB_gt = undef;
                    }
                }
                #handle cases where one is known and other is not by assuming contrast
                if (defined $parentA_gt && !defined $parentB_gt) {
                    if ($parentA_gt eq "0/0") {
                        $parentB_gt = "1/1";
                    }
                    elsif ($parentA_gt eq "1/1") {
                        $parentB_gt = "0/0";
                    }
                }
                elsif (defined $parentB_gt && !defined $parentA_gt) {
                    if ($parentB_gt eq "0/0") {
                        $parentA_gt = "1/1";
                    }
                    elsif ($parentB_gt eq "1/1") {
                        $parentA_gt = "0/0";
                    }
                }
                #in this case, we assume that the reference was one parent; for historical
                #reasons, ref => B (because ALT => A)
                if (!defined $parentA_idx && !defined $parentB_idx) {
                    $parentA_gt = "1/1";
                    $parentB_gt = "0/0";
                }
                for (my $i=1; $i < @c_gt; $i++) {
                    my $gt;
                    if ($c_gt[$i] eq "B") {
                        $gt = $parentB_gt;
                    }
                    elsif ($c_gt[$i] eq "A") {
                        $gt = $parentA_gt;
                    }
                    elsif ($c_gt[$i] eq "H") {
                        $gt = "1/0";
                    }
                    $data[$i-1+GT_START_IDX] =~ s/^[^:]+/$gt/;
                }
                print join("\t",@data),"\n";
            }
        }
}


