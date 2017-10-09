#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $uncorrected=0;
#when comparing before/after correction, it is problematic for the simple 
#"join" if they are imbalanced. This default seeks to address that
my $skip_uncorrected_if_corrected_empty=1;
GetOptions(
    "uncorrected" => \$uncorrected,
    "skip_uncorrected_if_corrected_empty!" => \$skip_uncorrected_if_corrected_empty
    );
my $suffix = ($uncorrected ? ".gt_table" : ".mpr_corrected");
my $wrote_header=0;
my $comment_headers="";
while (<>) {
    chomp;
    if ($uncorrected && $skip_uncorrected_if_corrected_empty) {
        if (-z $_.".mpr_corrected") {
            next;
        }
    }
    open(F,$_.$suffix) || die $!;
    while (<F>) {
        if (/^#/) {
           next if $wrote_header;
           $comment_headers .= $_; 
        }
        elsif (/^CHROM/) {
            next if $wrote_header;
            print $comment_headers;
            my @headers = split /\t/;
            print "CHROM_POS\t", join("\t", @headers[2..$#headers]);
            $wrote_header=1;
        }
        else {
            chomp;
            my ($chrom, $pos, @data) = split /\t/;
            my $data = join "", @data;
            $data =~ tr/123/BAH/;
            print join("\t", "${chrom}_${pos}", split(//, $data)), "\n";
        }
    }
    close F;
}
