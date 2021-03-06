#!/usr/bin/env perl
use strict;

my $correction_output = shift;
open(C, $correction_output) || die $!;
#for now, we can assume the files can be read in parallel (ie rows and columns are consistent between the two, modulo some missing headers and non-essential columns)
#skip head
<C>;

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
                next;
            }
            else {
                chomp;
                my @data = split /\t/;
                my ($chr, $pos) = split /_/, $c_gt[0];
                next unless $chr eq $data[0] && $pos == $data[1];
                $matched_it = 1;
                for (my $i=1; $i < @c_gt; $i++) {
                    my $gt;
                    if ($c_gt[$i] == 0) {
                        $gt = "./.";
                    }
                    #FIXME: this logic only works when the reference is always identical to parent1
                    elsif ($c_gt[$i] == 1) {
                        $gt = "0/0";
                    }
                    elsif ($c_gt[$i] == 2) {
                        $gt = "1/1";
                    }
                    elsif ($c_gt[$i] == 3) {
                        $gt = "1/0";
                    }
                    $data[$i-2+GT_START_IDX] =~ s/^[^:]+/$gt/;
                }
                print join("\t",@data),"\n";
            }
        }
}


