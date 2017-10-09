#!/usr/bin/env perl
use strict;
my $total;
my %class_symbols = (0 => "NA", 1 => "HOM_P1", 2 => "HOM_P2", 3 => "HET");
my %class_counts;
<>;
while (<>) {
    next if /^#/;
    chomp;
    my @data = split /\t/;
    for (my $i = 2; $i < @data; $i++) {
        $total++;
        $class_counts{$data[$i]}++;
    }
}

foreach my $class (sort keys %class_counts) {
    print $class_symbols{$class} , "\t", sprintf("%0.8f", $class_counts{$class}/$total), "\n";
}
