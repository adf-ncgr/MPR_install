#!/usr/bin/env perl

use strict;
use Getopt::Long;
use constant GT_START_IDX => 9;
my @parents;
my $assume_ref_diffs_are_parent2=0;
my $rescue_sites_with_single_parent_coverage=0;
my $coordinate_map;
my $min_ns;
my $min_dp;
GetOptions(
    "parent=s" => \@parents,
    "rescue_sites_with_single_parent_coverage!" => \$rescue_sites_with_single_parent_coverage,
    "coordinate_map=s" => \$coordinate_map,
    "min_ns=i" => \$min_ns,
    "min_dp=i" => \$min_dp,
);
if (@parents == 0) {
    #this implies ref is parent1
    $assume_ref_diffs_are_parent2=1;
}
elsif (@parents > 2) {
    die "must supply <= 2 --parent\n";
}
my %coordinate_map;
if (defined $coordinate_map) {
    open(CM, $coordinate_map) || die $!;
    %coordinate_map = map {chomp; /([^\t]+)\t(.*)/;} <CM>;
    close CM;
}
my $line;
while ($line = <>) {
    last if $line =~ /^#CHROM/;
}
$line =~ s/^#//;
chomp $line;
my @headers = split /\t/, $line;
print join("\t", @headers[0,1,GT_START_IDX..$#headers]), "\n";
my $parent1_idx = -1;
my $parent2_idx = -1;
if (@parents > 0) {
    for (my $i=GT_START_IDX; $i < @headers; $i++) {
        if ($headers[$i] eq $parents[0]) {
            $parent1_idx = $i;
        }
        elsif ($headers[$i] eq $parents[1]) {
            $parent2_idx = $i;
        }
    }
}
while (<>) {
    next if /^#/;
    chomp;
    my @data = split /\t/;
    if (defined $min_ns) {
        my ($ns) = ($data[GT_START_IDX-2] =~ /NS=(\d+)/);
        if ($ns < $min_ns) {
            next;
        }
    }
    if (defined $min_dp) {
        my ($dp) = ($data[GT_START_IDX-2] =~ /DP=(\d+)/);
        if ($dp < $min_dp) {
            next;
        }
    }
    my $genotype_parent1;
    my $genotype_parent2;
    if ($assume_ref_diffs_are_parent2) {
        $genotype_parent1 = "0/0";
        $genotype_parent2 = "1/1";
    }
    else {
        #parent1 has to be assigned (by def)
        ($genotype_parent1) = split /:/, $data[$parent1_idx];
        $genotype_parent1 =~ s/\|/\//;
        if ($parent2_idx == -1) {
            #$genotype_parent2 = "0/0";
            $genotype_parent2 = "./.";
        }
        else {
            ($genotype_parent2) = split /:/, $data[$parent2_idx];
            $genotype_parent2 =~ s/\|/\//;
        }
    }
    if ($rescue_sites_with_single_parent_coverage) {
        if ($genotype_parent1 eq "./.") {
            if ($genotype_parent2 eq "1/1") {
                $genotype_parent1 = "0/0";
            }
            elsif ($genotype_parent2 eq "0/0") {
                $genotype_parent1 = "1/1";
            }
        }
        elsif ($genotype_parent2 eq "./.") {
            if ($genotype_parent1 eq "1/1") {
                $genotype_parent2 = "0/0";
            }
            elsif ($genotype_parent1 eq "0/0") {
                $genotype_parent2 = "1/1";
            }
        }
    }
    #for now, at least, let's limit to cases in which parents are homozygous and contrasted
    if (!(join("/", sort($genotype_parent1, $genotype_parent2)) eq "0/0/1/1")) {
        #not printing will cause problems with vcf correction, but correcting would be problematic too, and these seem
        #to be making the algorithm for gt correction behave badly in any case
    #print $data[0],"\t",$data[1];
        #print "\t0"x(scalar(@data)-9);
    }
    else {
        if (defined $coordinate_map) {
            my $mapped_coord = $coordinate_map{$data[0].":".$data[1]};
            #comment out this check due to the chimera trimming, but possibly worth revisiting
            if (! defined $mapped_coord) {
                #die "could not find coordinate mapping for $data[0]:$data[1]\n" unless defined $mapped_coord;
                warn "could not find coordinate mapping for $data[0]:$data[1]\n" unless defined $mapped_coord;
                next;
            }
            else {
                print $mapped_coord;
            }
        }
        else {
            print $data[0],"\t",$data[1];
        }
        #there may be cases in which both parents differ from the ref
        #TODO: will ordering sometimes be an issue?
        my $het = join("/", sort {$a <=> $b} (substr($genotype_parent1,0,1), substr($genotype_parent2,0,1)));
        for (my $i = 9; $i < @data; $i++) {
            my ($gt) = split /:/, $data[$i];
            $gt =~ s/\|/\//;
            if ($gt eq $genotype_parent1) {
                print "\t1";
            }
            elsif ($gt eq $genotype_parent2) {
                print "\t2";
            }
            elsif ($gt eq $het) {
                print "\t3";
            }
            else { 
                #not sure about this, but the hmm function errors when NAs are supplied
                print "\t0";
            }
        }
        print "\n";
    }
}
