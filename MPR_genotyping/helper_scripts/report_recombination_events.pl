#!/usr/bin/env perl
use strict;
my $header;
while ($header = <>) {
    next if $header =~ /^#/;
    last;
}
chomp $header;
my (undef, @headers) = split /\t/, $header;

my ($last_genotypes, $last_sequence, %recombinations, %line_recombinations, @last_genotypes, $last_recombination_pos, %interrecomb_dists, %genotypes2contigs);

while (<>) {
    chomp;
    my ($sequence, $pos, $genotypes) = /^([^\t]+)_(\d+)\t(.*)/;
    if (!($sequence eq $last_sequence)) {
        if (defined $last_sequence) {
            $genotypes2contigs{&convert_to_binary($last_genotypes)}->{$last_sequence.":".$last_recombination_pos."-END"} = 1;
        }
        $recombinations{$sequence} = 0;
        $interrecomb_dists{$sequence} = [];
        $last_recombination_pos = 1;
        for (my $i=0; $i < @headers; $i++) {
            $line_recombinations{$sequence}->{$headers[$i]} = 0;
        }
        @last_genotypes = split /\t/, $genotypes;
        $last_genotypes = $genotypes;
    }
    elsif (!($genotypes eq $last_genotypes)) {
        if (!defined $last_recombination_pos) {
            $genotypes2contigs{&convert_to_binary($last_genotypes)}->{$sequence.":BEGIN-".$pos} = 1;
        }
        else {
            $genotypes2contigs{&convert_to_binary($last_genotypes)}->{$sequence.":".$last_recombination_pos."-".$pos} = 1;
        }
        #push @{$interrecomb_dists{$sequence}}, ($pos-$last_recombination_pos);
        #$last_recombination_pos = $pos;
        $recombinations{$sequence}++;
        my @genotypes = split /\t/, $genotypes;
        for (my $i=0; $i < @headers; $i++) {
            if ($genotypes[$i] cmp $last_genotypes[$i]) {
                push @{$interrecomb_dists{$sequence}}, ($pos-$last_recombination_pos);
                $last_recombination_pos = $pos;
                $line_recombinations{$sequence}->{$headers[$i]}++;
            }
        }
        @last_genotypes = split /\t/, $genotypes;
    }
    $last_sequence = $sequence;
    $last_genotypes = $genotypes;
}

print "#per-sequence recombination counts\n";
print join("\n", map {$_."\t".$recombinations{$_};} sort {$recombinations{$b} <=> $recombinations{$a}} keys %recombinations);
print "\n";
print "#per-sequence inter-recombination distances\n";
print join("\n", map {$_."\t".join(" ", @{$interrecomb_dists{$_}});} sort {$recombinations{$b} <=> $recombinations{$a}} keys %recombinations);
print "\n";
print "#per-line recombination counts\n";
print $header, "\n";
foreach my $sequence (sort {$recombinations{$b} <=> $recombinations{$a}} keys %recombinations) {
    print join("\t", $sequence, map {$line_recombinations{$sequence}->{$_};} @headers), "\n";
}

print "#suggested linkages\n";
foreach my $g (keys %genotypes2contigs) {
    my @contigs = keys %{$genotypes2contigs{$g}};
    if (scalar(@contigs) > 1) {
        print join("\t", @contigs), "\n";
    }
}

my @genotypes = keys %genotypes2contigs;
#(@genotypes) = &convert_to_binary(@genotypes);
my $threshold_distance = 2;
for (my $i = 0; $i < @genotypes; $i++) {
    my $min_dist = ~0;
    my @contigs = ();
    for (my $j = $i+1; $j < @genotypes; $j++) {
        my $distance = &genotype_distance($genotypes[$i], $genotypes[$j]);
        if ($distance < $min_dist) {
            $min_dist = $distance;
            @contigs = ();
        }
        if ($distance == $min_dist) {
            push @contigs, keys %{$genotypes2contigs{$genotypes[$j]}};
        }
        if ($distance <= $threshold_distance) {
            print "at distance $distance: " . join(",", keys %{$genotypes2contigs{$genotypes[$i]}}) . " with " . join(",", keys %{$genotypes2contigs{$genotypes[$j]}}) . "\n";
        }
    }
    print "min_dist=$min_dist for ".join(",", keys %{$genotypes2contigs{$genotypes[$i]}})
        .  " " . ($min_dist > $threshold_distance ? join(",", @contigs) : "")
        ."\n";
}

sub convert_to_binary() {
    #return map {
        my ($genotypes) = @_;
        my @genotypes = split //, $genotypes;
        my $retval = '';
        for (my $i = 0; $i < @genotypes; $i++) {
            if ($genotypes[$i] eq "A") {
                vec($retval,$i,2) = 0b00;
            }
            elsif ($genotypes[$i] eq "B") {
                vec($retval,$i,2) = 0b11;
            }
            else {
                vec($retval,$i,2) = 0b01;
            }
        }
        return $retval;
    #} @_;
}

sub genotype_distance() {
    my ($g1, $g2) = @_;
    #my @bits = split(//, unpack("b*", ($g1^$g2)));
    #my $distance = 0;
    #map {$distance += $_;} @bits;
    my $bits = unpack("b*", ($g1^$g2));
    my $distance =()= $bits =~ /1/g;
    return $distance;
}
