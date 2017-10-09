#!/usr/bin/env perl
use strict;
my $header;
while ($header = <>) {
    next if $header =~ /^#/;
    last;
}
chomp $header;
my (undef, @headers) = split /\t/, $header;

my ($last_genotypes, $last_sequence, %recombinations, %line_recombinations, @last_genotypes);

while (<>) {
    chomp;
    my ($sequence, $genotypes) = /([^\t]*)\t(.*)/;
    $sequence =~ s/_\d+$//;
    if (!($sequence eq $last_sequence)) {
        $recombinations{$sequence} = 0;
        for (my $i=0; $i < @headers; $i++) {
            $line_recombinations{$sequence}->{$headers[$i]} = 0;
        }
        @last_genotypes = split /\t/, $genotypes;
    }
    elsif (!($genotypes eq $last_genotypes)) {
        $recombinations{$sequence}++;
        my @genotypes = split /\t/, $genotypes;
        for (my $i=0; $i < @headers; $i++) {
            if ($genotypes[$i] cmp $last_genotypes[$i]) {
                $line_recombinations{$sequence}->{$headers[$i]}++;
            }
        }
        @last_genotypes = split /\t/, $genotypes;
    }
    $last_sequence = $sequence;
    $last_genotypes = $genotypes;
}

print join("\n", map {$_."\t".$recombinations{$_};} sort {$recombinations{$b} <=> $recombinations{$a}} keys %recombinations);
print "\n";
print $header, "\n";
foreach my $sequence (sort {$recombinations{$b} <=> $recombinations{$a}} keys %recombinations) {
    print join("\t", $sequence, map {$line_recombinations{$sequence}->{$_};} @headers), "\n";
}
