#!/usr/bin/env perl
use strict;
my $header;
while ($header = <>) {
    next if $header =~ /^#/;
    last;
}
chomp $header;
my (undef, @headers) = split /\t/, $header;

my ($last_genotypes, $last_sequence, $last_pos, %recombinations, %line_recombinations, @last_genotypes, $last_recombination_pos, %interrecomb_dists, %genotypes2blocks, %blocks2markers, @block_marker_positions);

while (<>) {
    chomp;
    my ($sequence, $pos, $genotypes) = /^([^\t]+)_(\d+)\t(.*)/;
    if (!($sequence eq $last_sequence)) {
        if (defined $last_sequence) {
            my $block = $last_sequence.":".$last_recombination_pos."-END";
            $genotypes2blocks{&convert_to_binary($last_genotypes)}->{$block} = 1;
            $blocks2markers{$block} = [$block_marker_positions[0], $block_marker_positions[$#block_marker_positions/2], $block_marker_positions[$#block_marker_positions]];
            @block_marker_positions = ();
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
            my $block = $sequence.":BEGIN-".$last_pos;
            $genotypes2blocks{&convert_to_binary($last_genotypes)}->{$block} = 1;
            $blocks2markers{$block} = [$block_marker_positions[0], $block_marker_positions[$#block_marker_positions/2], $block_marker_positions[$#block_marker_positions]];
            @block_marker_positions = ();
        }
        else {
            my $block = $sequence.":".$last_recombination_pos."-".$last_pos;
            $genotypes2blocks{&convert_to_binary($last_genotypes)}->{$block} = 1;
            $blocks2markers{$block} = [$block_marker_positions[0], $block_marker_positions[$#block_marker_positions/2], $block_marker_positions[$#block_marker_positions]];
            @block_marker_positions = ();
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
    $last_pos = $pos;
    push @block_marker_positions, $pos;
}
my $block = $last_sequence.":".$last_recombination_pos."-END";
$genotypes2blocks{&convert_to_binary($last_genotypes)}->{$block} = 1;
$blocks2markers{$block} = [$block_marker_positions[0], $block_marker_positions[$#block_marker_positions/2], $block_marker_positions[$#block_marker_positions]];

print "#VERSION postprocessing_recombination_initial-20-g99ff0e4-26;\n";
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

my %good_blocks;
print "#suggested linkages\n";
foreach my $g (keys %genotypes2blocks) {
    my @blocks = keys %{$genotypes2blocks{$g}};
    if (scalar(@blocks) > 1) {
        print join("\t", @blocks), "\n";
        #FIXME: maybe too naive- should look for presence of inter-contig blocks?
        map {$good_blocks{$_} = 1;} @blocks;
    }
}

my @genotypes = keys %genotypes2blocks;
#(@genotypes) = &convert_to_binary(@genotypes);
my $threshold_distance = 2;
my %min_dists;
for (my $i = 0; $i < @genotypes; $i++) {
    #may have previously been compared in earlier rounds (when it was a j to some other i)
    my @blocks;
    my $min_dist = $min_dists{$genotypes[$i]}->{dist};
    if (! defined $min_dist) {
        $min_dist = ~0;
        @blocks = ();
    }
    else {
        @blocks = @{$min_dists{$genotypes[$i]}->{blocks}};
    }
    for (my $j = $i+1; $j < @genotypes; $j++) {
        my $distance = &genotype_distance($genotypes[$i], $genotypes[$j]);
        if ($distance < $min_dist) {
            $min_dist = $distance;
            @blocks = ();
            $min_dists{$genotypes[$i]}->{dist} = $min_dist;
            $min_dists{$genotypes[$i]}->{blocks} = \@blocks;
        }
        #this will kick in if the above conditional reset min_dist or if we happen to equal a previously set min_dist
        if ($distance == $min_dist) {
            push @blocks, keys %{$genotypes2blocks{$genotypes[$j]}};
        }
        if ($distance <= $threshold_distance) {
            print "at distance $distance: " . join(",", keys %{$genotypes2blocks{$genotypes[$i]}}) . " with " . join(",", keys %{$genotypes2blocks{$genotypes[$j]}}) . "\n";
            map {$good_blocks{$_} = 1;} keys %{$genotypes2blocks{$genotypes[$i]}};
            map {$good_blocks{$_} = 1;} keys %{$genotypes2blocks{$genotypes[$j]}};
        }
    }
    print "min_dist=$min_dist for ".join(",", keys %{$genotypes2blocks{$genotypes[$i]}})
        .  " " . ($min_dist > $threshold_distance ? join(",", @blocks) : "")
        ."\n";
}

print "#good markers\n";
print $header,"\n";
foreach my $genotype (keys %genotypes2blocks) {
    my $genotype_string;
    foreach $block (keys %{$genotypes2blocks{$genotype}}) {
        next unless $good_blocks{$block};
        if (! defined $genotype_string) {
            $genotype_string = substr(&convert_from_binary($genotype),0,scalar(@headers));
            $genotype_string = join("\t", split(//, $genotype_string));
        }
        my ($contig) = ($block =~ /^([^:]+)/);
        my $genotypes = $block;
        ##I think this was based on a request from Joe Curley, but I don't think her cares anymore so changing it
        #for (my $i = 0; $i < 3; $i++) {
            #print "${contig}_$blocks2markers{$block}->[$i]\t$genotype_string\n";
        #}
        print "${contig}_${block}\t$genotype_string\n";
    }
}

sub convert_to_binary() {
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
}

sub convert_from_binary() {
    my ($genotypes) = @_;
    my @genotype_codes = ("A", "H", "H", "B");
    my @bits = split(//, unpack("b*", $genotypes));
    my $retval = "";
    for (my $i = 0; $i < @bits; $i += 4) {
        $retval .= $genotype_codes[oct("0b".join("", $bits[$i], $bits[$i+1]))];
    }
    return $retval;
}

sub genotype_distance() {
    my ($g1, $g2) = @_;
    #my $bits = unpack("b*", ($g1^$g2));
    #my $distance =()= $bits =~ /1/g;
    my $distance = unpack("%32b*", ($g1^$g2));
    return $distance;
}
