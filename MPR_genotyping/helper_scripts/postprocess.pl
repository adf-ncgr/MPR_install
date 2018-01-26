#!/usr/bin/env perl
my @headers;
my $count=0;
my @genotypes;
#FIXME: parameterize
my $min_block_size = 50000;
my $chr;
while (<>) {
    print && next if /^#/;
    chomp;
    if (/^CHROM\tPOS/) {
        (undef, undef, @headers) = split /\t/;
        print $_, "\n";
    }
    else {
        my @data;
        my $pos;
        ($chr, $pos, @data) = split /\t/;
        $positions[$count] = $pos;
        for (my $i=0; $i < @data; $i++) {
            $genotypes[$i]->[$count] = $data[$i];
        }
        $count++;
    }
}

for (my $i=0; $i < @genotypes; $i++) {
    my %blocks;
    my $gtype = $genotypes[$i]->[0];
    my $prev_block_id;
    for (my $j=0; $j < @{$genotypes[$i]}; $j++) {
        my $start_block = $j;
        while ($genotypes[$i]->[$j] == $gtype) {
            $j++;
        }
        my $block_id = $positions[$start_block]."-".$positions[$j-1];
        $blocks{$block_id}->{gtype} = $gtype;
        $blocks{$block_id}->{size} = $positions[$j-1]-$positions[$start_block]+1;
        $blocks{$block_id}->{start_idx} = $start_block;
        $blocks{$block_id}->{stop_idx} = $j-1;
        $blocks{$block_id}->{prev_block} = $prev_block_id;
        $gtype = $genotypes[$i]->[$j];
        if (defined $prev_block_id) {
            $blocks{$prev_block_id}->{next_block} = $block_id;
        }
        $prev_block_id = $block_id;
    }
    my @size_sorted_blocks = sort {$b->{size} <=> $a->{size}} values %blocks;
    my $look_again = 1;
    while ($look_again) {
        #TODO: should I throw an error if no block initially exceeds the min size? 
        $look_again = 0;
        for (my $j=0; $j < @size_sorted_blocks; $j++) {
            my $prev_block = $size_sorted_blocks[$j]->{prev_block};
            if (defined $prev_block && $blocks{$prev_block}->{size} < $min_block_size) {
                my $new_gtype = $size_sorted_blocks[$j]->{gtype};
                my $start_idx = $blocks{$prev_block}->{start_idx};
                my $stop_idx = $blocks{$prev_block}->{stop_idx};
                #overwrite gtypes
                @{$genotypes[$i]}[$start_idx..$stop_idx] = split //, $new_gtype x ($stop_idx-$start_idx+1);
                my $new_prev_block = $blocks{$prev_block}->{prev_block};
                #extend to neighboring block if overwritten block was between two of same gtype
                if ($blocks{$new_prev_block}->{gtype} eq $new_gtype) {
                    delete $blocks{$prev_block};
                    $prev_block = $new_prev_block;
                }
                $size_sorted_blocks[$j]->{prev_block} = $new_prev_block;
                $size_sorted_blocks[$j]->{start_idx} = $blocks{$prev_block}->{start_idx};
                $size_sorted_blocks[$j]->{size} = $positions[$size_sorted_blocks[$j]->{stop_idx}] - $positions[$size_sorted_blocks[$j]->{start_idx}] + 1;
                delete $blocks{$prev_block};
                $look_again = 1;
            }
            my $next_block = $size_sorted_blocks[$j]->{next_block};
            if (defined $next_block && $blocks{$next_block}->{size} < $min_block_size) {
                my $new_gtype = $size_sorted_blocks[$j]->{gtype};
                my $start_idx = $blocks{$next_block}->{start_idx};
                my $stop_idx = $blocks{$next_block}->{stop_idx};
                #overwrite gtypes
                @{$genotypes[$i]}[$start_idx..$stop_idx] = split //, $new_gtype x ($stop_idx-$start_idx+1);
                my $new_next_block = $blocks{$next_block}->{next_block};
                #extend to neighboring block if overwritten block was between two of same gtype
                if ($blocks{$new_next_block}->{gtype} eq $new_gtype) {
                    delete $blocks{$next_block};
                    $next_block = $new_next_block;
                }
                $size_sorted_blocks[$j]->{next_block} = $new_next_block;
                $size_sorted_blocks[$j]->{stop_idx} = $blocks{$next_block}->{stop_idx};
                $size_sorted_blocks[$j]->{size} = $positions[$size_sorted_blocks[$j]->{stop_idx}] - $positions[$size_sorted_blocks[$j]->{start_idx}] + 1;
                delete $blocks{$next_block};
                $look_again = 1;
            }
            if ($look_again) {
                @size_sorted_blocks = sort {$b->{size} <=> $a->{size}} values %blocks;
                last;
            }
        }
    }
}

for (my $i = 0; $i < @positions; $i++) {
    print $chr,"\t",$positions[$i];
    for (my $j = 0; $j < @genotypes; $j++) {
           print "\t", $genotypes[$j]->[$i];
    }
    print "\n";
}

