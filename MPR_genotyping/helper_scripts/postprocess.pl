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
        $blocks{$block_id}->{block_id} = $block_id;
        $blocks{$block_id}->{gtype} = $gtype;
        $blocks{$block_id}->{size} = $positions[$j-1]-$positions[$start_block]+1;
        $blocks{$block_id}->{num_obs} = $j-1-$start_block+1;
        $blocks{$block_id}->{start_idx} = $start_block;
        $blocks{$block_id}->{stop_idx} = $j-1;
        $blocks{$block_id}->{prev_block} = $prev_block_id;
        $gtype = $genotypes[$i]->[$j];
        if (defined $prev_block_id) {
            $blocks{$prev_block_id}->{next_block} = $block_id;
        }
        $prev_block_id = $block_id;
        $j--;
    }
    my $look_again = 1;
    while ($look_again) {
        #TODO: should I throw an error if no block initially exceeds the min size? 
        $look_again = 0;
        #my @size_sorted_blocks = sort {$b->{size} <=> $a->{size}} values %blocks;
        my @size_sorted_blocks = sort {$b->{num_obs} <=> $a->{num_obs}} values %blocks;
        for (my $j=0; $j < @size_sorted_blocks; $j++) {
            my $this_block = $size_sorted_blocks[$j];
            my $prev_block = $this_block->{prev_block};
            if (defined $prev_block && $blocks{$prev_block}->{size} < $min_block_size) {
                $prev_block = $blocks{$prev_block};
                my $new_gtype = $this_block->{gtype};
                my $prev_prev_block = $prev_block->{prev_block};
                if (defined $prev_prev_block && $blocks{$prev_prev_block}->{gtype} eq $this_block->{gtype}) {
                    $prev_prev_block = $blocks{$prev_prev_block};
                    my $start_idx = $prev_block->{start_idx};
                    my $stop_idx = $prev_block->{stop_idx};
                    #overwrite gtypes
                    @{$genotypes[$i]}[$start_idx..$stop_idx] = split //, $new_gtype x ($stop_idx-$start_idx+1);
                    $this_block->{prev_block} = $prev_prev_block->{prev_block};
                    $this_block->{start_idx} = $prev_prev_block->{start_idx};
                    $this_block->{size} = $positions[$this_block->{stop_idx}] - $positions[$this_block->{start_idx}] + 1;
                    $this_block->{num_obs} += $prev_prev_block->{num_obs} + $prev_block->{num_obs};
                    delete $blocks{$this_block->{block_id}};
                    $this_block->{block_id} = $positions[$this_block->{start_idx}]."-".$positions[$this_block->{stop_idx}];
                    if (defined $this_block->{prev_block}) {
                        $blocks{$this_block->{prev_block}}->{next_block} = $this_block->{block_id};
                    }
                    if (defined $this_block->{next_block}) {
                        $blocks{$this_block->{next_block}}->{prev_block} = $this_block->{block_id};
                    }
                    $blocks{$this_block->{block_id}} = $this_block;
                    delete $blocks{$prev_block->{block_id}};
                    delete $blocks{$prev_prev_block->{block_id}};
                    $look_again = 1;
                    last;
                }
            }
            my $next_block = $this_block->{next_block};
            if (defined $next_block && $blocks{$next_block}->{size} < $min_block_size) {
                $next_block = $blocks{$next_block};
                my $new_gtype = $this_block->{gtype};
                my $next_next_block = $next_block->{next_block};
                if (defined $next_next_block && $blocks{$next_next_block}->{gtype} eq $this_block->{gtype}) {
                    $next_next_block = $blocks{$next_next_block};
                    my $start_idx = $next_block->{start_idx};
                    my $stop_idx = $next_block->{stop_idx};
                    #overwrite gtypes
                    @{$genotypes[$i]}[$start_idx..$stop_idx] = split //, $new_gtype x ($stop_idx-$start_idx+1);
                    $this_block->{next_block} = $next_next_block->{next_block};
                    $this_block->{stop_idx} = $next_next_block->{stop_idx};
                    $this_block->{size} = $positions[$this_block->{stop_idx}] - $positions[$this_block->{start_idx}] + 1;
                    $this_block->{num_obs} += $next_next_block->{num_obs} + $next_block->{num_obs};
                    delete $blocks{$this_block->{block_id}};
                    $this_block->{block_id} = $positions[$this_block->{start_idx}]."-".$positions[$this_block->{stop_idx}];
                    if (defined $this_block->{next_block}) {
                        $blocks{$this_block->{next_block}}->{prev_block} = $this_block->{block_id};
                    }
                    if (defined $this_block->{prev_block}) {
                        $blocks{$this_block->{prev_block}}->{next_block} = $this_block->{block_id};
                    }
                    $blocks{$this_block->{block_id}} = $this_block;
                    delete $blocks{$next_block->{block_id}};
                    delete $blocks{$next_next_block->{block_id}};
                    $look_again = 1;
                    last;
                }
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

