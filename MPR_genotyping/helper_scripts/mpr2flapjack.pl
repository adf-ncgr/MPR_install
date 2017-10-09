#!/usr/bin/env perl
use strict;
use Getopt::Long;
use constant ALT_ALLELE => "A";
use constant REF_ALLELE => "B";
use constant HET => "H";
my $phenotype1_re;
my $phenotype2_re;
GetOptions(
    "phenotype1_re=s" => \$phenotype1_re,
    "phenotype2_re=s" => \$phenotype2_re,
);
my $outfile_root = shift || die "supply outfile_root\n";
my @headers;
my $header = <>;
chomp $header;
@headers = split /\t/, $header;
my $refseq_idx;
my $pos_idx;
my %sample_idx;
my %phenotype1_samples;
my %phenotype2_samples;
for (my $i=0; $i < @headers; $i++) {
    if ($headers[$i] eq "CHROM") {
        $refseq_idx = $i;
    }
    elsif ($headers[$i] eq "POS") {
        $pos_idx = $i;
    }
    else {
        $sample_idx{$headers[$i]} = $i;
        if (defined $phenotype1_re) {
            if ($headers[$i] =~ /$phenotype1_re/) {
                $phenotype1_samples{$headers[$i]} = 1;
            }
            elsif ($headers[$i] =~ /$phenotype2_re/) {
                $phenotype2_samples{$headers[$i]} = 1;
            }
        }
    }
}
my %markers;
my %samples;
while (<>) {
    next if /^CHROM\tPOS/;
    chomp;
    my @data = split /\t/;
    my ($pos) = ($data[$pos_idx] =~ /(\d+)/);
    my $marker = $data[$refseq_idx] . ":" . $pos;
    $markers{$marker} = {
        refseq => $data[$refseq_idx],
        pos => $pos ,
    };
    $samples{REFERENCE}->{$marker} = REF_ALLELE;
    for my $sample (keys %sample_idx) {
        my $sample_idx = $sample_idx{$sample};
        #FIXME: should parse FORMAT to be sure we have this
        my $gt = $data[$sample_idx];
        if ($gt == 1) {
            $samples{$sample}->{$marker} = REF_ALLELE;
        }
        elsif ($gt == 2) {
            $samples{$sample}->{$marker} = ALT_ALLELE;
        }
        elsif ($gt == 3) {
            $samples{$sample}->{$marker} = HET;
        }
    }
}

open MAP, ">$outfile_root.map";
open GT, ">$outfile_root.gt";
if (defined $phenotype1_re) {
    open GRAPH, ">$outfile_root.graph";
    print GRAPH "MarkerName\tGraphName\tValue\n";
}
foreach my $marker (keys %markers) {
    print MAP "$marker\t$markers{$marker}->{refseq}\t$markers{$marker}->{pos}\n";
    if (defined $phenotype1_re) {
        my $phenotype1_count=0;
        my $phenotype2_count=0;
        foreach my $sample (keys %samples) {
            my $allele = $samples{$sample}->{$marker};
            #FIXME: assumes trait is dominant
            if ($phenotype1_samples{$sample}) {
                if ($allele eq ALT_ALLELE || $allele eq HET) {
                    $phenotype1_count++;
                }
            }
            elsif ($phenotype2_samples{$sample}) {
                if ($allele eq ALT_ALLELE || $allele eq HET) {
                    $phenotype2_count++;
                }
            }
        }
        print GRAPH "$marker\tCOUNT\t",abs($phenotype1_count-$phenotype2_count),"\n";
    }
}
print GT join("\t","",sort keys %markers),"\n";
foreach my $sample (keys %samples) {
    print GT $sample;
    foreach my $marker (sort keys %markers) {
        my $allele = $samples{$sample}->{$marker};
        $allele = "N" unless defined $allele;
        print GT "\t$samples{$sample}->{$marker}";
    }
    print GT "\n";
}
