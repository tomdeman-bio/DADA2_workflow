#! /anaconda2/bin/perl

use strict;

#RDP default output + species info
my $input = shift;
my $fh;

open FILE, "$input" or die "cannot open $input for reading \n";
open ($fh, '>', "$input.c80_BIOMversion.tsv");
print $fh "#OTU ID\ttaxonomy\n";

while (<FILE>) {
    chomp;
    my @out = split ("\t", $_);
    #check bacteria
    if ($out[7] < 0.8) {
        print $fh "$out[0]\tUnclassified;Unclassified;Unclassified;Unclassified;Unclassified;Unclassified;Unclassified\n";
        next;
    }
    #check phylum
    elsif ($out[10] < 0.8) {
        print $fh "$out[0]\t$out[5];Unclassified;Unclassified;Unclassified;Unclassified;Unclassified;Unclassified\n";
        next;
    }
    #check class
    elsif ($out[13] < 0.8) {
        print $fh "$out[0]\t$out[5];$out[8];Unclassified;Unclassified;Unclassified;Unclassified;Unclassified;\n";
        next;
    }
    #check order
    elsif ($out[16] < 0.8) {
        print $fh "$out[0]\t$out[5];$out[8];$out[11];Unclassified;Unclassified;Unclassified;Unclassified\n";
        next;
    }
    #check family
    elsif ($out[19] < 0.8) {
        print $fh "$out[0]\t$out[5];$out[8];$out[11];$out[14];Unclassified;Unclassified;Unclassified\n";
        next;
    }
    #check genus
    elsif ($out[22] < 0.8) {
        print $fh "$out[0]\t$out[5];$out[8];$out[11];$out[14];$out[17];Unclassified;Unclassified\n";
        next;
    }
    #check species
    elsif ($out[25] < 0.8) {
        print $fh "$out[0]\t$out[5];$out[8];$out[11];$out[14];$out[17];$out[20];$out[23]<0.8\n";
        next;
    }
    elsif ($out[25] >= 0.8) {
        print $fh "$out[0]\t$out[5];$out[8];$out[11];$out[14];$out[17];$out[20];$out[23]\n";
        next;
    }
}
