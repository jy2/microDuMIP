#!/usr/bin/perl

# 585     chr1    62202   62203   C/T     2       0,0     64,64
# 585     chr1    101685  101686  A/G     2       0,0     67,67
# http://moma.ki.au.dk/genome-mirror/cgi-bin/hgTables?db=hg18&hgta_group=varRep&hgta_track=pgSnp&hgta_table=pgNA12878&hgta_doSchema=describe+table+schema

# This program is to select pgSNP in BED file. 
# Developed by Jung-Ki Yoon (neododari@gmail.com)

use strict; 
# example:
# chr7    55086965        55087063        EGFR.1     <-- exon+splicing(5bp)
# chr7    55209973        55210135        EGFR.2

my $synt_prt = "./pgSNPinBED.pl [pgSNP file] [BED] \n";
my @chunks;
my $pgSNP = $ARGV[0];
my $BED_fname = $ARGV[1];
if ($#ARGV != 1) {
    print $synt_prt;
    exit;
}
my %BED;

open(IN, "<$BED_fname") or die "Can't open $BED_fname.\n";
while (my $line=<IN>){
    @chunks = split(/\s+/, $line);
    $BED{$chunks[3]}->{'chr'} = $chunks[0];
    $BED{$chunks[3]}->{'begin'} = $chunks[1];
    $BED{$chunks[3]}->{'end'} = $chunks[2];
}
close(IN);

open(IN, "<$pgSNP") or die "Can't open $pgSNP.\n";
while (my $line=<IN>){
    @chunks = split(/\s+/, $line);
    foreach my $key (keys(%BED)){
	if ($BED{$key}->{'chr'} eq $chunks[1]){
	    if (($BED{$key}->{'begin'}<=$chunks[2])&&($BED{$key}->{'end'}>=$chunks[3])){
		print join("\t", @chunks, $key, $BED{$key}->{'begin'}, $BED{$key}->{'end'}), "\n";
	    }
	}
    }
}
close(IN);
exit;
