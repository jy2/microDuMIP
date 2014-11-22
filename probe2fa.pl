#!/usr/bin/perl
use strict; 

my $synt_prt = "./probearm_seq2fasta.pl CP230.duMIP\n";
my $in_fname = $ARGV[0];
my $len_armID = 21; #min_length of arm seq

if (($#ARGV != 0)) {
    print $synt_prt;
    exit;
}
# #num    gene.exon.MIP   chr     begin   end     MIP_type        cover_begin     cover_end       gap_begin       gap_end seq_EXT Tm_EXT GC_EXT LEN_EXT seq_LIG Tm_LIG GC_LIG LEN_LIG
# 0       ABL1.0.0        chr9    133589706       133589847       normal          133589706       133589809       133589698       133589809     TTCTGGAAAGGGGTACCTATT  58.49  42 4_21 GGGGTCCACACTGCAATGTTTT 64.55 50 0_22

open(IN, "<$in_fname") or die "Can't open $in_fname.\n";
while(<IN>){
    if ($_ =~ "#"){ next;}
    my @chunks = split(/\s+/, $_);
    print ">$chunks[1]\n", substr($chunks[10], 0, $len_armID), substr(rev_comp($chunks[14]), 0, $len_armID), "\n";
}
close(IN);
exit;

sub rev_comp
{
    my $in_ = uc(reverse(shift @_));
    my $re_="";
    my $cur_;
    for (my $i_=0; $i_<length($in_); $i_++){
	$cur_ = substr($in_, $i_, 1);
	if ($cur_ eq "A"){
	    $re_ .= "T";
	}elsif($cur_ eq "T"){
	    $re_ .= "A";
	}elsif($cur_ eq "C"){
	    $re_ .= "G";
	}elsif($cur_ eq "G"){
	    $re_ .= "C";
	}elsif($cur_ eq "N"){
	    $re_ .= "N";
	}else{
	    print "Err : can't distinguish $cur_ in $in_\n";
	    exit;
	}
    }
    return $re_;
}
