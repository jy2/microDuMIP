#!/usr/bin/perl
use strict;

# This program is to extract arm sequence from reads and make fasta format file for novoalign.
# Developed by JK Yoon : 2013. 9. 12.   (for duMIP)
# Modified; A tail is now considered: 2013. 10. 11
# Modified; truncated read is now considered: 2013. 10. 13. 

my $len_armID = 21;

my $synt_str = "./armseq2fa.pl fwd(input) rev(input) sliding(0 or 1) output(fasta)\n";

my $fwd_fname = $ARGV[0];
my $rev_fname = $ARGV[1];
my $sliding = $ARGV[2];
my $out_fname = $ARGV[3];

if ($#ARGV != 3){
    print $synt_str;
    exit;
}

# duMIP
# Forward(1) : barcode(6mer)-T-ampF(21mer)-ext_arm(21~25mer)-targets
# Reverse(2) : T-ampR(19mer)-randomN(10mer)-lig_arm(21~25mer)-targets
# if Sliding != 0, the first bases before ampF and ampR were additional spacers(1bp)
# INPUT  : @HISEQ20:84:H0949ADXX:2:1101:5190:2171 1:N:0:GATCAG

my @headers_F;
my @headers_R;
my $line_F;
my $line_R;

#@HISEQ20:84:H0949ADXX:2:1101:5190:2171 1:N:0:GATCAG
#TCAAGTGCAGATGTTATCGAGGTCCGACATGGACCTCATCGTGCTTTCTGTTGCACATACCTGCATATGCCATGATCGCATGAGCCATGCCCACTGCACATGTTGGGGCACTGCTGCCCAATGTAAATGCTGTCCAAAGCCCACTCGTCAA
#+
#@@?DFFFFHHGHHIHJJJJJJIJJJIJJJJIGIJIIJJIIGGGHGIJCGHJJIIJJJJJIHJIIGHHHHGFFFFFFEDDDBBCCCDDDCDDDDDBCDCCCDCDEDDBDDDDDA@A@CBDDDDDACBDDDCCC@CCDDDCCDDDBBABDDD5

# AmpF: CAGATGTTATCGAGGTCCGAC
# AmpR: GGAACGATGAGCCTCCAAC

open(INF, "<$fwd_fname") or die "Can't open $fwd_fname\n";
open(INR, "<$rev_fname") or die "Can't open $rev_fname\n";
open(OUT, ">$out_fname");
my $num_line=0;
while($line_R = <INR>){
    $num_line++;
    if ($num_line % 10000 == 0){
	print "$num_line were processed.\n";
    }
    $line_F = <INF>;
    @headers_F = split(/\s+/, $line_F);
    @headers_R = split(/\s+/, $line_R);
    if ((substr($headers_F[1], 0, 1) != 1)||(substr($headers_R[1], 0, 1) != 2)){
	print "Err : wrong direction code!\n$synt_str\nFWD:$line_F\nREV:$line_R";
	exit;
    }
    if ($headers_F[0] ne $headers_R[0]){
	print "Err : headers are not same\n FWD:$line_F\nREV:$line_R";
	exit;
    }
    substr($headers_F[0], 0, 1) = ""; # remove the first letter @

    $line_F = <INF>;
    $line_R = <INR>;
    chop($line_F);
    chop($line_R);

    if($sliding>0) {
	substr($line_F, 0, 1)="";
	substr($line_R, 0, 1)="";
    }

    substr($line_F, 0, 6+1) = ""; # remove barcode and "T"
    substr($line_R, 0, 1) = ""; 

    my $seq_F = substr($line_F, 0, 30);
    my $seq_R = substr($line_R, 0, 30);
    my @temp;

    my $bF = my $bR = 0;

    if (($seq_F =~ "GGTCCGAC")&&($seq_R =~ "CCTCCAAC")){
	@temp = split("GGTCCGAC", $seq_F);
	substr($line_F, 0, 30) = "";
	$line_F = $temp[1] . $line_F;
	$bF = 6 + 1 + length($temp[0]) + 8;
	@temp = split("CCTCCAAC", $seq_R);
	substr($line_R, 0, 30) = "";
	$line_R = $temp[1] . $line_R;
	$bR = 1 + length($temp[0]) + 8;
	print OUT ">>_$bF\_$bR\_$headers_F[0]\n";
	print OUT substr($line_F, 0, $len_armID), substr($line_R, 10, $len_armID), "\n";
    }elsif(($seq_R =~ "GGTCCGAC")&&($seq_F =~ "CCTCCAAC")){
	# adapter ligated on opposite direction.
	@temp = split("GGTCCGAC", $seq_R);
	substr($line_R, 0, 30) = "";
	$line_R = $temp[1] . $line_R;
	$bR = 1 + length($temp[0]) + 8;
	@temp = split("CCTCCAAC", $seq_F);
	substr($line_F, 0, 30) = "";
	$line_F = $temp[1] . $line_F;
	$bF = 6 + 1 + length($temp[0]) + 8;
	print OUT "><_$bF\_$bR\_$headers_F[0]\n";
	print OUT substr($line_R, 0, $len_armID), substr($line_F, 10, $len_armID), "\n";
    }else{
	# Cannot distinguish arm seq....
	print OUT ">-_0_0_$headers_F[0]\n";
	print OUT substr($line_F, 28, $len_armID), substr($line_R, 31, $len_armID), "\n";
    }

    $line_F = <INF>; 
    $line_R = <INR>; 
    $line_F = <INF>; 
    $line_R = <INR>; 
}
close(INF);
close(INR);
close(OUT);
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
