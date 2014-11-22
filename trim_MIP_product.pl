#!/usr/bin/perl
use strict;

# This program is to trim NGS reads of MIP capture. It will identify MIP_ID and target seq from arm sequences. 
# Developed by JK Yoon
# trim_JK_MIPcode.pl: : 2013. 3. 22.
# -> trim_MIP_product.pl : 2013. 8. 1.   (for duMIP)
# -> trim_MIP_product_w_novo.pl : 2013. 9. 12. (using novoalign to identify arm seq)
# -> updated script for modified ligation method; considered A tails.: 2013. 10. 11. 
# -> modified for char. of new ligation method; some read might be truncated... : 2013. 10. 13.

my $target_size = 116;
my $len_armID = 21;
#my $max_anneal = 25;
#my $check_size = 12;
my $len_ranN = 10;

#my $synt_str = "./trim_MIP_product.pl CP230.duMIP CP230.duMIP.stat(output) fwd(input) rev(input) fwd(output) rev(output) sliding(0 or 1)\n";
my $synt_str = "./trim_MIP_product_w_novo.pl CP230.duMIP CP230.duMIP.stat(output) fwd(input) rev(input) arminfo_novo(input) fwd(output) rev(output) sliding(0=no or 1=yes)\n";

my $MIP_fname = $ARGV[0];
my $MIP_out_fname = $ARGV[1];
my $fwd_fname = $ARGV[2];
my $rev_fname = $ARGV[3];
my $arminfo_fname = $ARGV[4];
#my $out_fwd_fname = $ARGV[2] . ".trim";
#my $out_rev_fname = $ARGV[3] . ".trim";
#my $sliding = $ARGV[4];
my $out_fwd_fname = $ARGV[5];
my $out_rev_fname = $ARGV[6];
my $sliding = $ARGV[7];

my $remain_fwd_fname = $ARGV[5] . ".dropped";
my $remain_rev_fname = $ARGV[6] . ".dropped";

my $num_MIP=0;
my $num_succ=0;

if (($#ARGV != 7) || (($sliding!=0)&&($sliding!=1))){
    print $synt_str;
    exit;
}
#if ($sliding > 0){ $sliding = 1;}

# duMIP
# Forward(1) : barcode(6mer)-T-ampF(21mer)-ext_arm(21~25mer)-targets
# Reverse(2) : T-ampR(19mer)-randomN(10mer)-lig_arm(21~25mer)-targets
# if Sliding != 0, the first bases before ampF and ampR were additional spacers(1bp)
# INPUT  : @HISEQ20:84:H0949ADXX:2:1101:5190:2171 1:N:0:GATCAG
# OUTPUT : @randomN_armID:84:H0949ADXX:2:1101:5190:2171 1:N:0:GATCAG

# AmpF: CAGATGTTATCGAGGTCCGAC
# AmpR: GGAACGATGAGCCTCCAAC

my %MIP;

#@HISEQ20:84:H0949ADXX:2:1101:5190:2171 1:N:0:GATCAG
#TCAAGTGCAGATGTTATCGAGGTCCGACATGGACCTCATCGTGCTTTCTGTTGCACATACCTGCATATGCCATGATCGCATGAGCCATGCCCACTGCACATGTTGGGGCACTGCTGCCCAATGTAAATGCTGTCCAAAGCCCACTCGTCAA
#+
#@@?DFFFFHHGHHIHJJJJJJIJJJIJJJJIGIJIIJJIIGGGHGIJCGHJJIIJJJJJIHJIIGHHHHGFFFFFFEDDDBBCCCDDDCDDDDDBCDCCCDCDEDDBDDDDDA@A@CBDDDDDACBDDDCCC@CCDDDCCDDDBBABDDD5

##num    gene.exon.MIP   chr     begin   end     MIP_type        cover_begin     cover_end       gap_begin       gap_end seq_EXT Tm_EXT GC_EXT LEN_EXT seq_LIG Tm_LIG GC_LIG LEN_LIG
#0       ABL1.0.0        chr9    133589706       133589847       normal          133589706       133589809       133589698       133589809     TTCTGGAAAGGGGTACCTATT  58.49  42 4_21 GGGGTCCACACTGCAATGTTTT 64.55  50 0_22

open(IN, "<$MIP_fname") or die "Can't open $MIP_fname\n";
while(<IN>){
    if (substr($_, 0, 1) eq "#"){next;}
    my @chunks = split(/\s+/, $_);
#    $MIP{$chunks[1]}->{'ext'} = substr($chunks[10], 0, $len_armID);  
#    $MIP{$chunks[1]}->{'lig'} = substr($chunks[14], 0, $len_armID);
#    $MIP{$chunks[1]}->{'line'} = $_;
    my $t;
    ($t, $MIP{$chunks[1]}->{'e_len'}) = split(/_/, $chunks[13]);
    ($t, $MIP{$chunks[1]}->{'l_len'}) = split(/_/, $chunks[17]);
    $MIP{$chunks[1]}->{'gap_len'} = $chunks[9] - $chunks[8];
    $MIP{$chunks[1]}->{'count'} = 0;
    $MIP{$chunks[1]}->{'perfect'} = 0;
    $num_MIP++;
}
close(IN);
print "$num_MIP MIPs loaded\n";

open(INF, "<$fwd_fname") or die "Can't open $fwd_fname\n";
open(INR, "<$rev_fname") or die "Can't open $rev_fname\n";
open(ARM, "<$arminfo_fname") or die "Can't open $arminfo_fname\n";
# @HD     VN:1.0  SO:unsorted
# @PG     ID:novoalign    PN:novoalign    VN:V2.07.18     CL:novoalign -d CP230.ndx -f both100_001.fa -o SAM
# >_15_23_HWI-ST1303:146:H12K4ADXX:1:1101:19927:12255    0       NOTCH1.3.3      1       70      42M     *       0       0       TGAGTAGCGGGCGGCCAGGTGTCCTGCAGCGGGGGCGGCCTG      *       PG:Z:novoalign  AS:i:0  UQ:i:0  NM:i:0  MD:Z:42
my $num_line=0;
my @line_F; 
my @line_R;
open(OUTF, ">$out_fwd_fname");
open(OUTR, ">$out_rev_fname");
open(REMF, ">$remain_fwd_fname");
open(REMR, ">$remain_rev_fname");
my $num_succ_orient=0;

while(my $line_A = <ARM>){
    if (substr($line_A, 0, 1) eq "@") {next;}

    $num_line++;
    if ($num_line % 10000 == 0){
	print "$num_line were processed. $num_succ were saved.\n";
    }
    my @arm = split(/\s+/, $line_A);
#    my $code = substr($arm[0], 0, 1); # >: proper orientation, < : opposite orientation, - : cannot distinguish
#    substr($arm[0], 0, 1) = "";
    (my $code, my $bF, my $bR, my $r_name,) = split(/_/, $arm[0]);
    if (($code eq ">")||($code eq "<")){
	$num_succ_orient++;
    }

#    $arm[0] = "@" . $arm[0];
    $r_name = "@" . $r_name;
    $line_F[0] = <INF>;
    $line_F[1] = <INF>;
    $line_F[2] = <INF>;
    $line_F[3] = <INF>;
    $line_R[0] = <INR>;
    $line_R[1] = <INR>;
    $line_R[2] = <INR>;
    $line_R[3] = <INR>;

    my @headers_F = split(/\s+/, $line_F[0]);
    my @headers_R = split(/\s+/, $line_R[0]);
    if ((substr($headers_F[1], 0, 1) != 1)||(substr($headers_R[1], 0, 1) != 2)){
	print "Err : wrong direction code!\n$synt_str\nFWD:$line_F[0]\nREV:$line_R[0]";
	exit;
    }
    if (!(($headers_F[0] eq $headers_R[0])&&($headers_F[0] eq $r_name))){
#    if (!(($headers_F[0] eq $headers_R[0])&&($headers_F[0] eq $arm[0]))){
	print "Err : headers are not same\n FWD:$line_F[0]\nREV:$line_R[0]\nARMINFO:$line_A";
	print "$headers_F[0]\n$headers_R[0]\n$r_name\n";
	exit;
    }

    if (($arm[2] ne "*")&($code ne "-")){ # arm seq. and orientation are identified!
	$num_succ++;
	$MIP{$arm[2]}->{'count'}++;
	if ($arm[5] eq "42M") {$MIP{$arm[2]}->{'perfect'}++;}
	chop($line_F[1]);
	chop($line_R[1]);
	chop($line_F[3]);
	chop($line_R[3]);
	substr($line_F[1], 0, $sliding+$bF) = ""; # trim barcode and A tail
	substr($line_F[3], 0, $sliding+$bF) = "";
	substr($line_R[1], 0, $sliding+$bR) = ""; # trim A tail
	substr($line_R[3], 0, $sliding+$bR) = "";

	if ($code eq "<"){
	    my $temp = $line_F[1];
	    $line_F[1] = $line_R[1];
	    $line_R[1] = $temp;
	    $temp = $line_F[3];
	    $line_F[3] = $line_R[3];
	    $line_R[3] = $temp;
	}

	my $cur_begin = $MIP{$arm[2]}->{'e_len'}; # retrieve ext_arm
	$line_F[1] = substr($line_F[1], $cur_begin, $MIP{$arm[2]}->{'gap_len'}) . "\n";
	$line_F[3] = substr($line_F[3], $cur_begin, $MIP{$arm[2]}->{'gap_len'}) . "\n";

	$cur_begin = 0; # retrieve random N 10mer
	my $cur_ranN = substr($line_R[1], $cur_begin, $len_ranN);

	$cur_begin = 10 + $MIP{$arm[2]}->{'l_len'}; # retrieve lig_arm
	$line_R[1] = substr($line_R[1], $cur_begin, $MIP{$arm[2]}->{'gap_len'}) . "\n";
	$line_R[3] = substr($line_R[3], $cur_begin, $MIP{$arm[2]}->{'gap_len'}) . "\n";
	my @temp = split(/:/, $headers_F[0]);
	shift(@temp);
	$headers_F[0] = $headers_R[0] = join(":", ("@" . $cur_ranN . "_" . $arm[2]), @temp); 
	$line_F[0] = "$headers_F[0] $headers_F[1]\n";
	$line_R[0] = "$headers_R[0] $headers_R[1]\n";
	print OUTF @line_F;
	print OUTR @line_R;
    } else{
	print REMF @line_F;
	print REMR @line_R;
    }
}
close(INF);
close(INR);
close(ARM);
close(OUTF);
close(OUTR);
close(REMF);
close(REMR);
(my $codename, )=split(/\./, $arminfo_fname); 
print ">> $codename: Among $num_line -> $num_succ_orient has orientation ->  $num_succ were passed.\nRemaining reads data were stored in $remain_fwd_fname and $remain_rev_fname.\n";

open(OUT, ">$MIP_out_fname") or die "Can't open $MIP_out_fname\n";
print OUT join("\t", "MIP_ID", "count", "perfect_match", ("len_armID=" . $len_armID*2)), "\n";
foreach my $key (keys %MIP){
    print OUT join("\t", $key, $MIP{$key}->{'count'}, $MIP{$key}->{'perfect'}), "\n";
}
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
