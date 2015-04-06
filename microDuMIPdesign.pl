#!/usr/bin/perl

# This program is to design duplex molecular inversion probe(duMIP) to capture givin genomic regions (BED format)
# Developed by Jung-Ki Yoon (dr.jkyoon@gmail.com)

# Note : all starting points in this program is 0-based. http://genome.ucsc.edu/FAQ/FAQtracks#tracks1
# The real starting point is 1 base next to the position. 
 
# Version log
# ver1 2012/11/28
# ver1.01 2012/12/17 changed get_anneal_seq()
# ver1.1 2012/12/18  changed get_anneal_seq(); search anneal seq which having Tm closest to opt_Tm, anneal seq are excluded target region. 
# ver1.2 changed get_anneal_seq(); check uniqueness on genome of anneal seq. 
# ver1.3 2013/06/04 modified for duMIP, add EarI seq and flanking seq. Search the best arm seq within 5bps, splicing is removed (user must add splicing sites in input BED file)


# Before running this script, run gfServer first. 
# ./gfServer start (server IP) (port ex:6666) -stepSize=5 -log=untrans.log hg19.2bit
# ./gfClient (server IP) (port) ~/JK/blat/i386/ test2.fa test2_server.psl -minScore=0 -minIdentity=0

# $seq_path = "/home/schwarzwald/reference/hg19/";     reference genomic sequeunces "chr*.fa" 

use strict; 
use Cornerstone;
my $NNParam = "./param";  
#my $NNParam = "./param_Sugimoto";

Cornerstone::readNNParam($NNParam);

# input file format (BED format)
# example:
# chr7    55086965        55087063        EGFR.1     <-- exon+splicing(5bp)
# chr7    55209973        55210135        EGFR.2

my $synt_prt = "./MIP_design.pl [gene_BED] [hg19 path]\n";
my $in_fname = $ARGV[0];
my $seq_path = $ARGV[1];
my $out_fname = "$ARGV[0].MIP";
my $log_fname = "$ARGV[0].log";
if ($#ARGV != 1) {
    print $synt_prt;
    exit;
}
my $temp_fa = "$ARGV[0].temp_fa";
my $temp_psl = "$ARGV[0].temp_psl";
my $temp_log = "$ARGV[0].temp_log";
my $sys_str = "./gfClient 165.132.29.249 6666 ../blat/i386/ -minScore=0 -minIdentity=0 $temp_fa $temp_psl > $temp_log";


open(LOG, ">$log_fname");

my $ver_str = "duMIP_design.pl ver1.3";
my $filling_size = 100;
my $overlapping_size = 50; 
#my $splicing = 5; # definition of splicing region +- 5bp 
my $min_anneal = 21; # min. length of annealing arm
my $max_anneal = 25; # max. length of annealing arm
my $len_sliding = 5; # len. of arm sliding
my $opt_Tm = 60; # Optimal Tm of annealing arm
my $min_Tm = 58; # min. Tm of annealing arm 
my $max_Tm = 65; # max. Tm of annealing arm
my $opt_GC = 50;
my $min_GC = 10;
my $max_GC = 90;
my $homopolymer = 8; # > $homopolymer will removed
my @ATCG = qw(A T C G);
for (my $i=0; $i<4; $i++){
    for (my $j=0; $j<$homopolymer; $j++){
	$ATCG[$i] .= substr($ATCG[$i], 0, 1);
    }
}
#my $MIP_common_seq = "NNNNNNNNNN" . "GTTGGAGGCTCATCGTTCC" . "TATTCAGG" . "CAGATGTTATCGAGGTCCGAC"; 
my $MIP_common_seq = "NNNNNNNNNN" . "GTTGGAGGCTCATCGTTCC" . "CAGATGTTATCGAGGTCCGAC"; 
my @flanking_L;
my @flanking_E; 
$flanking_L[0] = "GGTAGCAAAGTGCAGATGTG";
$flanking_E[0] = "TTCAGAGCAGTGTGAGTTCA";
$flanking_L[1] = "CTATGAGCATGTTCTTCAGG";
$flanking_E[1] = "CAGCAAGCGTAATTAACTGC";
my $EarI_L = "CTCTTC" . "A";
my $EarI_E = "" . "GAAGAG";

# flanking_L(20mer) + EarI_L(6+1mer) + ligation_arm(25mer) + barcode(10mer) + AmpR(19mer) + AmpF(21mer) + extension_arm(25mer) + EarI_E(6+4mer) + flanking_E(20mer)
# = 157mer
# AmpR 61.17    AmpF 60.68




#print "$ver_str\nInput file name = $in_fname\nCapture size = $filling_size\nSplicing regions = $splicing\nTm = $opt_Tm($min_Tm ~ $max_Tm)\nGC = $opt_GC($min_GC ~ $max_GC)\nMIP common seq = $MIP_common_seq\n";
#print LOG "$ver_str\nInput file name = $in_fname\nCapture size = $filling_size\nTm = $opt_Tm($min_Tm ~ $max_Tm)\nGC = $opt_GC($min_GC ~ $max_GC)\nMIP common seq = $MIP_common_seq\n";
print "$ver_str\nInput file name = $in_fname\nCapture size = $filling_size\nOverlapping size = $overlapping_size\nSliding size = $len_sliding\nArm leng. = $min_anneal ~ $max_anneal\nOpt.Tm = $opt_Tm ($min_Tm ~ $max_Tm)\nhomopolymer <= $homopolymer\nGC% = $min_GC ~ $max_GC\nMIP common seq = $MIP_common_seq\n\n";
print LOG "$ver_str\nInput file name = $in_fname\nCapture size = $filling_size\nOverlapping size = $overlapping_size\nSliding size = $len_sliding\nArm leng. = $min_anneal ~ $max_anneal\nOpt.Tm = $opt_Tm ($min_Tm ~ $max_Tm)\nhomopolymer <= $homopolymer\nGC% = $min_GC ~ $max_GC\nMIP common seq = $MIP_common_seq\n";

my $in_chr; # input data
my $loaded_chr="blank";
my $in_begin;
my $in_end;

my $cur_length; # capture region in exon
my $cur_begin;
my $cur_end;

my $MIP_begin; #MIP target region = gap-filling + annealing regions + slidings
my $MIP_end;
my $MIP_length; #cover length of each MIP : annealing parts will be removed before analyzed. 
$MIP_length = $filling_size+2*$max_anneal;
my $gap_begin;
my $gap_end;
my $cover_begin;
my $cover_end;

my $cur_name; # gene.exon name
my $cur_seq = ""; # genomic seq of input chromosome 
my $cur_MIP; # MIP sequence
my $num_MIP4exon; # number of probe for each gene.exon
my $num_MIP=0;

open(OUT, ">$out_fname");
print "Input file : $in_fname\nOutput file : $out_fname\nLog file : $log_fname\n";
open(IN, "<$in_fname") or die "Can't open $in_fname.\n";
my $cur_target; 
print LOG join ("\t", "gene.exon", "chr", "begin", "end", "target_length", "target_seq"), "\n";
print LOG join ("\t", "MIP_num", "gene.exon.MIP", "remaining_seq_len", "MIP_type", "cover_begin", "cover_end", "gap_begin", "gap_end", "seq_EXT", "Tm_EXT", "GC_EXT", "LEN_EXT", "seq_LIG", "Tm_LIG", "GC_LIG", "LEN_LIG", "target_seq"), "\n\n";
print OUT join ("\t", "#num", "gene.exon.MIP", "chr", "begin", "end", "MIP_type", "cover_begin", "cover_end", "gap_begin", "gap_end", "seq_EXT", "Tm_EXT", "GC_EXT", "LEN_EXT", "seq_LIG", "Tm_LIG", "GC_LIG", "LEN_LIG", "seq_preMIP"), "\n";

my $MIP_type;
my $num_lines==0;

while (my $line=<IN>){
    $num_lines++;
    if ($num_lines % 100 == 0){
	print "$num_lines were processed. $num_MIP probes were made.\n";
    }
    ($in_chr, $in_begin, $in_end, $cur_name, ) = split(/\s+/, $line);
    if ($in_chr ne $loaded_chr) {
	$cur_seq = load_seq($in_chr, $seq_path);
    }
    if ($in_begin>=$in_end){
	print "Err: begin >= end; $line";
	exit;
    }
    $num_MIP4exon = 0;
    $cur_begin = $in_begin;  
    $cur_length = $in_end-$in_begin;
    $cur_end = $in_end;
    print LOG join ("\t", $cur_name, $in_chr, $in_begin, $in_end, $cur_length, substr($cur_seq, $cur_begin, $cur_length)), "\n";


    if ($cur_length < $filling_size){   # if the exon is short
	$MIP_type = "short";
	$MIP_begin = $cur_begin - int(abs($filling_size - $cur_length)/2+0.5) - $max_anneal;
	$MIP_end = $MIP_begin + $MIP_length;
	get_anneal_seq();
    } else{
	while ($cur_length > 0){
	    if ($cur_length >= $filling_size){
		$MIP_type = "normal";
		$MIP_begin = $cur_begin - $max_anneal;
		$MIP_end = $MIP_begin + $MIP_length;
		get_anneal_seq();
		$cur_length -= $overlapping_size;
		$cur_begin += $overlapping_size;
	    } else{
		$MIP_type = "last";
		$MIP_begin = $cur_begin - $max_anneal;
		$MIP_end = $MIP_begin + $MIP_length; # extend to intron
		get_anneal_seq();
		$cur_length = 0;
	    }
	    $num_MIP4exon++;
	}
    }
    print LOG "\n";
}
close(IN);
close(LOG);
close(OUT);
print "$num_lines were completed.\n";
exit;

sub get_anneal_seq  # Make arm sequences
{
    my $in_target = shift @_; # calc. Tm, GC of all combinations. 
    my %EXT;
    my %LIG;
    my $code;
    open(TEMP, ">$temp_fa");
    for (my $i=0; $i<$len_sliding; $i++){
	for (my $j=$min_anneal; $j<=$max_anneal; $j++){
	    $code = $i . "_" . $j;
	    $EXT{$code}->{'seq'} = substr($cur_seq, $MIP_begin-$i, $j);
	    $LIG{$code}->{'seq'} = substr($cur_seq, $MIP_end+$i-$j, $j);
	    ($EXT{$code}->{'status'}, $LIG{$code}->{'status'}) = check_EarI($EXT{$code}->{'seq'}, $LIG{$code}->{'seq'});
	    $EXT{$code}->{'status'} = check_homopolymer($EXT{$code}->{'seq'}, $EXT{$code}->{'status'});
	    $LIG{$code}->{'status'} = check_homopolymer($LIG{$code}->{'seq'}, $LIG{$code}->{'status'});
	    
	    $EXT{$code}->{'redun'} = 0;
	    $LIG{$code}->{'redun'} = 0;
	    $EXT{$code}->{'Tm'} = Cornerstone::shortOligoTm($EXT{$code}->{'seq'});
	    $LIG{$code}->{'Tm'} = Cornerstone::shortOligoTm($LIG{$code}->{'seq'});
	    $EXT{$code}->{'GC'} = Cornerstone::shortOligoGC($EXT{$code}->{'seq'});
	    $LIG{$code}->{'GC'} = Cornerstone::shortOligoGC($LIG{$code}->{'seq'});
	    $EXT{$code}->{'status'} = check_GC($EXT{$code}->{'GC'}, $EXT{$code}->{'status'});
	    $LIG{$code}->{'status'} = check_GC($LIG{$code}->{'GC'}, $LIG{$code}->{'status'});
	    print TEMP ">EXT_$code\n", $EXT{$code}->{'seq'}, "\n>LIG_$code\n", $LIG{$code}->{'seq'}, "\n";
	}
    }
    close(TEMP);
    system($sys_str);
    open(INTEMP, "<$temp_psl");
    my @temp_names;
    while(my $temp_line = <INTEMP>){
	my @temp_chunks = split(/\s+/, $temp_line);
	if (($temp_chunks[0] == $temp_chunks[10])&&($temp_chunks[0] == $temp_chunks[12])){ #perfect match
	    @temp_names = split(/_/, $temp_chunks[9]);
	    $code = $temp_names[1] . "_" . $temp_names[2]; 
	    if ($temp_names[0] eq "EXT"){
		$EXT{$code}->{'redun'}++;
		if ($EXT{$code}->{'redun'}>1){
		    $EXT{$code}->{'status'} = "redun";
		}
	    } elsif ($temp_names[0] eq "LIG"){
		$LIG{$code}->{'redun'}++;
		if ($LIG{$code}->{'redun'}>1){
		    $LIG{$code}->{'status'} = "redun";
		}
	    }
	}
    }
    close(INTEMP);
    
    my $best_Tm = 1000;
    my $best_GC = 1000;
    my $idx = "Tm";
    my $EXT_key = -1;
    my $LIG_key = -1;
    foreach my $key (sort(keys(%EXT))){
	if ($EXT{$key}->{'status'} eq "pass"){
	    if (($EXT{$key}->{'Tm'}>=$min_Tm)&&($EXT{$key}->{'Tm'}<=$max_Tm)){
		$idx = "GC";
		if (abs($EXT{$key}->{'GC'}-$opt_GC)<abs($best_GC-$opt_GC)){
		    $best_GC = $EXT{$key}->{'GC'};
		    $EXT_key = $key;
		}
	    } else{
		if ($idx eq "GC"){
		    #ignore
		} elsif ($idx eq "Tm"){
		    if (abs($EXT{$key}->{'Tm'}-$opt_Tm)<abs($best_Tm-$opt_Tm)){
			$best_Tm = $EXT{$key}->{'Tm'};
			$EXT_key = $key;
		    }
		}
	    }
	}
	print OUT join("\t", "EXT_$key", $EXT{$key}->{'seq'}, $EXT{$key}->{'Tm'}, $EXT{$key}->{'GC'}, $EXT{$key}->{'redun'}, $EXT{$key}->{'status'}, $idx, $EXT_key), "\n";
    }

    $best_Tm = 1000;
    $best_GC = 1000;
    $idx = "Tm";
    foreach my $key (sort(keys(%LIG))){
	if ($LIG{$key}->{'status'} eq "pass"){
	    if ($min_Tm<=($LIG{$key}->{'Tm'})&&($LIG{$key}->{'Tm'}<=$max_Tm)){
		$idx = "GC";
		if (abs($LIG{$key}->{'GC'}-$opt_GC)<abs($best_GC-$opt_GC)){
		    $best_GC = $LIG{$key}->{'GC'};
		    $LIG_key = $key;
		}
	    } else{
		if ($idx eq "GC"){
		    #ignore
		} elsif ($idx eq "Tm"){
		    if (abs($LIG{$key}->{'Tm'}-$opt_Tm)<abs($best_Tm-$opt_Tm)){
			$best_Tm = $LIG{$key}->{'Tm'};
			$LIG_key = $key;
		    }
		}
	    }
	}
	print OUT join("\t", "LIG_$key", $LIG{$key}->{'seq'}, $LIG{$key}->{'Tm'}, $LIG{$key}->{'GC'}, $LIG{$key}->{'redun'}, $LIG{$key}->{'status'}, $idx, $LIG_key), "\n";
    }

    my $cur_flk;
    my $spacer_EarI;
    my $code_i;
    my $code_j;

    if (($EXT_key == -1)||($LIG_key == -1)){
	# fail to design
	$MIP_type = "fail:EXT$EXT_key,LIG:$LIG_key";
	print OUT join("\t", $num_MIP, ($cur_name . "." . $num_MIP4exon), $in_chr, $in_begin, $in_end, $MIP_type), "\n";
	print LOG join("\t", $num_MIP, ($cur_name . "." . $num_MIP4exon), $cur_length, $MIP_type, $cur_target), "\n";
    }else{

	($code_i, $code_j) = split(/_/, $EXT_key);
	$cur_flk = $num_MIP % 2;
	$spacer_EarI = substr($cur_seq, $MIP_begin - $code_i + $code_j, 4);

	$gap_begin = $MIP_begin - $code_i + $code_j;
	($code_i, $code_j) = split(/_/, $LIG_key);
	$gap_end = $MIP_end + $code_i - $code_j;

	if ($gap_begin < $in_begin){
	    $cover_begin = $in_begin;
	}else{
	    $cover_begin = $gap_begin;
	}
	if ($gap_end > $in_end){
	    $cover_end = $in_end;
	}else{
	    $cover_end = $gap_end;
	}

	# build precursor of duMIP
	$cur_MIP = $flanking_L[$cur_flk] . $EarI_L . $LIG{$LIG_key}->{'seq'} . $MIP_common_seq . $EXT{$EXT_key}->{'seq'} . $spacer_EarI . $EarI_E . $flanking_E[$cur_flk];
	print OUT join("\t", $num_MIP, ($cur_name . "." . $num_MIP4exon), $in_chr, $in_begin, $in_end, $MIP_type, $cover_begin, $cover_end, $gap_begin, $gap_end, $EXT{$EXT_key}->{'seq'}, $EXT{$EXT_key}->{'Tm'}, $EXT{$EXT_key}->{'GC'}, $EXT_key, $LIG{$LIG_key}->{'seq'}, $LIG{$LIG_key}->{'Tm'}, $LIG{$LIG_key}->{'GC'}, $LIG_key, $cur_MIP), "\n";
	print LOG join("\t", $num_MIP, ($cur_name . "." . $num_MIP4exon), $cur_length, $MIP_type, $cover_begin, $cover_end, $gap_begin, $gap_end, $EXT{$EXT_key}->{'seq'}, $EXT{$EXT_key}->{'Tm'}, $EXT{$EXT_key}->{'GC'}, $EXT_key, $LIG{$LIG_key}->{'seq'}, $LIG{$LIG_key}->{'Tm'}, $LIG{$LIG_key}->{'GC'}, $LIG_key, $cur_target), "\n";
    }

    $num_MIP++;
    print "\n";
    return;
}
sub check_GC
{
    my $GC_ = shift @_;
    my $re_ = shift @_;
    if (($GC_ < $min_GC)||($GC_ > $max_GC)){
	$re_ = "badGC";
    }
    return $re_;
}
sub check_homopolymer
{
    my $str_ = shift @_;
    my $re_ = shift @_;
    for(my $c_=0; $c_<4; $c_++){
	if ($str_ =~ $ATCG[$c_]){
	    $re_ = "homo";
	}
    }
    return $re_;
}
sub check_EarI
{
    my $ext_ = shift @_;
    my $lig_ = shift @_;
    my $re_ext = "pass";
    my $re_lig = "pass";
#my $EarI_L = "CTCTTC" . "A";
#my $EarI_E = "" . "GAGAAG";;

    $ext_ = "C" . $ext_; #from ampF
    if (($ext_ =~ "CTCTTC")||($ext_ =~ "GAAGAG")){
	$re_ext = "EarI";
    }
    if (($lig_ =~ "CTCTTC")||($lig_ =~ "GAAGAG")){
	$re_lig = "EarI";
    }
    return ($re_ext, $re_lig);
}

sub load_seq   # Load reference bases
{
    my $seq = "";
    my $chr = shift @_;
    my $path = shift @_;
    my $fname = "$path$chr.fa";
    open (REF, "<$fname") or die "Can't open $fname.\n";
    while (my $line = <REF>){
	chop($line);
	if ($line =~ ">"){
	}else{
	    $seq .= $line;
	}
    }
    close(REF);
    return $seq;
}
