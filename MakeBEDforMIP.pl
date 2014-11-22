#!/usr/bin/perl

# This program is to make input file (BED format) for duplex molecular inversion probe(duMIP) from refSeq tract in UCSC Genome Browser
# Developed by Jung-Ki Yoon (neododari@gmail.com)

# Note : all starting points in this program is 0-based. http://genome.ucsc.edu/FAQ/FAQtracks#tracks1
# The real starting point is 1 base next to the position. 
 
use strict; 
# output file format (BED format)
# example:
# chr7    55086965        55087063        EGFR.1     <-- exon+splicing(5bp)
# chr7    55209973        55210135        EGFR.2

my $synt_prt = "./build_BED4MIP.pl [refSeq table] [GeneName] \n";
my $Gene_name = $ARGV[1];
my $refSeq = $ARGV[0];
my $out_fname = $Gene_name . ".BED";
if ($#ARGV != 1) {
    print $synt_prt;
    exit;
}

my $splicing = 5; # definition of splicing region +- 5bp 

my @chunks;
my @pos_b;
my @pos_e;
my $cur_chr;
my $num=0;
my $ID=0;
my @b;
my @e;
my @judge;

# 125     NM_201284       chr7    +       55086724        55238738        55086970        55238237        16      55086724,55209978,55210997,55214298,55218986,55220238,55221703,55223522,55224225,55224451,55225355,55227831,55229191,55231425,55232972,55237999,       55087058,55210130,55211181,55214433,55219055,55220357,55221845,55223639,55224352,55224525,55225446,55228031,55229324,55231516,55233130,55238738,       0       EGFR    cmpl    cmpl    0,1,0,1,1,1,0,1,1,2,1,2,1,2,0,2,
# 984     NM_020327       chr12   +       52347163        52390863        52369113        52387894        9       52347163,52369048,52370110,52374752,52377782,52378975,52380601,52385646,52387768,      52347283,52369288,52370359,52374983,52377950,52379132,52380726,52385777,52390863,       0       ACVR1B  cmpl    cmpl    -1,0,1,1,1,1,2,1,0,

my $coding_start;
my $coding_end;
open(OUT, ">$out_fname");
open(IN, "<$refSeq") or die "Can't open $refSeq.\n";
while (my $line=<IN>){
    @chunks = split(/\s+/, $line);
    if ((substr($chunks[1], 0, 2) eq "NM")&&($Gene_name eq $chunks[12])&&(!($chunks[2]=~ "_"))){
	$cur_chr = $chunks[2];
	@pos_b = split(/\,/, $chunks[9]);
	@pos_e = split(/\,/, $chunks[10]);
	# identifying coding exons
	for (my $i=0; $i<=$#pos_e; $i++){
#	    print join("\t", $pos_e[$i], $chunks[6]), "\n";
	    if ($pos_e[$i]>$chunks[6]){
		$coding_start=$i;
#		print $coding_start, "\n";
		last;
	    }
	}
	for (my $i=$#pos_b; $i>=0; $i--){
#	    print join("\t", $pos_b[$i], $chunks[7]), "\n";
	    if ($pos_b[$i]<$chunks[7]){
		$coding_end=$i;
#		print $coding_end, "\n";
		last;
	    }
	}

	for(my $i=$coding_start; $i<=$coding_end; $i++){
	    $b[$num] = $pos_b[$i] - $splicing;
	    $e[$num] = $pos_e[$i] + $splicing;
	    if ($i==$coding_start){   # coding bases only
		$b[$num] = $chunks[6];
	    }elsif($i==$coding_end){
		$e[$num] = $chunks[7];
	    }
	    $judge[$num] = 0;
	    $num++;
	}
    }
}
close(IN);

for (my $i=0; $i<$num; $i++){
    if ($judge[$i]==0){
	for (my $j=0; $j<$num; $j++){
	    if (($i!=$j)&&($judge[$j]==0)){ 
		if ( ($e[$i]<$b[$j])||($e[$j]<$b[$i]) ){
		    # separate regions
		}else{
		    # merge!
		    $b[$i] = min($b[$i], $b[$j]);
		    $e[$i] = max($e[$i], $e[$j]);
		    $judge[$j]=1;
		}
	    }		
	}
    }
}
my %out;

for (my $i=0; $i<$num; $i++){
    if ($judge[$i]==0){
	$out{$b[$i]}=$e[$i];
    }
}

foreach my $key (sort(keys(%out))){
    print OUT join("\t", $cur_chr, $key, $out{$key}, ($Gene_name . "." . $ID)), "\n";
    $ID++;
}
close(OUT);
exit;

sub max($$){ $_[$_[0] < $_[1]] }
sub min($$){ $_[$_[0] > $_[1]] }
