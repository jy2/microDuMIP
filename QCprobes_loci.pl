#!/usr/bin/perl
use strict; 

# #num gene.exon.MIP  chr     begin   end     MIP_type        cover_begin     cover_end       gap_begin       gap_end seq_EXT Tm_EXT  GC_EXT  LEN_EXT seq_LIG Tm_LIG  GC_LIGLEN_LIG  seq_preMIP
# 1 ARID1A.0.0      chr1    27022894        27024036        normal  27022894        27022994        27022886        27022994        GTCTCTCCGCGGACGAGACAG   67.16   66    4_21     AGGCGGGGGGCGAGGCGGCGGCGGC       86.79   92      0_25    GGTAGCAAAGTGCAGATGTGCTCTTCAAGGCGGGGGGCGAGGCGGCGGCGGCNNNNNNNNNNGTTGGAGGCTCATCGTTCCCAGATGTTATCGAGGTCCGACGTCTCTCCGCGGACGAGACAGCGGGGAGAAGTTCAGAGCAGTGTGAGTTCA

my $synt_prt = "./QC_design.pl [duMIP_design]\n";
my $in_fname = $ARGV[0];
if ($#ARGV != 0) {
    print $synt_prt;
    exit;
}

my %Exon;
my @chr;
my @b;
my @e;
my @jud;
my $num=0;
my $exon_len=0;
my $temp;

my @chunks;
my $cur_EXT;
my $cur_LIG;
my $Tm_interval = 5;
my $GC_interval = 10;
my $fail=0;
my $pass=0;
my %Tm;
my %GC;
my $gene;

open(IN, "<$in_fname") or die "Can't open $in_fname.\n";
while (my $line=<IN>){
    @chunks = split(/\s+/, $line);
    if ($chunks[0] =~ "#"){
	next;
    }
    ($gene, ) = split(/\./, $chunks[1]);
    $temp = $chunks[2] . ":" . $chunks[3];
    if ($Exon{$gene}->{$temp} eq ""){
	$exon_len += ($chunks[4] - $chunks[3]);
	$Exon{$gene}->{'exon_len'} += ($chunks[4] - $chunks[3]);
	$Exon{$gene}->{$temp} = "!";
    }
    if ($chunks[5] ne "fail"){
	$pass++;
	$chr[$num] = $chunks[2];
	$b[$num] = $chunks[6];
	$e[$num] = $chunks[7];
	$jud[$num] = 0;
	$num++;

	$cur_EXT = int($chunks[11]/$Tm_interval + 0.5);  # ex: 60 w/ 5 => 12 = 57.5~62.5
	$cur_LIG = int($chunks[15]/$Tm_interval + 0.5);
	$Tm{$cur_EXT}->{$cur_LIG}++;
	$cur_EXT = int($chunks[12]/$GC_interval + 0.5); 
	$cur_LIG = int($chunks[16]/$GC_interval + 0.5);
	$GC{$cur_EXT}->{$cur_LIG}++;
    }else{
	$fail++;
    }
}
close(IN);
print "MIP Pass: $pass\tFail: $fail\t", $pass/($pass+$fail)*100, "\n";

for (my $i=0; $i<$num; $i++){
    if ($jud[$i]==0){
	for (my $j=0; $j<$num; $j++){
	    if (($i!=$j)&&($jud[$j]==0)&&($chr[$i] eq $chr[$j])){ 
		if ( ($e[$i]<$b[$j])||($e[$j]<$b[$i]) ){
		    # separate regions
		}else{
		    # merge!
		    $b[$i] = min($b[$i], $b[$j]);
		    $e[$i] = max($e[$i], $e[$j]);
		    $jud[$j]=1;
		}
	    }		
	}
    }
}

my $pass_len=0;
for(my $i=0; $i<$num; $i++){
    if ($jud[$i]==0){
	$pass_len += ($e[$i]-$b[$i]);
    }
}
print "Target bp: $exon_len, Covered bp: $pass_len\t", $pass_len/$exon_len*100, "\n\n";

print "Tm\t";
for(my $j=0; $j<=int(100/$Tm_interval); $j++){
    print $j*$Tm_interval, "\t";
}
print "\n";
for(my $i=0; $i<=int(100/$Tm_interval); $i++){
    print $i*$Tm_interval, "\t";
    for(my $j=0; $j<=int(100/$Tm_interval); $j++){
	if ($Tm{$i}->{$j} eq ""){
	    $Tm{$i}->{$j} = 0;
	}
	print $Tm{$i}->{$j}, "\t";
    }    
    print "\n";
}
print "\n\n";

print "GC\t";
for(my $j=0; $j<=int(100/$GC_interval); $j++){
    print $j*$GC_interval, "\t";
}    
print "\n";
for(my $i=0; $i<=int(100/$GC_interval); $i++){
    print $i*$GC_interval, "\t";
    for(my $j=0; $j<=int(100/$GC_interval); $j++){
	if ($GC{$i}->{$j} eq ""){
	    $GC{$i}->{$j} = 0;
	}
	print $GC{$i}->{$j}, "\t";
    }    
    print "\n";
}
exit;

sub max($$){ $_[$_[0] < $_[1]] }
sub min($$){ $_[$_[0] > $_[1]] }
