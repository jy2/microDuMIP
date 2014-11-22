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

my %Gene;
my @chr;
my @b;
my @e;
my @jud;
my $exon_len=0;
my $temp;

my @chunks;
my $cur_gene;
my $old_gene = "Start";
my $num_gene=0;
my $num_exon=-1;
my $num_probe=0;
my $pass_len;

open(IN, "<$in_fname") or die "Can't open $in_fname.\n";
while (my $line=<IN>){
    @chunks = split(/\s+/, $line);
    if ($chunks[0] =~ "#"){
	next;
    }
    ($cur_gene, $temp, ) = split(/\./, $chunks[1]);
    if ($temp > $num_exon){
	$num_exon = $temp;
    }

    if ($cur_gene ne $old_gene){
	if ($old_gene ne "Start"){
       # calc. and print the result
	    $pass_len = calc_region();
	    print join("\t", $cur_gene, $num_exon+1, $num_probe, $exon_len, $pass_len, int($pass_len/$exon_len*100+0.5)), "\n";
	}else{
	    print join("\t", "Gene", "#EXON", "#PROBE", "TARGET(bp)", "COVERED(bp)", "percent"), "\n";
	}
	$old_gene = $cur_gene;
	$num_gene++;
	$num_exon=-1;
	$num_probe=0;
	$exon_len=0;
    }

    $temp = $chunks[2] . ":" . $chunks[3];
    if ($Gene{$cur_gene}->{$temp} eq ""){
	$exon_len+=($chunks[4] - $chunks[3]);
	$Gene{$cur_gene}->{$temp} = "!";
    }

    if ($chunks[5] ne "fail"){
	$chr[$num_probe] = $chunks[2];
	$b[$num_probe] = $chunks[6];
	$e[$num_probe] = $chunks[7];
	$jud[$num_probe] = 0;
	$num_probe++;
    }else{

    }
}
close(IN);
#print "MIP Pass: $pass\tFail: $fail\t", $pass/($pass+$fail)*100, "\n";
exit;

sub max($$){ $_[$_[0] < $_[1]] }
sub min($$){ $_[$_[0] > $_[1]] }
sub calc_region()
{
    for (my $i_=0; $i_<$num_probe; $i_++){
	if ($jud[$i_]==0){
	    for (my $j_=0; $j_<$num_probe; $j_++){
		if (($i_!=$j_)&&($jud[$j_]==0)&&($chr[$i_] eq $chr[$j_])){ 
		    if ( ($e[$i_]<$b[$j_])||($e[$j_]<$b[$i_]) ){
			# separate regions
		    }else{
			# merge!
			$b[$i_] = min($b[$i_], $b[$j_]);
			$e[$i_] = max($e[$i_], $e[$j_]);
			$jud[$j_]=1;
		    }
		}		
	    }
	}
    }
    
    my $return_=0;
    for(my $i_=0; $i_<$num_probe; $i_++){
	if ($jud[$i_]==0){
	    $return_ += ($e[$i_]-$b[$i_]);
	}
    }
    return $return_;
}
