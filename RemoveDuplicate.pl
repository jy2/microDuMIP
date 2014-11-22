#!/usr/bin/perl
use strict;

# This program is to remove PCR duplicates of MIP capture reads by NNNNNNNNNN (10mer)
# Developed by JK Yoon (2013. 3. 23.)
# Modified for duMIP (2013. 8. 7.)
# samtools view -h *.bam -o *.sam   --> sort *.sam --> run this script
# This program will consider pair-end reads when both pairs have same NNNNNNNNNN and same align location on genome.
# and remove the read pair with lower mapping QS. If the mapping QS are identical, the read pair with lower base QS will be discarded. (Same approach to Picard MarkDuplicates)
# Modified the criteria of identical reads: MIP_ID:randomN:chr:posF:posR -> MIP_ID:randomN:chr (2013.10.13)

my $len_random = 10;

my $synt_str = "./RemoveDuplicate.pl Gene_MIP.anno input.sam output.sam\n";
my $anno_fname = $ARGV[0];
my $sam_fname = $ARGV[1];
my $out_fname = $ARGV[2];
my $log_fname = $ARGV[2] . ".log";

if ($#ARGV != 2){
    print $synt_str;
    exit;
}

my %MIP;
my @chunks;
my $num_MIP=0;
#num    gene.exon.MIP   chr     begin   end     MIP_type        cover_begin     cover_end       gap_begin       gap_end seq_EXT Tm_EXT GC_EXT LEN_EXT seq_LIG Tm_LIG GC_LIG LEN_LIG
#0       ABL1.0.0        chr9    133589706       133589847       normal          133589706       133589809       133589698       133589809     TTCTGGAAAGGGGTACCTATT  58.49  42 4_21 GGGGTCCACACTGCAATGTTTT 64.55 50 0_22

open(IN, "<$anno_fname") or die "Can't open $anno_fname\n";
while(my $line = <IN>){
    if (substr($line, 0, 1) eq "#"){
	next;
    }
    @chunks = split(/\s+/, $line);
    $MIP{$chunks[1]}->{'num'}=0;
    $MIP{$chunks[1]}->{'unique'}=0;
    $num_MIP++;
}
close(IN);
print "$num_MIP MIP were loaded.\n";

open(OUT, ">$out_fname");
open(LOG, ">$log_fname");
open(IN, "<$sam_fname") or die "Can't open $sam_fname\n";

print LOG "Log file of RemoveDuplicate\nlength of random seq = $len_random\n";
print LOG "Annotated MIP: $anno_fname\nSam input: $sam_fname\nOutput: $out_fname\nLog: $log_fname\n";
print "Annotated MIP: $anno_fname\nSam input: $sam_fname\nOutput: $out_fname\nLog: $log_fname\n";

#ATGAACCAAA_CREBBP.0.20:52:H78V4ADXX:2:1101:11230:2129   99      16      3778715 70      101M    =       3778721 108     ATTGGCCACGTACTTGGCTGTGCGCTGTTTGATGAAAGCTGCCATTAGCTGCGGGTTTGATTTGAGAATGTTCAGCACCTGCTGTTGCTGCTGAGGGGAGC   EHFCDCEEDCAEDBDDDDDDACDDD?DBCDDCDCACCCCDDDDDDCCDDCCD@BB0<@BDCCDCDDBCDD:CCDDDDDDDDDDDCCDDDDCDDDDDDD9BB   PG:Z:novoalign  RG:Z:both500    PU:Z:HiSeq2500  LB:Z:duMIP      AS:i:0  UQ:i:0  NM:i:0  MD:Z:101        PQ:i:2  SM:i:70 AM:i:70
#ATGAACCAAA_CREBBP.0.20:52:H78V4ADXX:2:1101:11230:2129   147     16      3778721 70      102M    =       3778715 -108    CACGTACTTGGCTGTGCGCTGTTTGATGAAAGCTGCCATTAGCTGCGGGTTTGATTTGAGAATGTTCAGCACCTGCTGTTGCTGCTGAGGGGAGCTGGGCGA  A?DDCCB?8DDC><BBBBCBDCCCA@DCA>A?DDBCACCCCBDBDDBBA@C@CCCA:CCDCCACACCCC:9BBCDDDDDDBDCDDDDCDFFFFFHFHGIGHG  PG:Z:novoalign  RG:Z:both500    PU:Z:HiSeq2500  LB:Z:duMIP      AS:i:0  UQ:i:0  NM:i:0  MD:Z:102        PQ:i:2  SM:i:70 AM:i:70

my %lines;
my %UNIQ;
my $line_F;
my $line_R;
my @chunks_F;
my @chunks_R;
my $check_str;
my $ID;
my $cur_MQ;
my $cur_BQ;
my $cur_ran;
my $total_unique=0;
my $num_line = 0;
while(my $line = <IN>){
    if (substr($line, 0, 1) eq "@"){
	print OUT $line;
	next;
    }
    ($ID, ) = split(/\s+/, $line);
    if($lines{$ID}->{'idx'} eq ""){
	$lines{$ID}->{'idx'}=1;
	$lines{$ID}->{'fwd'} = $line;
    }elsif($lines{$ID}->{'idx'} eq "1"){
	$lines{$ID}->{'idx'}=2;
	$lines{$ID}->{'rev'} = $line;
    }else{
	print "Err $line, $lines{$ID}->{'idx'}\n";
	exit;
    }
}
close(IN);

foreach my $head (keys %lines){
    $num_line++;
    if ($num_line % 10000 == 0){
	print "$num_line read pairs were processed. $total_unique were unique.\n";
    }
    $line_F = $lines{$head}->{'fwd'};
    $line_R = $lines{$head}->{'rev'};
    if ($lines{$head}->{'idx'} ne "2"){
	print "Err $lines{$head}->{'idx'}\n";
	next;
#	exit;
    }

    @chunks_F = split(/\s+/, $line_F);
    @chunks_R = split(/\s+/, $line_R);
    if ($chunks_F[0] ne $chunks_R[0]){
	print "Err : headers of pair were not same\n$line_F$line_R";
	exit;
    }

    ($ID, ) = split(/:/, $chunks_F[0]);
    $cur_ran = substr($ID, 0, $len_random); 
    substr($ID, 0, $len_random+1) = "";
    $check_str = "$ID:$cur_ran:$chunks_F[2]";
#    $check_str = "$ID:$cur_ran:$chunks_F[2]:$chunks_F[3]:$chunks_R[3]";
    $MIP{$ID}->{'num'}++;
    $cur_MQ = $chunks_F[4] + $chunks_R[4];
    $cur_BQ = calc_BQ($chunks_F[10]) + calc_BQ($chunks_R[10]);
    if($UNIQ{$check_str}->{'num'} eq ""){ # Same ID, randomN, chr, pos_F, pos_R
	$MIP{$ID}->{'unique'}++; # count number of unique reads
	$total_unique++;
	$UNIQ{$check_str}->{'num'}=1;
	$UNIQ{$check_str}->{'line'} = $line_F . $line_R;
	$UNIQ{$check_str}->{'MQ'} = $cur_MQ;
	$UNIQ{$check_str}->{'BQ'} = $cur_BQ; 
    }else{
	# pick remaining reads
	$UNIQ{$check_str}->{'num'}++;
	if ($cur_MQ > $UNIQ{$check_str}->{'MQ'}){
	    # exchange
	    $UNIQ{$check_str}->{'line'} = $line_F . $line_R;
	    $UNIQ{$check_str}->{'MQ'} = $cur_MQ;
	    $UNIQ{$check_str}->{'BQ'} = $cur_BQ;
	}elsif($cur_MQ == $UNIQ{$check_str}->{'MQ'}){
	    if ($cur_BQ > $UNIQ{$check_str}->{'BQ'}){
		# exchange
		$UNIQ{$check_str}->{'line'} = $line_F . $line_R;
		$UNIQ{$check_str}->{'MQ'} = $cur_MQ;
		$UNIQ{$check_str}->{'BQ'} = $cur_BQ;
	    }
	}
    }
}
print "REPORT: $num_line read pairs were processed. $total_unique were unique.\n";
print LOG join("\t", "MIP_ID", "read_pairs", "unique_read_pairs"), "\n";
foreach $ID (sort(keys %MIP)){
    print LOG join("\t", $ID, $MIP{$ID}->{'num'}, $MIP{$ID}->{'unique'}), "\n";
}
foreach $check_str (sort(keys %UNIQ)){
    print OUT $UNIQ{$check_str}->{'line'};
}
close(IN);
close(OUT);
close(LOG);
exit;

sub calc_BQ
{
    #Illumina 1.8+  phred+33 (0,41)
    my $in_ = shift @_;
    my $cur_BQ_;
    my $min_BQ_ = 15;
    my $sum_BQ_ = 0;
    for (my $i=0; $i<length($in_); $i++){
	$cur_BQ_ = ord(substr($in_, $i, 1))-33;
	if (($cur_BQ_ < 0)||($cur_BQ_>41)){
	    print "Err: Base quality score is not in Illumina 1.8+ format \n($cur_BQ_ <-- $i:$in_)\n";
	    exit;
	}
	if ($min_BQ_ <= $cur_BQ_){
	    $sum_BQ_ = $sum_BQ_ + $cur_BQ_;
	}
    }
    return $sum_BQ_;
}
