#!/usr/bin/perl
use strict; 

my $synt_prt = "./auto_qsub_trim_novoN.pl list\n";
my $in_fname = $ARGV[0];
#sample barcode index
#even50  AGTGTC  ATCACG
#odd50   GTGCTA  ATCACG
#both50  TCGATG  ATCACG

if (($#ARGV != 0)) {
    print $synt_prt;
    exit;
}

my $name;
open(IN, "<$in_fname") or die "Can't open $in_fname.\n";
while(<IN>){
    if (substr($_, 0, 1) eq "#"){next;}
    ($name, ) = split(/\s+/, $_);
    my $out_fname = "trimN_$name\.sh";

    open(OUT, ">$out_fname");
    print OUT "#!/bin/bash\n";
    print OUT "./armseq2faN.pl /BiO/duMIP/raw0203/$name\_R1.fastq /BiO/duMIP/raw0203/$name\_R2.fastq 0 /BiO/duMIP/temp0203/trimN/$name\.armseq.fa\n";
    print OUT "novoalign -d /BiO/duMIP/CP230/CP230.ndx -f /BiO/duMIP/temp0203/trimN/$name\.armseq.fa -o SAM > /BiO/duMIP/temp0203/trimN/$name\.arminfo.sam\n";
    print OUT "./trim_MIP_product_w_novoN.pl /BiO/duMIP/CP230/CP230.duMIP /BiO/duMIP/temp0203/trimN/CP230.$name /BiO/duMIP/raw0203/$name\_R1.fastq /BiO/duMIP/raw0203/$name\_R2.fastq /BiO/duMIP/temp0203/trimN/$name\.arminfo.sam /BiO/duMIP/temp0203/trimN/$name\_R1.trim.fastq /BiO/duMIP/temp0203/trimN/$name\_R2.trim.fastq 0\n";
    close(OUT);
    
    my $sys_str = "qsub -cwd -j y -l hostname=compute-0-0.local $out_fname";
    print $sys_str, "\n";
    system($sys_str);
}
close(IN);
exit;
