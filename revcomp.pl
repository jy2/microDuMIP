#!/usr/bin/perl
use strict;

my $in = shift;
my $out = rev_comp($in);
print join("\t", $in, "->", $out),"\n";
exit;

sub rev_comp
{
    my $in_seq = shift;
    my $cur_let;
    my $out_seq = "";
    $in_seq = reverse($in_seq);
    for (my $i_=0; $i_<length($in_seq); $i_++){
	$cur_let = uc(substr($in_seq, $i_, 1));
	if ($cur_let eq "G"){
	    $out_seq .= "C";
	}elsif($cur_let eq "C"){
	    $out_seq .= "G";
	}elsif($cur_let eq "A"){
	    $out_seq .= "T";
	}elsif($cur_let eq "T"){
	    $out_seq .= "A";
	}elsif($cur_let eq "N"){
	    $out_seq .= "N";
	}else{
	    print "Err : Wrong base $cur_let in $in_seq\n";
#	    print_out();
	}
    }
    return $out_seq;
}
