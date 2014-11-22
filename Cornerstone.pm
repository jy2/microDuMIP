package Cornerstone;

# parameters for sub shortOligoTm($)
my $R = 1.987;
my %deltaH;
my %deltaS;

sub del_space{
	my $string_ = @_[0];
	chomp($string_);
	for (my $k=0;$k<length($string_);$k++) {
		if (!(substr($string_,$k,1) eq " ")) {
			my $return_str = substr($string_, $k, length($string_)-$k+1);
			return $return_str;
		}
	}
}


# generate the complementary sequences (capital letters only)
sub reverse_seq($){
	my $cur_seq = @_[0];
	my $cur_len_seq = length($cur_seq);

	my $seq_rev = "";
	for(my $temp_i=0;$temp_i<$cur_len_seq;$temp_i++){
		if (substr($cur_seq, $temp_i, 1) eq "A") {
			substr($seq_rev, $temp_i, 1) = "T";
		} elsif (substr($cur_seq, $temp_i, 1) eq "T") {
			substr($seq_rev, $temp_i, 1) = "A";
		} elsif (substr($cur_seq, $temp_i, 1) eq "C") {
			substr($seq_rev, $temp_i, 1) = "G";
		} elsif (substr($cur_seq, $temp_i, 1) eq "G") {
			substr($seq_rev, $temp_i, 1) = "C";
		} else{
			print "Err: $cur_seq[$temp_i]|\n";
#			exit;
		}
	}
	return $seq_rev;
}

# generate random sequences with AT percentage
sub random_seq_generator_AT($$){
	my $len_random_seq = @_[0];

	#calculate the percentages of each nucleotide
	my $percentage_AT = sprintf("%.f", (@_[1]*10));
	my $percentage_A = sprintf("%.f", ($percentage_AT/2));
	my $percentage_GC = 1000 - $percentage_AT;
	my $percentage_G = sprintf("%.f", ($percentage_GC/2));
	my $percentage_ATG = $percentage_AT+$percentage_G;

	#prepare the random pool
	my @random_pool=();
	for(my $temp_i=0; $temp_i<$percentage_A; $temp_i++){
		$random_pool[$temp_i] = "A";
	}
	for (my $temp_i=$percentage_A;$temp_i<$percentage_AT;$temp_i++) {
		$random_pool[$temp_i] = "T";
	}
	for (my $temp_i=$percentage_AT;$temp_i<$percentage_ATG;$temp_i++) {
		$random_pool[$temp_i] = "G";
	}
	for (my $temp_i=$percentage_ATG;$temp_i<1000 ;$temp_i++) {
		$random_pool[$temp_i] = "C";
	}
	$random_pool[1000] = "unused";
#	print "A: $percentage_A, AT: $percentage_AT, G: $percentage_G, ATG: $percentage_ATG\n";

	#generate random seq.
	my $return_ = "";
	for (my $temp_i=0;$temp_i<$len_random_seq;$temp_i++) {
		$return_ .= $random_pool[int rand($#random_pool)]; 
	}
	return $return_;
}

#load the parameter file for calculating Tm
sub readNNParam($){
	my $NNParamFile = @_[0];
	open(NNFILE, $NNParamFile) || die("Error in opening \"$NNParamFile\" file!");
	my $line = <NNFILE>;
	while($line = <NNFILE>){
		chop($line);
		my ($seqF, $seqR, $dH, $dS, $mismatch) = split(/[:\t]/, $line);
		if(!$mismatch){
			$deltaH{'pm'}->{$seqF} = $dH;
			$deltaS{'pm'}->{$seqF} = $dS;
		}else{
			$deltaH{'mm'}->{$seqF}->{$seqR} = $dH;
			$deltaS{'mm'}->{$seqF}->{$seqR} = $dS;			
		}		
	}
	close(NNFILE);
}

sub shortOligoGC($){
    my $seq = shift @_;
    my @num_base; #AT, GC, others
    my $cur_base;
    my $GC_percent;
    for (my $i_=0; $i_<=length($seq); $i_++){
	$cur_base = uc(substr($seq, $i_, 1));
	if (($cur_base eq "A")||($cur_base eq "T")){
	    $num_base[0]++;
	} elsif(($cur_base eq "C")||($cur_base eq "G")){
	    $num_base[1]++;
	} else{
	    $num_base[2]++;
	}
    }
    # neglect other bases
    if ($num_base[0]+$num_base[1]==0){
	print "Err: No ATCG in $seq.";
	exit;
    }else{
	$GC_percent = int(100*$num_base[1]/($num_base[0]+$num_base[1]));
    }
    return $GC_percent;
}


#Calculate the Tm of given sequences.
sub shortOligoTm($){
	my $seq = @_[0];

	my $C_primer = 200;
	my $C_Mg = 1.5;
	my $C_MonovalentIon = 50;
	my $C_dNTP = 0.4;
	my $percentage_DMSO = 0;
	my $percentage_annealed = 50;       
	    $seq =~ s/[ \t\n]+//g;
        $percentage_annealed = 50.0 if (!$percentage_annealed);
        $percentage_annealed /= 100.0;

        my $C_SodiumEquivalent = $C_MonovalentIon + 120 * sqrt($C_Mg-$C_dNTP);
        my $seqLength = length($seq);
        my $dH = $deltaH{'pm'}->{substr($seq, 0, 1)} + $deltaH{'pm'}->{substr($seq, $seqLength-1, 1)};
        my $dS = $deltaS{'pm'}->{substr($seq, 0, 1)} + $deltaS{'pm'}->{substr($seq, $seqLength-1, 1)};
        $seq = uc($seq);
        for(my $temp_i = 0; $temp_i < $seqLength - 1; $temp_i ++){
                $dH += $deltaH{'pm'}->{substr($seq, $temp_i, 2)};
                $dS += $deltaS{'pm'}->{substr($seq, $temp_i, 2)};
#                print "$seq $temp_i ", substr($seq, $temp_i, 2), "\n" if(!$deltaS{'pm'}->{substr($seq, $temp_i, 2)});
        }
        $dS += 0.368 * $seqLength * log($C_SodiumEquivalent/1000.0);
        my $Tm = sprintf("%5.2f", ($dH * 1000) / ($dS + $R * (log($C_primer*(1-$percentage_annealed)/$percentage_annealed)-21.4164)) - 273.15 - 0.75*$percentage_DMSO);
        return $Tm;
}

#Calculate the dG at 37'C of given sequences.
sub shortOligodG($){
	my $seq = @_[0];

	my $C_primer = 200;
	my $C_Mg = 1.5;
	my $C_MonovalentIon = 50;
	my $C_dNTP = 0.4;
	my $percentage_DMSO = 0;
	my $percentage_annealed = 50;       
	    $seq =~ s/[ \t\n]+//g;
        $percentage_annealed = 50.0 if (!$percentage_annealed);
        $percentage_annealed /= 100.0;

        my $C_SodiumEquivalent = $C_MonovalentIon + 120 * sqrt($C_Mg-$C_dNTP);
        my $seqLength = length($seq);
        my $dH = $deltaH{'pm'}->{substr($seq, 0, 1)} + $deltaH{'pm'}->{substr($seq, $seqLength-1, 1)};
        my $dS = $deltaS{'pm'}->{substr($seq, 0, 1)} + $deltaS{'pm'}->{substr($seq, $seqLength-1, 1)};
        $seq = uc($seq);
        for(my $temp_i = 0; $temp_i < $seqLength - 1; $temp_i ++){
                $dH += $deltaH{'pm'}->{substr($seq, $temp_i, 2)};
                $dS += $deltaS{'pm'}->{substr($seq, $temp_i, 2)};
#                print "$seq $temp_i ", substr($seq, $temp_i, 2), "\n" if(!$deltaS{'pm'}->{substr($seq, $temp_i, 2)});
        }
        $dS += 0.368 * $seqLength * log($C_SodiumEquivalent/1000.0);
        my $dG = sprintf("%5.2f", $dH - (273.15 + 37.0)*$dS/1000.0);
        return $dG;
}


1;
