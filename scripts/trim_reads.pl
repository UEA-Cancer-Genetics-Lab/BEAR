#!/usr/bin/env perl 

use warnings;
use strict;
use Bio::SeqIO;
use Getopt::Long;
use List::Util qw(max);

my $input="";
my $fasta="";
my $output="";
my $err_rate="";
my $err_qual="";
my $matrix_file="";
my $h_ins=0;
my $h_del=0;

my @lengths;
my @nucleotides = ('A', 'T', 'G', 'C');
my $trainingLength = 200; # there is 1200 used which is too long
my @seeds; #first quality value of each sequence

GetOptions ("i=s" => \$input,
	    "f=s" => \$fasta,
	    "o=s" => \$output,
	    "r=s" => \$err_rate,
	    "q=s" => \$err_qual,
	    "m=s" => \$matrix_file,
	    "hi=f" => \$h_ins,
	    "hd=f" => \$h_del);

# error rates processing
my %rates;
if($err_rate ne "0"){
	system("error_models.py $err_rate > $err_rate.r.out");
	open(ERR_RATE_FILE, "$err_rate.r.out") or die "Could not open $err_rate.r.out\n"; 
	while(<ERR_RATE_FILE>){
		my ($m, $b);
		my $next;
		chomp;
		if($_ =~ /(.*)([ATGCX])(\.)(nonzero)(.*)/){
			my $base = $2;
			$next = <ERR_RATE_FILE>;
			chomp $next;
			if ($next =~ /(\s*)(-?\d*\.*\d*)(\s*)(-?\d*\.*\d*)/){
				($b, $m) = ($2, $4);
				$rates{$base}{'slope'} = $m;
				$rates{$base}{'intercept'} = $b;
				print $base."\t".$m."\t".$b."\n";
			}
		}
	}
	close(ERR_RATE_FILE);
}

# base quality error rate processing
my %quals;
my %indel_rates;
my %hash;
my $min_q = 1;
my $max_q = 40;

# initializing indel rates and uninformed priors in one loop
# assume insertion/deletion rates are equal unless otherwise specified
for my $pos (0..$trainingLength){
	$indel_rates{'D'}{$pos} = 0.5;
	$indel_rates{'I'}{$pos} = 1.0; # cumulative probability
	for (my $row = $min_q; $row <= $max_q; $row++){
		for (my $col = $min_q; $col <= $max_q; $col++){
			$hash{$pos}{$row}{$col} = 1; # uninformed prior
		}
	}
}

if($err_qual ne "0"){
	open(DRISEE_QUAL_FILE, $err_qual);
	my $i = 0;
	while(<DRISEE_QUAL_FILE>){
		($indel_rates{'D'}{$i}, $indel_rates{'I'}{$i}) = (split("\t", $_))[6, 7];
		$indel_rates{'I'}{$i} += $indel_rates{'D'}{$i}; # cumulative probability
		$i++;
	}
	close(DRISEE_QUAL_FILE);

	system("error_quality.py $err_qual > $err_qual.r.out");
	open(ERR_QUAL_FILE, "$err_qual.r.out") or die "Could not open $err_qual.r.out\n";
	while(<ERR_QUAL_FILE>){
	        my ($inter, $coef1, $coef2);
	        my $next;
	        chomp;
	        if( $_ =~ /(.*)([ATGCX])(\.)(nonzero)(.*)/ ){
	        	my $base = $2;
	        	$next = <ERR_QUAL_FILE>;
	        	chomp $next;
	        	if($next =~ /(\s*)(-?\d*\.*\d*)(\s*)(-?\d*\.*\d*)(\s*)(-?\d*\.*\d*)/){
					($inter, $coef1, $coef2) = ($2, $4, $6);
					$quals{$base}{'intercept'} = $inter;
					$quals{$base}{'coef1'} = $coef1;
					$quals{$base}{'coef2'} = $coef2;	        		
	        	}
	        }
	}
	close(ERR_QUAL_FILE);
}

my %ins_matrix;
my %sub_matrix;
if($matrix_file ne "0"){
	my $mat_count = 0; #should only read in 8 lines of input (4 lines insertion matrix, 4 lines substitution matrix
	open(MATRIX_FILE, $matrix_file) or die "Could not open $matrix_file\n";
	while(<MATRIX_FILE>){
		chomp;
		if(($_ =~ /(\d).*/) && ($mat_count < 8)){
			my @line = split("\t", $_);
			my $index = $mat_count < 4 ? $mat_count : $mat_count - 4;
			my $hash_ref = $mat_count < 4 ? \%ins_matrix : \%sub_matrix;
			my $cumProb = 0;
			my $base1 = $nucleotides[$index];
			for my $index2 (0..3){
				$cumProb += $line[$index2];
				${$hash_ref}{$base1}{$nucleotides[$index2]} = $cumProb;
			}
			$mat_count++;
		}
	}
	close(MATRIX_FILE);
}


#populate the markov chain based on previous value + current state
my $counter=0;
my $prev_q1 = my $prev_q2 = $min_q;
my $seq_in= Bio::SeqIO->new( -format => 'fastq', -file => $input);
while( my $cur_seq = $seq_in->next_seq() ){
	my $prev_states = 0;
	my $prev_avg = 0;
	push( @lengths, length($cur_seq->seq()) ); 
	$counter += 1;
	if ($counter % 1000000 == 0){
		print "Training: $counter\n";
	}
	my $position = 0;
	foreach my $cur_qual (split / /, $cur_seq->qual_text()){
		if($prev_states == 0){ 
			push(@seeds, $cur_qual);
			$prev_states += 1;
		}else{
			$prev_q2 = $prev_q1;
			$prev_states += 1;	
		}
		$prev_q1 = $cur_qual;
		$hash{$position}{$prev_q2}{$cur_qual} += 10;
		$position += 1;
	}
}

my %rowsum;
my $length = max(@lengths);

for (my $pos = 0; $pos <= $length; $pos++){
	for (my $row = $min_q; $row <= $max_q; $row++){
		$rowsum{$pos}{$row} = 0;
		for ( my $col = $min_q; $col <= $max_q; $col++){
			$rowsum{$pos}{$row} += $hash{$pos}{$row}{$col};
		}
		if ($rowsum{$pos}{$row} == 0){
			$rowsum{$pos}{$row} = 1;
		}
		my $cumulateProb = 0;
		for (my $col = $min_q; $col <= $max_q; $col++){
			$cumulateProb += $hash{$pos}{$row}{$col}/$rowsum{$pos}{$row};
			$hash{$pos}{$row}{$col} = $cumulateProb;
		}
	}
}


#generate sequences
$counter=0;
open(MYFILE, ">$output");
my $header = ' ';
my $readSeq = ' ';
open(MYFILE3, $fasta);
while(<MYFILE3>){
	my $line = $_;
	chomp($line);
	if($line =~ s/>/@/){
		$header = $line;
		next;
	}else{
		$readSeq = $line;
	}
	$counter += 1;
	if ($counter % 1000000 == 0){
		print "Generating $counter\n";
	}

	my $state_num = 0;
	my $prev1 = my $prev2 = my $prev3 = my $prev4 = my $prev5 = my $prev6 = my $avg = $min_q;
	my $markov_qual = 0;
	$markov_qual = $seeds[rand @seeds];
	my @qual_string=();
	# first base quality set here, so below 'for' loop starts from 2ndbase
	push(@qual_string, chr($markov_qual+33));

	my $read_length = length($readSeq);
	for (my $i = 1; $i < $read_length; $i++){
		my $avg2 = $avg;
		my $rand_per = rand();
		my $cur_nuc = substr($readSeq, $i, 1);
		my $sub_check = rand();
		my $indel_check = rand();
		if ($state_num == 0){
			$avg = $min_q;
			$state_num += 1;
		}elsif ($state_num == 1){
			$prev2 = $prev1;
			$avg = $prev1;
			$state_num += 1;
		}elsif ($state_num == 2){
			$prev3 = $prev2;
			$prev2 = $prev1;
			$avg = int( ($prev2 + $prev3 ) / 2.0);
			$state_num += 1;
		}elsif ($state_num == 3){
			$prev4 = $prev3;
			$prev3 = $prev2;
			$prev2 = $prev1;
			$avg = int ( ($prev2 + $prev3 + $prev4) / 3.0);
			$state_num += 1;
		}elsif($state_num == 4){
			$prev5 = $prev4;
			$prev4 = $prev3;
			$prev3 = $prev2;
			$prev2 = $prev1;
			$avg = int ( ($prev2 + $prev3 + $prev4 + $prev5) / 4.0);
			$state_num += 1;
		}else{
			$prev6 = $prev5;
			$prev5 = $prev4;
			$prev4 = $prev3;
			$prev3 = $prev2;
			$prev2 = $prev1;
			$avg = int( ($prev2*.38) + ($prev3*.29) + ($prev4*.20) + ($prev5*.11) + ($prev6*.02));
		}

		#handle ambiguous nucleotides, smartmatch requires perl 5.10
		my $nuc_check;
		if( $cur_nuc ~~ @nucleotides){
			$nuc_check = $cur_nuc;
		}else{
			$nuc_check = $nucleotides[rand @nucleotides];
		}
		$prev1 = $markov_qual;

		# setting sub_rate and indel_rate to 0 for "N" nucleotide position
		my $sub_rate = 0;
		my $indel_rate = 0;
		unless ($cur_nuc eq "N") {
			my $sub_rate_check = $err_rate ne "0" ? exp( ($rates{$cur_nuc}{'slope'} * $i) + $rates{$cur_nuc}{'intercept'}) : 10**((-$markov_qual)/10);
			#my $sub_rate_check = $err_rate ne "0" ? $rates{$cur_nuc}{'slope'} * exp( exp($rates{$cur_nuc}{'intercept'}) * $i) : 10**((-$markov_qual)/10);
			$sub_rate = $sub_rate_check > 1.0 ? 1.0 : $sub_rate_check;

			#my $sub_rate = $err_rate ne "0" ? exp( ($rates{$cur_nuc}{'slope'} * $i) + $rates{$cur_nuc}{'intercept'})/100 : 10**((-$markov_qual)/10);
			#my $indel_rate = $err_rate ne "0" ? exp( ($rates{'X'}{'slope'} * $i) + $rates{'X'}{'intercept'})/100 : 0;
			#my $indel_rate_check = $err_rate ne "0" ? $rates{'X'}{'slope'} * exp( exp($rates{'X'}{'intercept'}) * $i) : 0;
			my $indel_rate_check = $err_rate ne "0" ? exp( ($rates{'X'}{'slope'} * $i) + $rates{'X'}{'intercept'}) : 0;
			$indel_rate = $indel_rate_check > 1.0 ? 1.0 : $indel_rate_check;
		}

		# sanity checking - base position and  generated qulaity length are equal
		my $seq_chunk = substr($readSeq, 0, $i);
		if($i != scalar(@qual_string)){
			my $els = scalar(@qual_string);
			print "POS $i\tELS $els\n";
			print "$readSeq\n",@qual_string,"\n";
			die;
		}

		$avg = $prev1;
		my $type=rand();
		my $new_base = $nuc_check;
		if($sub_check < $sub_rate){
			if($err_rate ne "0"){
				my $nuc_prob = rand();
				for my $base_check (@nucleotides){
					if($nuc_prob < $sub_matrix{$nuc_check}{$base_check}){
						$new_base = $base_check;
						last;
					}
				}
				# substitute new base in the sequence
				substr($readSeq, $i, 1) = $new_base;
				#my $sub_qual = int( ($quals{$nuc_check}{'slope'} * $i) + $quals{$nuc_check}{'intercept'});
				my $new_sub_qual = int($quals{$nuc_check}{'intercept'} + ($quals{$nuc_check}{'coef1'} * $i) + (($quals{$nuc_check}{'coef2'}) * ($i**2)));
				my $sub_qual;
				if ($new_sub_qual < $min_q){
					$sub_qual = $min_q;
				}elsif($new_sub_qual > $max_q){
					$sub_qual = $max_q;
				}else{
					$sub_qual = $new_sub_qual;
				}
				push(@qual_string, chr($sub_qual+33));
				$markov_qual = $sub_qual;
			}else{
				substr($readSeq, $i, 1) = $nucleotides[rand @nucleotides];
				my $sub_qual = $qual_string[$#qual_string];
				push(@qual_string, $sub_qual);
				$markov_qual = ord($sub_qual)-33;
			}
		}else{
			for(my $col = $max_q; $col >= ($min_q); $col--){
				$avg2 = $avg;	
				if( ($rand_per == $hash{$i}{$avg}{$col}) ){
					push(@qual_string, chr($col+33));
					$markov_qual = $col;
					last;
				}elsif ( ($col != $min_q) && ($rand_per < $hash{$i}{$avg}{$col}) && ($hash{$i}{$avg}{$col-1} <= $rand_per)){
					push(@qual_string, chr($col+33));
					$markov_qual = $col;
					last;
				}elsif( ($col == $min_q+1) && ($hash{$i}{$avg}{$col-1} > $rand_per)) {
					push(@qual_string, chr($col-1+33));
					$markov_qual = $col-1;
					last;
				}elsif( ($col == $min_q)){ #&& ($hash{$i}{$avg}{$col} > $rand_per)){
					push(@qual_string, chr($col+33));
					$markov_qual = $col;
					last;
				}else{
					#print "RANDPER: $rand_per\n POS: $i\nAVG: $avg \n QUAL: $col $blah \n COL $hash{$i}{$avg}{$col}\n COL-1 $hash{$i}{$avg}{$col-1}\n";
				}
			}
			my $temp_els = scalar(@qual_string);
			my $prev_base = substr($readSeq, $i-1, 1);
			my $h_total = $h_ins + $h_del > 0 ? $h_ins + $h_del : 1;
			my $h_i_per = $h_ins / ($h_total);
			my $h_d_per = $h_del / ($h_total);
			#if($nuc_check eq $prev_base){
				my $h_check = rand();
				my $h_max = 12;
				my $h_cur = 0;
				while($h_check < $h_total && $h_cur < $h_max && $h_total != 1){
					my $h_err_type = rand();
					if($h_err_type < $h_i_per){
					#homopolymer insertion
						#print "OLD SEQ: $readSeq\n";
						substr($readSeq, $i+1, 0, $nuc_check);
						#print "NEW SEQ: $readSeq\n";
						my $new_ins_qual = int($quals{'X'}{'intercept'} + ($quals{'X'}{'coef1'} * $i) + (($quals{'X'}{'coef2'}) * ($i**2)));
						my $ins_qual;
						if($new_ins_qual > $max_q){
							$ins_qual = $max_q;
						}elsif($new_ins_qual < $min_q){
							$ins_qual = $min_q;
						}else{
							$ins_qual = $new_ins_qual;
						}
						push(@qual_string, chr($ins_qual+33));
						$i++;
						$read_length++;
						$markov_qual = $ins_qual;
					}else{
						substr($readSeq, $i, 1) = "";
						$i--;
						$read_length--;
						$markov_qual = ord(pop(@qual_string))-33;
					}
					$h_check = rand();
					$h_cur++;
				}
			#}


			if($indel_check < $indel_rate){
				if($type>$indel_rates{'D'}{$i}){ #insertion
					my $nuc_prob = rand();
					for my $base_check (@nucleotides){
						if($nuc_prob < $ins_matrix{$nuc_check}{$base_check}){
							$new_base = $base_check;
							last;
						}
					}
					substr($readSeq, $i+1, 0, $new_base);
					my $new_ins_qual = int($quals{'X'}{'intercept'} + ($quals{'X'}{'coef1'} * $i) + (($quals{'X'}{'coef2'}) * ($i**2)));
					my $ins_qual;
					if($new_ins_qual > $max_q){
						$ins_qual = $max_q;
					}elsif($new_ins_qual < $min_q){
						$ins_qual = $min_q;
					}else{
						$ins_qual = $new_ins_qual;
					}
					#my $ins_qual = int( ($quals{'X'}{'slope'} * $i) + $quals{'X'}{'intercept'});
					push(@qual_string, chr($ins_qual+33));
					$i++;
					$read_length++;
					$markov_qual = $ins_qual;
				}else{ #deletion
					substr($readSeq, $i, 1) = "";
					$i--;
					$read_length--;
					$markov_qual = ord(pop(@qual_string))-33;
				}
			}
		}
	}
	$prev1 = 0;
	$prev2 = 0;
	$prev3 = 0;
	$prev4 = 0;
	$prev5 = 0;
	$avg = 0;

	my $trimqual = join "", @qual_string;
	print MYFILE "$header\n$readSeq\n+\n$trimqual\n";
}
close(MYFILE);
