#!/usr/bin/perl -w

#
# FASTA Parser
# Rohit Gopal
#
# A local alignment tool for project part 3
# based off of part 2
#
# Accepts a HOXD Matrix and uses affine gap penalties
# Reports all optimal alignments


use warnings;
use strict;
use Bio::SeqIO;

my $file1 = $ARGV[0];
my $file2 = $ARGV[1];
my $hoxd = $ARGV[2];
my $gap_open = $ARGV[3];
my $gap_extend = $ARGV[4];


#read files into seq objects

my $bio_file1 = Bio::SeqIO->new(-file => $file1);
my $bio_file2 = Bio::SeqIO->new(-file => $file2);

my $hoxd_handle;
open($hoxd_handle, $hoxd);
my @hoxd_lines = <$hoxd_handle>;
close($hoxd_handle);

# read in the HOXD matrix to a hash
my %hoxd_scores = (	A => {},
			T => {},
			C => {},
			G => {});

foreach(@hoxd_lines){ 
	if($_ =~ /([ACTG]),([ACTG])=(-?\d+)/){
		$hoxd_scores{$1}{$2} = $3; 	
		$hoxd_scores{$2}{$1} = $3;
	}
}

#perform the local alignment

my $seq_obj1 = $bio_file1->next_seq;
my $seq_obj2 = $bio_file2->next_seq;	

my $seq1 = $seq_obj1->seq;
my $seq2 = $seq_obj2->seq;

	# we add one to make room for the initial row / col
my $length1 = length($seq1);
my $length2 = length($seq2);

	#doing the global alignment is easier with the strings being arrays

my @seq1 = split(//, $seq1);
my @seq2 = split(//, $seq2);


	#the traceback will have the following values
	#	1 	left		gap to the right
	#	2	middle		match / mismatch
	#	4	up		gap down
	#	3	left/middle	gap to right or match/mismatch
	#	5	left/up		gap right or gap down
	#	6	middle/up	match/mismatch or gap down
	#	7	all		all have same score
	# use the 1 2 4 method to keep track of traceback

my $right_trace = 1;
my $diag_trace = 2;
my $down_trace = 4;

my @matrix = undef;
my @traceback = undef;

#initialize the matrices
# this could honestly just be done in a single matrix, but it's easier
# to read
for my $index (0..$length1) {
	for my $index2 (0..$length2) {
		$matrix[$index][$index2] = 0;
		$traceback[$index][$index2] = 0;
	}
}	
	# O(n)
for my $index (0..$length1) {
	$traceback[$index][0] =  $right_trace;
}

	# O(n)
for my $index (1..$length2) {
	$traceback[0][$index] = $down_trace; 
}




my @max = ({ i => 0,
		j => 0,
		score => 0 });

	#O(n^2) no way around this..? ;(
for my $i (1..$length1) {
	for my $j (1..$length2) {
	
		my $score = $hoxd_scores{$seq1[$i-1]}{$seq2[$j-1]};
		my $diagonal = $matrix[$i-1][$j-1] + $score;
		my $down_gap = $matrix[$i][$j-1] + $gap_open;
		my $right_gap = $matrix[$i-1][$j] + $gap_open;
	
		if ($diagonal < 0 && $down_gap < 0 && $right_gap < 0 ) {
			$matrix[$i][$j] = 0;
			$traceback[$i][$j] = 0;
			next; 
		}
		if ($right_gap >= $diagonal && $right_gap >= $down_gap) {
			$matrix[$i][$j] = $right_gap; 
			$traceback[$i][$j] = $right_trace;
			#$traceback[$i][$j] += $diag_trace if $right_gap == $diagonal;
			#$traceback[$i][$j] += $down_trace if $right_gap == $down_gap;
		} elsif ($diagonal > $down_gap && $diagonal > $right_gap) { 
			$matrix[$i][$j] = $diagonal; 
			$traceback[$i][$j] = $diag_trace;
			#$traceback[$i][$j] += $down_trace if $diagonal == $down_gap;
		} else { 
			$matrix[$i][$j]=$down_gap; 
			$traceback[$i][$j] = $down_trace;
		}
		if ( $max[0]{score} < $matrix[$i][$j] ) {
			@max = undef;
			@max = ( { i => $i,
				j => $j,
				score => $matrix[$i][$j]} );
		} elsif ($max[0]{score} == $matrix[$i][$j]){
			push(@max, { i => $i, j => $j, score => $matrix[$i][$j]});

		}
		
	}
}

#do the traceback for each valid sub section

foreach(@max) {
	my ($i, $j, $score) = ($_->{i}, $_->{j}, $_->{score});	
	my @pointers = undef; # tracks the indexes of the sequences
	my @result = undef;
	my $index = 0;
	$pointers[1] = $i;
	$pointers[3] = $j;	
		#we need to print backwards into the array since we don't know how many gaps there are
	while ( ($matrix[$i][$j] != 0) && ($i != 0 || $j != 0)) {
		$result[1][$index] = " ";
		if ($traceback[$i][$j] == 2) { #diagonal match / mismatch
			$result[1][$index] = "|" if ($seq1[$i-1] eq $seq2[$j-1]);
			$result[0][$index] = $seq1[--$i];
			$result[2][$index] = $seq2[--$j];
		} elsif ($traceback[$i][$j] == 1) {  
			$result[0][$index] = $seq1[--$i];
			$result[2][$index] = "-";
		} else {
			$result[0][$index] = "-";
			$result[2][$index] = $seq2[--$j];
		}
		$index++;
	}
	$pointers[0] = $i;
	$pointers[2] = $j;	

	#print out the score and alignment
	print_result($index, $score, $pointers[0]+1, $pointers[1]+1, $pointers[2]+1, $pointers[3]+1,  @result);
}

#######################################
#           extra functions           #
#######################################

sub print_2d {
	my @array_2d=@_;
	for my $j (0..$length1) {
		for my $i (0..$length2){
			print "$array_2d[$j][$i]\t";
	   }
	   print "\n";
	}
}

sub print_result {
	my $length = shift;
	my $score = shift;
	my $start1 = shift;
	my $end1 = shift;
	my $start2 = shift;
	my $end2 = shift;
	my @array = @_;
	my $header1 = $seq_obj1->id . " " . $start1 . " ";
	my $header2 = " " x length($header1);
	my $header3 = $seq_obj2->id . " " . $start2 . " ";
	my $footer1 = " " . $end1;
	my $footer3 = " " . $end2;


	my $identities_count = 0;

	# print the edit distance


	my ($result_line1, $result_line2, $result_line3) = (undef, undef, undef);
	
	
	for my $i (1..$length) {
		$result_line1 .= $array[0][0-$i];
		$result_line2 .= $array[1][0-$i];
		$identities_count++ if ($array[1][0-$i] eq "|");
		$result_line3 .= $array[2][0-$i];
	}
	my $longest = ($end1 - $start1) > ($end2 - $start2) ? ($end1 - $start1) : ($end2 - $start2);
	my $percent = $identities_count / ($longest) * 100;
	printf "Score = %d, Identities = %d\/%d (%2.0f%%)\n", $score, $identities_count, $longest, $percent;	
	my $remaining_length = length($header1) > length($header3) ? 80 - length($header1) - length($footer1) : 80 - length($header3) - length($footer3);
	
	# make sure that the headers are the same length so that everything lines up
	
	while( length($header1) < length($header3) ) {
		$header1 .= " ";
	}

	while( length($header3) < length($header1) ) { 
		$header3 .= " ";
	}

	while (length($header1 . $result_line1) > 80) {
		print $header1;
		print substr( $result_line1, 0, $remaining_length);
		print $footer1 . "\n";
		$result_line1 = substr $result_line1, $remaining_length;
		print $header2;
		print substr( $result_line2, 0, $remaining_length) . "\n";
		$result_line2 = substr $result_line2, $remaining_length;
		print $header3;
		print substr( $result_line3, 0, $remaining_length);
		print $footer3 . "\n";
		$result_line3 = substr $result_line3, $remaining_length;
	}

	print $header1 . $result_line1 . $footer1 . "\n";
	print $header2 . $result_line2 . "\n";
	print $header3 . $result_line3 . $footer3 . "\n";
	
}	
