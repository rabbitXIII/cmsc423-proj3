#!/usr/bin/perl -w

#
# FASTA Parser
# Rohit Gopal
#
# A local alignment tool for project part 3
# based off of part 2

use warnings;
use strict;
use Bio::SeqIO;

my $file1 = $ARGV[0];
my $file2 = $ARGV[1];
my $hoxd = $ARGV[2];
my $match = $ARGV[3];
my $mismatch = $ARGV[4];
my $gap_open = $ARGV[5];
my $gap_extend = $ARGV[6];


#read files into seq objects

my $bio_file1 = Bio::SeqIO->new(-file => $file1);
my $bio_file2 = Bio::SeqIO->new(-file => $file2);



#perform the global alignment

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
	#	-1 	left		gap to the right
	#	0	middle		match / mismatch
	#	1	up		gap down
my @matrix = undef;
my @traceback = undef;

	# O(n)
for my $index (0..$length1) {
	$matrix[$index][0] = $gap_open * $index;
	$traceback[$index][0] =  -1;
}

	# O(n)
for my $index (1..$length2) {
	$matrix[0][$index] = $gap_open * $index;
	$traceback[0][$index] = 1; 
}

	#O(n^2) no way around this..? ;(
for my $i (1..$length1) {
	for my $j (1..$length2) {
		my $score = undef;	
		if ($seq1[$i-1] eq $seq2[$j-1]) {
			$score = $match;
		} else {
			$score = $mismatch;
		} 
		my $diagonal = $matrix[$i-1][$j-1] + $score;
		my $down_gap = $matrix[$i][$j-1] + $gap_open;
		my $right_gap = $matrix[$i-1][$j] + $gap_open;
		if ($right_gap >= $diagonal && $right_gap >= $down_gap) {
			$matrix[$i][$j]=$right_gap; 
			$traceback[$i][$j] = -1;
		} elsif ($diagonal > $down_gap && $diagonal > $right_gap) { 
			$matrix[$i][$j]=$diagonal; 
			$traceback[$i][$j] = 0;
		} else { 
			$matrix[$i][$j]=$down_gap; 
			$traceback[$i][$j] = 1;
		}
	}
}

#print_2d(@matrix);
#print_2d(@traceback);

#do the traceback

my $i = $length1;
my $j = $length2;
my @result = undef;
my $index = 0;
	#we need to print backwards into the array since we don't know how many gaps there are
while ($i != 0 || $j != 0) {
	$result[1][$index] = " ";
	if ($traceback[$i][$j] == 0) { #diagonal match / mismatch
		$result[1][$index] = "|" if ($seq1[$i-1] eq $seq2[$j-1]);
		$result[0][$index] = $seq1[--$i];
		$result[2][$index] = $seq2[--$j];
	} elsif ($traceback[$i][$j] == -1) {  
		$result[0][$index] = $seq1[--$i];
		$result[2][$index] = "-";
	} else {
		$result[0][$index] = "-";
		$result[2][$index] = $seq2[--$j];
	}
	$index++;
}
my $score = $matrix[-1][-1];

#print out the score and alignment
print_result($index, $score,  @result);

#######################################
#           extra functions           #
#######################################

sub print_2d {
	my @array_2d=@_;
	for(my $i = 0; $i <= $length1; $i++){
	   for(my $j = 0; $j <= $length2 ; $j++){
	      print "$array_2d[$j][$i]\t";
	   }
	   print "\n";
	}
}

sub print_result {
	my $length = shift;
	my $score = shift;
	my @array = @_;
	my $header1 = $seq_obj1->id . " ";
	my $header2 = " " x length($header1);
	my $header3 = $seq_obj2->id . " ";

	# print the edit distance

	print "Edit Distance = $score\n";

	my ($result_line1, $result_line2, $result_line3) = (undef, undef, undef);
	
	
	for $i (1..$length) {
		$result_line1 .= $array[0][0-$i];
		$result_line2 .= $array[1][0-$i];
		$result_line3 .= $array[2][0-$i];
	}
	
	my $remaining_length = 80 - length($header1);

	while (length($header1 . $result_line1) > 80) {
		print $header1;
		print substr( $result_line1, 0, $remaining_length) . "\n";
		$result_line1 = substr $result_line1, $remaining_length;
		print $header2;
		print substr( $result_line2, 0, $remaining_length) . "\n";
		$result_line2 = substr $result_line2, $remaining_length;
		print $header3;
		print substr( $result_line3, 0, $remaining_length) . "\n";
		$result_line3 = substr $result_line3, $remaining_length;
	}
	print $header1 . $result_line1 . "\n";
	print $header2 . $result_line2 . "\n";
	print $header3 . $result_line3 . "\n";
	
}	
