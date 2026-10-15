#!/usr/bin/perl

###

use List::Util qw[max shuffle];
use strict;
#use warnings;

my $inputSeqFile = $ARGV[0];
my ($out_prefix) = $inputSeqFile =~ m/tmp_(.*)_tmp\.txt$/;
$out_prefix = "clean_".$out_prefix."_all.txt";




open(COORD,$inputSeqFile) || die "no such file\n";

open(OUT1, ">log_$out_prefix");
open(OUT2, ">$out_prefix");
#open(OUT3, ">$out_prefix"."_allInfo.txt");


#my $test_Seq = "GGTACCTGACCTCGCTAGCTATTTGTCTTGCTTGGAGGATCAGGGGTCAGGCCCTGTCCGTGACCTAGTTCAGGTTAAATAAAGCTGAGCTGGGAGGCTGTTTACTGCCTGTAGACCGCTGAGCAGAGTCCATGGCCCTGCCGGGCTGCCCTGGCTGGGGGCAGGAGCACGCCTGGGCATGGGGAAGCTTCTCGAGGCATGCACAGTTATGCGAGCTCTCGCTTCGAGCA";
#my $test_Seq = "GAGCTCGGTACCTGACCTCGCTAGCTATTTGTCTTGCTTGGAGGATCAGGGGTCAGGCCCTGTCCGTGACCTAGTTCAGGTTAAATAAAGCTGAGCTGGGAGGCTGTTTACTGCCTGTAGACCGCTGAGCAGAGTCCATGGCCCTGCCGGGCTGCCCTGGCTGGGGGCAGGAGCACGCCTGGGCATGGGGAAGCTTCTCGAGGCATGCACAGTTATGCGAGCTCTCGCTTCGAGCA";
#count_restriction($test_Seq);

my %dedup_hash;

### others
while (my $curr_rna = <COORD>) {
  chomp $curr_rna;
  my @tmp = split /\t/, $curr_rna;
  my $seq_here = uc($tmp[0]);
  if ((scalar @tmp) != 2 || length($seq_here) != 165) {
  	# body...
  	die "something went wrong: @tmp\n";
  }
	if ($seq_here =~/GCTAGC|AAGCTT|CTCGAG|GCATGC|GAGCTC|TCGCTTCGAGCA|GGTACCTGACCTC|^ATGC|^TAGC|^TCGAG/ || rc($seq_here) =~ /TCGCTTCGAGCA|GGTACCTGACCTC|^ATGC|^TAGC|^TCGAG/) { 
		print OUT1 "$tmp[1] has restriction site; give up\n";
		next;
	}
	if(exists($dedup_hash{$seq_here})){
		$dedup_hash{$seq_here} = $dedup_hash{$seq_here}.":::".$tmp[1];
	}else{
		$dedup_hash{$seq_here} = $tmp[1];
	}
}

my @sorted_keys = sort { $dedup_hash{$a} cmp $dedup_hash{$b} } keys(%dedup_hash);

foreach my $seq (@sorted_keys){
	print OUT2 $seq."\t$dedup_hash{$seq}\n";
}



sub count_restriction{
	my $string = shift;
	my @sites = ("GCTAGC","AAGCTT","CTCGAG","GCATGC","GAGCTC","TCGCTTCGAGCA","GGTACCTGACCTC");
	foreach my $one (@sites){
		my $count = () = $string =~ /\Q$one/g;
		#print "$one $count\n";
		if ($count != 1) {
			# body...
			die "something went wrong $string $one\n";
		}
	}
}


sub rc {
  my $string = shift;
  $string =~ tr/acgtACGT/tgcaTGCA/;
  $string = reverse($string);
  return($string);
}





