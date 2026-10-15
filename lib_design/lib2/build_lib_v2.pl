#!/usr/bin/perl

###

use List::Util qw[max shuffle];
use strict;
#use warnings;

my $inputSeqFile = $ARGV[0];
#my $out_prefix = $ARGV[1];
my ($out_prefix) = $inputSeqFile =~ m/clean_cat_(.*)_allInfo_all\.txt$/;


open(USED,"</scratch/bclab/yuning/sequencing/pft_v1_dict_name_bc_cre.txt");

my %used_barcodes;

while(my $line=<USED>){
	chomp $line;
	my @cols = split(/\t/, $line);
	$used_barcodes{$cols[1]}=1;
}
close(USED);


my @barcodes_array;
open(BAR,"</scratch/bclab/yuning/libDesign/pfPilot/workplace/myBarcodes_10mer.txt");


### filter used barcodes

while(my $oneBarcode=<BAR>){
	chomp $oneBarcode;
	if (exists($used_barcodes{$oneBarcode})) {
		# body...
		next;
	}
	push(@barcodes_array, $oneBarcode);
#	print "$oneBarcode\n";
}
close(BAR);

@barcodes_array = shuffle(@barcodes_array);



open(COORD,$inputSeqFile) || die "no such file $inputSeqFile\n";


open(OUT1, ">$out_prefix"."_agilent_submission.txt");
open(OUT2, ">$out_prefix"."_agilent_alldata.txt");
#open(OUT3, ">$out_prefix"."_allInfo.txt");


#my $test_Seq = "GGTACCTGACCTCGCTAGCTATTTGTCTTGCTTGGAGGATCAGGGGTCAGGCCCTGTCCGTGACCTAGTTCAGGTTAAATAAAGCTGAGCTGGGAGGCTGTTTACTGCCTGTAGACCGCTGAGCAGAGTCCATGGCCCTGCCGGGCTGCCCTGGCTGGGGGCAGGAGCACGCCTGGGCATGGGGAAGCTTCTCGAGGCATGCACAGTTATGCGAGCTCTCGCTTCGAGCA";
#my $test_Seq = "GAGCTCGGTACCTGACCTCGCTAGCTATTTGTCTTGCTTGGAGGATCAGGGGTCAGGCCCTGTCCGTGACCTAGTTCAGGTTAAATAAAGCTGAGCTGGGAGGCTGTTTACTGCCTGTAGACCGCTGAGCAGAGTCCATGGCCCTGCCGGGCTGCCCTGGCTGGGGGCAGGAGCACGCCTGGGCATGGGGAAGCTTCTCGAGGCATGCACAGTTATGCGAGCTCTCGCTTCGAGCA";
#count_restriction($test_Seq);


### basal
my @basal_probes = makeBasal(101);
foreach my $x (0..$#basal_probes) {
	#print OUT1 "$probeseq1[$x]\n";
	count_restriction($basal_probes[$x]);
	print OUT1 "$basal_probes[$x]\n";
	print OUT2 "$basal_probes[$x]\t"."v2_basal_hsp68_rep".($x+1)."\n";
}

### others
while (my $curr_rna = <COORD>) {
    chomp $curr_rna;
    my @tmp = split /\t/, $curr_rna;
    if ((scalar @tmp) != 2 || length($tmp[0]) != 165) {
    	# body...
    	die "something went wrong: @tmp\n";
    }
    my $bc_num = 3;

    my @probeseq1 = makeProbe($tmp[0],$bc_num);
	foreach my $x (0..$#probeseq1) {
		#print OUT1 "$probeseq1[$x]\n";
		count_restriction($probeseq1[$x]);
		print OUT1 "$probeseq1[$x]\n";
		print OUT2 "$probeseq1[$x]\t"."$tmp[1]"."_rep".($x+1)."\n";
	}

}

sub makeBasal{
	my ($numOfBC) = @_;
	#print "$CREseq\t$numOfBC\n";
	my @probe_array;
	foreach my $i (1..$numOfBC) {
		my $tmp_bc = pop @barcodes_array;
		#print "$i\t$tmp_bc\n";
		#$probe_array[$i-1] = "GGTACCTGACCTCGCTAGC ".$CREseq." AAGCTTCTCGAGGCATGC ".$tmp_bc." GAGCTCTCGCTTCGAGCA";
		$probe_array[$i-1] = "GGTACCTGACCTCGCTAGCAAGCTTACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTCGAGGCATGC".$tmp_bc."GAGCTCTCGCTTCGAGCA";
	}
	return(@probe_array);
}


sub makeProbe{
	my ($CREseq, $numOfBC) = @_;
	#print "$CREseq\t$numOfBC\n";
	my @probe_array;
	foreach my $i (1..$numOfBC) {
		my $tmp_bc = pop @barcodes_array;
		#print "$i\t$tmp_bc\n";
		#$probe_array[$i-1] = "GGTACCTGACCTCGCTAGC ".$CREseq." AAGCTTCTCGAGGCATGC ".$tmp_bc." GAGCTCTCGCTTCGAGCA";
		$probe_array[$i-1] = "GGTACCTGACCTCGCTAGC".$CREseq."AAGCTTCTCGAGGCATGC".$tmp_bc."GAGCTCTCGCTTCGAGCA";
	}
	return(@probe_array);
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





