#!/usr/bin/perl
#

use Getopt::Long;
use List::Util qw[max sum];
use warnings;



my $motif_file = "/ref/bclab/software/yuning/refs/jaspar_2024/human755_20240529015908_JASPAR2024_combined_matrices_615921_meme.txt";

open(REF1,"<$motif_file") || die "no $motif_file\n";

my %motif_hash;
my $tmp_pos_count = 0;
my $tmp_tf;
my %tf_names_hash;

while(my $line = <REF1>){
	chomp $line;
	if($line =~ /^MOTIF\s(.*)\s.*\.(.*)$/){
		$tmp_tf = $1;
		$tf_names_hash{$2}=$tmp_tf;
		#print $tmp_tf."\t$2\n";
		$tmp_pos_count = 0;
	}
	if($line =~ /^\s/){
		my @row=split(" ",$line);
		$motif_hash{$tmp_tf}{"A"}[$tmp_pos_count] = $row[0];
		$motif_hash{$tmp_tf}{"C"}[$tmp_pos_count] = $row[1];
		$motif_hash{$tmp_tf}{"G"}[$tmp_pos_count] = $row[2];
		$motif_hash{$tmp_tf}{"T"}[$tmp_pos_count] = $row[3];
		$motif_hash{$tmp_tf}{"a"}[$tmp_pos_count] = $row[0];
		$motif_hash{$tmp_tf}{"c"}[$tmp_pos_count] = $row[1];
		$motif_hash{$tmp_tf}{"g"}[$tmp_pos_count] = $row[2];
		$motif_hash{$tmp_tf}{"t"}[$tmp_pos_count] = $row[3];

		$tmp_pos_count++;
		#print $motif_hash{$tmp_tf}[$tmp_pos_count][3]."\n";

	}
	#if($line =~ /^URL/){
		#$motif_hash{$tmp_tf} = \@tmp_matrix;
		#}
}


close(REF1);


my $sequence_file = "/scratch/bclab/yuning/sequencing/noRep_padded_pft_v1_name_bc.txt";
my %seq_hash;
open(REF3,"<$sequence_file") || die "no $sequence_file\n";
while(my $line = <REF3>){
	chomp $line;
	my @row=split("\t",$line);
	$seq_hash{$row[0]} = $row[1];
	#print $row[0]."\t".$row[1]."\n";
	#print @{$motif_hash{$row[0]}{"A"}};
}
close(REF3);

my $rna_expr = "/scratch/bclab/yuning/libDesign/tfValidate/fimoR5_pft_v1/just_humanTFs_k1f_r3_RNAgeneCount.txt";
my %expr_hash;
my %best_motif_hash;
open(REF5,"<$rna_expr") || die "no $rna_expr\n";
while(my $line = <REF5>){
	chomp $line;
	my @row=split("\t",$line);
	if (! ( exists($tf_names_hash{$row[1]}) ) ) {
		next;
		# body...
	}
	my $jaspar_id = $tf_names_hash{$row[1]};
	if (exists($expr_hash{$jaspar_id})) {
		if ( $row[9] > $expr_hash{$jaspar_id}) {
			$expr_hash{$jaspar_id} = $row[9]; 
		}
	}else{
		$expr_hash{$jaspar_id} = $row[9]; ## basal expr
	}
	#print $row[0]."\t".$row[1]."\n";
	#print @{$motif_hash{$row[0]}{"A"}};
	my $motif_length = scalar @{$motif_hash{$jaspar_id}{"A"}};
	my $best_motif = "";
	foreach my $i (0..($motif_length-1)){
		my ($score, $base) = max_pos($jaspar_id,$i);
		$best_motif = $best_motif.$base;
	}
	$best_motif_hash{$jaspar_id} = $best_motif;
}
close(REF5);




my @tfbs_list = sort { $expr_hash{$b} <=> $expr_hash{$a} } keys(%expr_hash);
my @tf_picked = @tfbs_list[0..106];
unshift @tf_picked, 'MA1494.2';
$best_motif_hash{'MA1494.2'} = 'CAAAGTCCA';

unshift @tf_picked, 'MA0148.5';
$best_motif_hash{'MA0148.5'} = 'ATGTTTAC';

#print $tfbs_number."\n";
#print join(",", @tf_picked)."\n";


### perl -le '{$a{"00"}="A";$a{"01"}="T";$a{"10"}="C";$a{"11"}="G";$n=5;for(0..(2**$n-1)){print join"::",map{$_} unpack("a1"x$n,sprintf("%0*b",$n,$_))}}' 


#my $screen_context ="CGAcgtgcttctgtcactatcagtatcgggcatgacccagccagatccgcacaggccctatttcgtgGaaaaggggacacaatgccgcgcaactggcTGAatcgaagTCTCCAAATcgtgcatagcttgtaattggGACAgcgagggttacgcccgGGAAGCAAT";
#my $screen_context2 ="TCAAactaagagaCTAAAAAAGGcatgagttgACGCTggcgagaacctccgTGTCACCgagcctcagcggctcccggttgcccacgaggcgagaggaaaaagcaGATGCTgcgggcggcgggaaaagcggagcccaaaggCTgcgagggcgggCAAGTCTGCGGC";

### try avoiding RE site
my $screen_context ="CGAcgtgcttctgtcactatcagtatcgggcatgacccagccagatccgcacaggccctatttcgtgGaaaaggggacacaatgccgcgcaactggcTGAatcgaagTCTCCAAATcgtgcatagcttgtaattggGACAgcgaCggttacgcccgGGAAGCAAT";
my $screen_context2 ="TCAAactaagagaCTAAAAAAGGcatgagttgACGCTggcgagaacctccgTGTCACCgagcctcagcggctGccggttCGccacgaggcgagaggaaaaagcaGATGCTgcgggcggcgggaTaTCcggagcccaaaggCTAcgagggcgggCAAGTCTGCGGC";

#$screen_context = uc($screen_context);
#$screen_context2 = uc($screen_context2);

my @all_context = ($screen_context, $screen_context2);

my $output_prefix_name = "tfscreen";
open(OUT3, ">$output_prefix_name"."_allInfo.txt");

my $F_site = "ATGTTTAC";
my $H_site = "CAAAGTCCA";

my @plugin_sites2 = (37,57,82,112,137);
foreach my $site (0..$#plugin_sites2){
	$plugin_sites2[$site] = $plugin_sites2[$site] - 9;
}
# print join(" ",@plugin_sites2)."\n";
### ++--+
### xFxFx


### make empty probe (just F/H, no other TFs)


foreach my $con (0..$#all_context){
	my $start_context = uc($all_context[$con]);
	my $focus_tf = $F_site;
	substr($start_context, $plugin_sites2[1], length($focus_tf), lc($focus_tf));
	substr($start_context, $plugin_sites2[3], length($focus_tf), lc(rc($focus_tf)));
	print OUT3 "$start_context\tscreenF_con$con"."_empty\n";

	$start_context = uc($all_context[$con]);
	$focus_tf = $H_site;
	substr($start_context, $plugin_sites2[1], length($focus_tf), lc($focus_tf));
	substr($start_context, $plugin_sites2[3], length($focus_tf), lc(rc($focus_tf)));
	print OUT3 "$start_context\tscreenH_con$con"."_empty\n";

}




foreach my $tf (@tf_picked){
	#print "$tf\t$expr_hash{$tf}\t$best_motif_hash{$tf}\n";
	foreach my $con (0..$#all_context){
		my $start_context = uc($all_context[$con]);
		my $othertf_bs = $best_motif_hash{$tf};
		my $focus_tf = $F_site;
		substr($start_context, $plugin_sites2[0], length($othertf_bs), lc($othertf_bs));
		substr($start_context, $plugin_sites2[1], length($focus_tf), lc($focus_tf));
		substr($start_context, $plugin_sites2[2], length($othertf_bs), lc(rc($othertf_bs)));
		substr($start_context, $plugin_sites2[3], length($focus_tf), lc(rc($focus_tf)));
		substr($start_context, $plugin_sites2[4], length($othertf_bs), lc($othertf_bs));
		print OUT3 "$start_context\tscreenF_con$con"."_"."$tf\n";

		$start_context = uc($all_context[$con]);
		$othertf_bs = $best_motif_hash{$tf};
		$focus_tf = $H_site;
		substr($start_context, $plugin_sites2[0], length($othertf_bs), lc($othertf_bs));
		substr($start_context, $plugin_sites2[1], length($focus_tf), lc($focus_tf));
		substr($start_context, $plugin_sites2[2], length($othertf_bs), lc(rc($othertf_bs)));
		substr($start_context, $plugin_sites2[3], length($focus_tf), lc(rc($focus_tf)));
		substr($start_context, $plugin_sites2[4], length($othertf_bs), lc($othertf_bs));
		print OUT3 "$start_context\tscreenH_con$con"."_"."$tf\n";

	}

}
close(OUT3);


foreach my $tf (@tfbs_list){
	#print "$tf\t$rna_expr{$tf}\n";
}




### print consensus and mutated TFBS;

#my ($t1, $t2) = max_pos("MA1142.2",3);
#print $t1.$t2;

sub return_maxPWM_loci{
	my ($tf, $sequence)= @_;
	my @individual_bases = split(//, $sequence);
	my $score =0;
	my $position=-1;
	foreach my $i (0..$#individual_bases){
		if ($motif_hash{$tf}{$individual_bases[$i]}[$i] > $score) {
			# body...
			$score = $motif_hash{$tf}{$individual_bases[$i]}[$i];
			$position = $i;
		}
		#print "$individual_bases[$i] $motif_hash{$tf}{$individual_bases[$i]}[$i] $score\n";
	}
	return(($score, $position));
}

sub min_pos {
	my ($tf, $pos) = @_;
#	print $pos."\t".$tf."\n";
	my $score = 111110;
	my $base = "N";
	foreach my $i ("A","C","G","T"){
		if ($motif_hash{$tf}{$i}[$pos] < $score) {
			# body...
			$base = $i;
			$score = $motif_hash{$tf}{$i}[$pos];
		}
		#print $i."\n";
	}
	return(($score, $base));
}

sub max_pos {
	my ($tf, $pos) = @_;
#	print $pos."\t".$tf."\n";
	my $score = 0;
	my $base = "N";
	foreach my $i ("A","C","G","T"){
		if ($motif_hash{$tf}{$i}[$pos] > $score) {
			# body...
			$base = $i;
			$score = $motif_hash{$tf}{$i}[$pos];
		}
		#print $i."\n";
	}
	return(($score, $base));
}


sub consensus {
	my $tf = shift;
	my $c_motif = "";
	foreach my $pos (0..(scalar @{$motif_hash{$tf}{"A"}} - 1)){
		my ($score, $base) = max_pos($tf,$pos);
		if ($score > 0.5) {
			$c_motif = $c_motif.$base;
		}else{
			$c_motif = $c_motif."-";
		}
		#print $row.",".$motif_hash{$tf}{"c"}[$row]."\t";
	}
	return($c_motif);

}


sub rc {
  my $string = shift;
  $string =~ tr/acgtACGT/tgcaTGCA/;
  $string = reverse($string);
  return($string);
}
