#!/usr/bin/perl  
#
#
use Getopt::Long;
use List::Util qw[max sum];
use warnings;
my ( $psam, $prob_cutoff, $outname, $isLogged, $ref_seq);  


#my $F_pwm_file = "/scratch/bclab/yuning/libDesign/pfPilot/pwms/foxa1_jasapr_forScoring.txt";
my $F_pwm_file = "/scratch/bclab/yuning/cutTagModel/0523_redo_foxa1/foxa1_pwm_nt9_v8v_v8_1e3_0611.txt";
	#$ref_seq = "AAGTAAACATGT";
my $F_ref_seq = "ATGTTTAC";
my $H_pwm_file = "/scratch/bclab/yuning/cutTagModel/0523_redo_hnf4a/hnf4a_pwm_nt15_v8v_v8_0611.txt";
my $H_ref_seq = "CAAAGTCCA";

my %F_psam_matrix;
my $f_count=0;

my %H_psam_matrix;
my $h_count=0;

my $probeFile = "/scratch/bclab/yuning/libDesign/tfValidate/eachCRE_tfbs_list/spacer_scramble_context.txt";
#my $out_prefix = $ARGV[1];

open(REF1,"<$F_pwm_file") || die "no $F_pwm_file\n";

while(<REF1>){
	chomp;
	my @row = split("\t");
	if($row[0] =~/^>/){
		print "$row[0]\n";
		next;
	}
	  $F_psam_matrix{"A"}[$f_count] = $row[0];
		$F_psam_matrix{"C"}[$f_count] = $row[1];
		$F_psam_matrix{"G"}[$f_count] = $row[2];
		$F_psam_matrix{"T"}[$f_count] = $row[3];
		$F_psam_matrix{"a"}[$f_count] = $row[0];
		$F_psam_matrix{"c"}[$f_count] = $row[1];
		$F_psam_matrix{"g"}[$f_count] = $row[2];
		$F_psam_matrix{"t"}[$f_count] = $row[3];
		$f_count++;
}
close(REF1);
#my $F_tf_conc =foxa1_pwmToKd($F_ref_seq);
my $F_tf_conc = 159.465;

open(REF2,"<$H_pwm_file") || die "no $H_pwm_file\n";

while(<REF2>){
	chomp;
	my @row = split("\t");
	if($row[0] =~/^>/){
		print "$row[0]\n";
		next;
	}
	  $H_psam_matrix{"A"}[$h_count] = $row[0];
		$H_psam_matrix{"C"}[$h_count] = $row[1];
		$H_psam_matrix{"G"}[$h_count] = $row[2];
		$H_psam_matrix{"T"}[$h_count] = $row[3];
		$H_psam_matrix{"a"}[$h_count] = $row[0];
		$H_psam_matrix{"c"}[$h_count] = $row[1];
		$H_psam_matrix{"g"}[$h_count] = $row[2];
		$H_psam_matrix{"t"}[$h_count] = $row[3];
		$h_count++;
}
close(REF2);
#my $H_tf_conc =hnf4a_pwmToKd($H_ref_seq);
my $H_tf_conc = 14.69;


#my %foxa1_mut_log;
#my %hnf4a_mut_log;

#my $f_count = length((keys %foxa1_mutation_lookup)[0]);
#my $h_count = length((keys %hnf4a_mutation_lookup)[0]);

my $out_prefix = "spacer_scramble_set";

#open(OUT2, ">$out_prefix".".log");
open(OUT3, ">$out_prefix"."_allInfo.txt");

my %context_list;

open(SEQ,"$probeFile") || die "no file $probeFile\n";
while(my  $curr_probe = <SEQ>){
	chomp $curr_probe;
	my @row = split("\t",$curr_probe);
	#print "$f_count $h_count\n";
	#foreach my $i (0..length($row[4])-$f_count){
	my $probe_to_mut = $row[0];
	my $name_prefix = $row[1];
	### F mut
	$context_list{"ss_ori_".$name_prefix} = $probe_to_mut;

	my ($tmpscr_prob_f, $tmpscr_prob_h) = (1,1);
	my $shuffled_seq;
	my $round_count=0;
	my $probes_generated_count = 0;
	my ($ori_prob_f, $ori_prob_h) = fh_prob($probe_to_mut);
	print "$name_prefix original prob: $ori_prob_f $ori_prob_h\n";
	if ($probe_to_mut =~/GCTAGC|AAGCTT|CTCGAG|GCATGC|GAGCTC|TCGCTTCGAGCA|GGTACCTGACCTC/ || $probe_to_mut =~/^ATGC|^TAGC|^TCGAG/ || rc($probe_to_mut) =~ /TCGCTTCGAGCA|GGTACCTGACCTC|^ATGC|^TAGC|^TCGAG/) { 
    die "wrong wrong wrong\n\n\n";
	}
	while($probes_generated_count <4){
		$round_count++;
		my $seed = int(rand(100)*1000);
		#my $seed =$round_count;
		$shuffled_seq = `python ~/pipelines/shuffle_oneSeq_dinu.py $probe_to_mut $seed`;
	  chomp $shuffled_seq;
	  if ($shuffled_seq =~/GCTAGC|AAGCTT|CTCGAG|GCATGC|GAGCTC|TCGCTTCGAGCA|GGTACCTGACCTC/ || $shuffled_seq =~/^ATGC|^TAGC|^TCGAG/ || rc($shuffled_seq) =~ /TCGCTTCGAGCA|GGTACCTGACCTC|^ATGC|^TAGC|^TCGAG/) { 
    	next;
		}
	  ($tmpscr_prob_f, $tmpscr_prob_h) = fh_prob($shuffled_seq);
		if ($tmpscr_prob_f>$ori_prob_f*1 || $tmpscr_prob_h>$ori_prob_h*1) {
			# body...
			next;
		}
		$probes_generated_count++;
	  #print "$round_count\t$tmpscr_prob_f\t$tmpscr_prob_h\n";
	  $context_list{"ss_scr".$probes_generated_count."_".$name_prefix} = $shuffled_seq;
	  #print "$round_count\t$probes_generated_count\t$tmpscr_prob_f\t$tmpscr_prob_h\n$probe_to_mut\n$shuffled_seq\n";
	}
#	print OUT2 "shuffle takes\t$round_count round\n$row[4]\n$shuffled_seq\n";
#print  "$round_count\t$row[0]_scr\t$row[1]\t$tmpscr_prob_f\t$tmpscr_prob_h\t$shuffled_seq\n";
}

my @plugin_sites1 = (67,87,107);
my @plugin_sites2 = (62,87,112);
my @plugin_sites3 = (47,87,135);

my $F_site = "ATGTTTAC";
my $H_site = "CAAAGTCCA";

foreach my $name (sort keys %context_list){
	#print "$name\t$context_list{$name}\n";

	my @tmp_site = @plugin_sites1;
	my $start_cre_ff = $context_list{$name};
	my $start_cre_hh = $context_list{$name};
	my $start_cre_fh = $context_list{$name};
	foreach my $site (@tmp_site){
		substr($start_cre_ff, $site-9, length($F_site), lc($F_site));
		substr($start_cre_hh, $site-9, length($H_site), lc($H_site));
	}
	substr($start_cre_fh, $tmp_site[0]-9, length($F_site), lc($F_site));
	substr($start_cre_fh, $tmp_site[1]-9, length($H_site), lc($H_site));
	substr($start_cre_fh, $tmp_site[2]-9, length($F_site), lc($F_site));
	print OUT3 "$start_cre_ff\tfff1_$name\n";
	print OUT3 "$start_cre_hh\thhh1_$name\n";
	print OUT3 "$start_cre_fh\tfhf1_$name\n";

	@tmp_site = @plugin_sites2;
	$start_cre_ff = $context_list{$name};
	$start_cre_hh = $context_list{$name};
	$start_cre_fh = $context_list{$name};
	foreach my $site (@tmp_site){
		substr($start_cre_ff, $site-9, length($F_site), lc($F_site));
		substr($start_cre_hh, $site-9, length($H_site), lc($H_site));
	}
	substr($start_cre_fh, $tmp_site[0]-9, length($F_site), lc($F_site));
	substr($start_cre_fh, $tmp_site[1]-9, length($H_site), lc($H_site));
	substr($start_cre_fh, $tmp_site[2]-9, length($F_site), lc($F_site));
	print OUT3 "$start_cre_ff\tfff2_$name\n";
	print OUT3 "$start_cre_hh\thhh2_$name\n";
	print OUT3 "$start_cre_fh\tfhf2_$name\n";

	@tmp_site = @plugin_sites3;
	$start_cre_ff = $context_list{$name};
	$start_cre_hh = $context_list{$name};
	$start_cre_fh = $context_list{$name};
	foreach my $site (@tmp_site){
		substr($start_cre_ff, $site-9, length($F_site), lc($F_site));
		substr($start_cre_hh, $site-9, length($H_site), lc($H_site));
	}
	substr($start_cre_fh, $tmp_site[0]-9, length($F_site), lc($F_site));
	substr($start_cre_fh, $tmp_site[1]-9, length($H_site), lc($H_site));
	substr($start_cre_fh, $tmp_site[2]-9, length($F_site), lc($F_site));
	print OUT3 "$start_cre_ff\tfff3_$name\n";
	print OUT3 "$start_cre_hh\thhh3_$name\n";
	print OUT3 "$start_cre_fh\tfhf3_$name\n";


}


close(SEQ);
#close(OUT2);
close(OUT3);


sub fh_prob{
	### F new prob
	my $sequence= shift;
	my $current_seq_p=0;
	foreach my $i (0..length($sequence)-$f_count){
		#print substr($probe_to_mut, $i, $f_count)."\n";
		my $tmp_score = foxa1_pwmToKa(substr($sequence, $i, $f_count));
		if ($tmp_score > $current_seq_p) {
			$current_seq_p = $tmp_score;
			# body...
		}
	}
	my $f_current_seq_final_p = ($F_tf_conc*$current_seq_p)/(1+$F_tf_conc*$current_seq_p);
	### H new prob
	$current_seq_p=0;
	foreach my $i (0..length($sequence)-$h_count){
		#print substr($probe_to_mut, $i, $f_count)."\n";
		my $tmp_score = hnf4a_pwmToKa(substr($sequence, $i, $h_count));
		if ($tmp_score > $current_seq_p) {
			$current_seq_p = $tmp_score;
			# body...
		}
	}
	my $h_current_seq_final_p = ($H_tf_conc*$current_seq_p)/(1+$H_tf_conc*$current_seq_p);
	return($f_current_seq_final_p, $h_current_seq_final_p);

}

##### modified to remove .25 normalization


sub foxa1_pwmToKa{
	my $sequence= shift;
	my $rc_seq = rc($sequence);
	my @individual_bases = split(//, $sequence);
	my @rc_individual_bases = split(//, $rc_seq);
	my $score =1;
	my $rc_score = 1;
	foreach my $i (0..$#individual_bases){
		$score = $score * $F_psam_matrix{$individual_bases[$i]}[$i];
		#print "$individual_bases[$i] $i $score\n";
	}
	foreach my $i (0..$#rc_individual_bases){
		$rc_score = $rc_score * $F_psam_matrix{$rc_individual_bases[$i]}[$i];
	}
	my $final_score =(max(($score,$rc_score)));
	return($final_score);
}

sub hnf4a_pwmToKa{
	my $sequence= shift;
	my $rc_seq = rc($sequence);
	my @individual_bases = split(//, $sequence);
	my @rc_individual_bases = split(//, $rc_seq);
	my $score =1;
	my $rc_score = 1;
	foreach my $i (0..$#individual_bases){
		$score = $score * $H_psam_matrix{$individual_bases[$i]}[$i];
	}
	foreach my $i (0..$#rc_individual_bases){
		$rc_score = $rc_score * $H_psam_matrix{$rc_individual_bases[$i]}[$i];
	}
	my $final_score =(max(($score,$rc_score)));
	return($final_score);
}


sub rc {
  my $string = shift;
  $string =~ tr/acgtACGT/tgcaTGCA/;
  $string = reverse($string);
  return($string);
}