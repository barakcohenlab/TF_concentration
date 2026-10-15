#!/usr/bin/perl  
#
#
use Getopt::Long;
use List::Util qw[max sum];
use warnings;
my ( $psam, $prob_cutoff, $outname, $isLogged, $ref_seq);  


my $F_pwm_file = "/scratch/bclab/yuning/libDesign/pfPilot/pwms/foxa1_jasapr_forScoring.txt";
	#$ref_seq = "AAGTAAACATGT";
my $F_ref_seq = "ACATGTTTACTT";
my $H_pwm_file = "/scratch/bclab/yuning/libDesign/pfPilot/pwms/hnf4a_jasapr_forScoring.txt";
my $H_ref_seq = "CAAAGTCCA";

my %F_psam_matrix;
my $F_pwm_length=0;

my %H_psam_matrix;
my $H_pwm_length=0;

my $probeFile = $ARGV[0];
my $out_prefix = $ARGV[1];

open(REF1,"<$F_pwm_file") || die "no $F_pwm_file\n";

while(<REF1>){
	chomp;
	my @row = split("\t");
	if($row[0] =~/^>/){
		print "$row[0]\n";
		next;
	}
	  $F_psam_matrix{"A"}[$F_pwm_length] = $row[0];
		$F_psam_matrix{"C"}[$F_pwm_length] = $row[1];
		$F_psam_matrix{"G"}[$F_pwm_length] = $row[2];
		$F_psam_matrix{"T"}[$F_pwm_length] = $row[3];
		$F_psam_matrix{"a"}[$F_pwm_length] = $row[0];
		$F_psam_matrix{"c"}[$F_pwm_length] = $row[1];
		$F_psam_matrix{"g"}[$F_pwm_length] = $row[2];
		$F_psam_matrix{"t"}[$F_pwm_length] = $row[3];
		$F_pwm_length++;
}
close(REF1);
my $F_tf_conc =foxa1_pwmToKd($F_ref_seq);

open(REF2,"<$H_pwm_file") || die "no $H_pwm_file\n";

while(<REF2>){
	chomp;
	my @row = split("\t");
	if($row[0] =~/^>/){
		print "$row[0]\n";
		next;
	}
	  $H_psam_matrix{"A"}[$H_pwm_length] = $row[0];
		$H_psam_matrix{"C"}[$H_pwm_length] = $row[1];
		$H_psam_matrix{"G"}[$H_pwm_length] = $row[2];
		$H_psam_matrix{"T"}[$H_pwm_length] = $row[3];
		$H_psam_matrix{"a"}[$H_pwm_length] = $row[0];
		$H_psam_matrix{"c"}[$H_pwm_length] = $row[1];
		$H_psam_matrix{"g"}[$H_pwm_length] = $row[2];
		$H_psam_matrix{"t"}[$H_pwm_length] = $row[3];
		$H_pwm_length++;
}
close(REF2);
my $H_tf_conc =hnf4a_pwmToKd($H_ref_seq);



open(MUT1,"</scratch/bclab/yuning/libDesign/pfPilot/pwms/0d0001_all12mer_mutateTable_foxa1.txt") || die "no file\n";
my %foxa1_mutation_lookup;
while(<MUT1>){
	chomp;
	my @row = split("\t");
	$foxa1_mutation_lookup{$row[0]}=$row[1];
}
close(MUT1);

open(MUT2,"</scratch/bclab/yuning/libDesign/pfPilot/pwms/0d0001_all12mer_mutateTable_hnf4a.txt") || die "no file\n";
my %hnf4a_mutation_lookup;
while(<MUT2>){
	chomp;
	my @row = split("\t");
	$hnf4a_mutation_lookup{$row[0]}=$row[1];
}
close(MUT2);

my %foxa1_mut_log;
my %hnf4a_mut_log;

my $f_count = length((keys %foxa1_mutation_lookup)[0]);
my $h_count = length((keys %hnf4a_mutation_lookup)[0]);

open(OUT2, ">$out_prefix".".log");
open(OUT3, ">$out_prefix"."_allInfo.txt");

open(SEQ,"$probeFile") || die "no file $probeFile\n";
while(my  $curr_probe = <SEQ>){
	chomp $curr_probe;
	my @row = split("\t",$curr_probe);
	#print "$f_count $h_count\n";
	#foreach my $i (0..length($row[4])-$f_count){
	my $probe_to_mut = $row[4];
	### F mut
	my $i = 0;
	while($i <=(length($row[4])-$f_count)){
		#print $seq_name."\t".$all_fa{$seq_name}."\t".substr($all_fa{$seq_name}, $i, 10)."\n";
		if (exists($foxa1_mutation_lookup{substr($probe_to_mut, $i, $f_count)})) {
			# body...
			my $tmp_new_kmer=$foxa1_mutation_lookup{substr($probe_to_mut, $i, $f_count)};
			substr($probe_to_mut, $i, $f_count, $tmp_new_kmer);
			#print "$row[0] $i $row[4] $probe_to_mut\n";
			$foxa1_mut_log{$row[0]}{$i}=1;
			$i=max(($i-$f_count,0));
		}else{
			$i=$i+1;
		}
	}
	my $probe_to_mut_onlyf = $probe_to_mut;
	my ($tmp_prob_f, $tmp_prob_h) = fh_prob($probe_to_mut_onlyf);
	if ($probe_to_mut_onlyf =~/GCTAGC|AAGCTT|CTCGAG|GCATGC|GAGCTC|TCGCTTCGAGCA|GGTACCTGACCTC/ || $probe_to_mut_onlyf =~/^ATGC|^TAGC|^TCGAG/ || rc($probe_to_mut_onlyf) =~ /TCGCTTCGAGCA|GGTACCTGACCTC|^ATGC|^TAGC|^TCGAG/) { 
    	print "warning: $row[0] mutF contains restriction site, ignored\n";
    	next;
	}
	print OUT3 "$row[0]_mutF\t$row[1]\t$tmp_prob_f\t$tmp_prob_h\t$probe_to_mut_onlyf\n";
	### only mutate H
	my $probe_to_mut_onlyh = $row[4];
	$i = 0;
	while($i <=(length($row[4])-$h_count)){
		#print $seq_name."\t".$all_fa{$seq_name}."\t".substr($all_fa{$seq_name}, $i, 10)."\n";
		if (exists($hnf4a_mutation_lookup{substr($probe_to_mut_onlyh, $i, $h_count)})) {
			# body...
			my $tmp_new_kmer=$hnf4a_mutation_lookup{substr($probe_to_mut_onlyh, $i, $h_count)};
			substr($probe_to_mut_onlyh, $i, $h_count, $tmp_new_kmer);
			#print "$row[0] $i $row[4] $probe_to_mut\n";
			$hnf4a_mut_log{$row[0]}{$i}=1;
			$i=max(($i-$h_count,0));
		}else{
			$i=$i+1;
		}
	}
	#print "$probe_to_mut_onlyh\n";
	if ($probe_to_mut_onlyh =~/GCTAGC|AAGCTT|CTCGAG|GCATGC|GAGCTC|TCGCTTCGAGCA|GGTACCTGACCTC/ || $probe_to_mut_onlyh =~/^ATGC|^TAGC|^TCGAG/ || rc($probe_to_mut_onlyh) =~ /TCGCTTCGAGCA|GGTACCTGACCTC|^ATGC|^TAGC|^TCGAG/) { 
    	print "$row[0] mutH contains restriction site, ignored\n";
    	next;
	}
	($tmp_prob_f, $tmp_prob_h) = fh_prob($probe_to_mut_onlyh);
	print OUT3 "$row[0]_mutH\t$row[1]\t$tmp_prob_f\t$tmp_prob_h\t$probe_to_mut_onlyh\n";
	### FH mut
	$i = 0;
	while($i <=(length($row[4])-$h_count)){
		#print $seq_name."\t".$all_fa{$seq_name}."\t".substr($all_fa{$seq_name}, $i, 10)."\n";
		if (exists($hnf4a_mutation_lookup{substr($probe_to_mut, $i, $h_count)})) {
			# body...
			my $tmp_new_kmer=$hnf4a_mutation_lookup{substr($probe_to_mut, $i, $h_count)};
			substr($probe_to_mut, $i, $h_count, $tmp_new_kmer);
			#print "$row[0] $i $row[4] $probe_to_mut\n";
			#$hnf4a_mut_log{$row[0]}{$i}=1;
			$i=max(($i-$h_count,0));
		}else{
			$i=$i+1;
		}
	}
	#print "$probe_to_mut\n";
	if ($probe_to_mut =~/GCTAGC|AAGCTT|CTCGAG|GCATGC|GAGCTC|TCGCTTCGAGCA|GGTACCTGACCTC/ || $probe_to_mut =~/^ATGC|^TAGC|^TCGAG/ || rc($probe_to_mut) =~ /TCGCTTCGAGCA|GGTACCTGACCTC|^ATGC|^TAGC|^TCGAG/) { 
    	print "$row[0] mutaFH contains restriction site, ignored\n";
    	next;
	}
	($tmp_prob_f, $tmp_prob_h) = fh_prob($probe_to_mut);
	print OUT3 "$row[0]_mutaFH\t$row[1]\t$tmp_prob_f\t$tmp_prob_h\t$probe_to_mut\n";

	print OUT2 "$row[0]:mutF:loci\t".scalar(keys %{$foxa1_mut_log{$row[0]}});
	foreach my $j (sort {$a <=> $b} (keys %{$foxa1_mut_log{$row[0]}})){
		print OUT2 "\t$j";
	}
	print OUT2 "\n$row[0] mutH loci\t".scalar(keys %{$hnf4a_mut_log{$row[0]}});
	foreach my $j (sort {$a <=> $b} (keys %{$hnf4a_mut_log{$row[0]}})){
		print OUT2 "\t$j";
	}
	print OUT2 "\n";

	#print "new prob F: $f_current_seq_final_p\t H: $h_current_seq_final_p\n";
	my ($tmpscr_prob_f, $tmpscr_prob_h) = (1,1);
	my $shuffled_seq;
	my $round_count=0;
	while($tmpscr_prob_f>0.0003 || $tmpscr_prob_h>0.0003){
		$round_count++;
		my $seed = int(rand(100)*1000);
		#my $seed =$round_count;
		$shuffled_seq = `python ~/pipelines/shuffle_oneSeq_dinu.py $row[4] $seed`;
	  chomp $shuffled_seq;
	  if ($shuffled_seq =~/GCTAGC|AAGCTT|CTCGAG|GCATGC|GAGCTC|TCGCTTCGAGCA|GGTACCTGACCTC/ || $shuffled_seq =~/^ATGC|^TAGC|^TCGAG/) { 
    	next;
		}
	  ($tmpscr_prob_f, $tmpscr_prob_h) = fh_prob($shuffled_seq);
	  #print "$seed\t$tmpscr_prob_f\t$tmpscr_prob_h\n$row[4]\n$shuffled_seq\n";
	}
	print OUT2 "shuffle takes\t$round_count round\n$row[4]\n$shuffled_seq\n";
	print OUT3 "$row[0]_scr\t$row[1]\t$tmpscr_prob_f\t$tmpscr_prob_h\t$shuffled_seq\n";
}
close(SEQ);
close(OUT2);
close(OUT3);


sub fh_prob{
	### F new prob
	my $sequence= shift;
	my $current_seq_p=1;
	foreach my $i (0..length($sequence)-$f_count){
		#print substr($probe_to_mut, $i, $f_count)."\n";
		my $tmp_score = foxa1_pwmToKd(substr($sequence, $i, $f_count));
		$current_seq_p = $current_seq_p * (1-1/(1+$tmp_score/$F_tf_conc));
	}
	my $f_current_seq_final_p = 1-$current_seq_p;
	### H new prob
	$current_seq_p=1;
	foreach my $i (0..length($sequence)-$h_count){
		#print substr($probe_to_mut, $i, $f_count)."\n";
		my $tmp_score = hnf4a_pwmToKd(substr($sequence, $i, $h_count));
		$current_seq_p = $current_seq_p * (1-1/(1+$tmp_score/$H_tf_conc));
	}
	my $h_current_seq_final_p = 1-$current_seq_p;
	return($f_current_seq_final_p, $h_current_seq_final_p);

}


sub foxa1_pwmToKd{
	my $sequence= shift;
	my $rc_seq = rc($sequence);
	my @individual_bases = split(//, $sequence);
	my @rc_individual_bases = split(//, $rc_seq);
	my $score =1;
	my $rc_score = 1;
	foreach my $i (0..$#individual_bases){
		$score = $score * $F_psam_matrix{$individual_bases[$i]}[$i]/0.25;
		#print "$individual_bases[$i] $i $score\n";
	}
	foreach my $i (0..$#rc_individual_bases){
		$rc_score = $rc_score * $F_psam_matrix{$rc_individual_bases[$i]}[$i]/0.25;
	}
	my $final_score =1/(max(($score,$rc_score)));
	return($final_score);
}

sub hnf4a_pwmToKd{
	my $sequence= shift;
	my $rc_seq = rc($sequence);
	my @individual_bases = split(//, $sequence);
	my @rc_individual_bases = split(//, $rc_seq);
	my $score =1;
	my $rc_score = 1;
	foreach my $i (0..$#individual_bases){
		$score = $score * $H_psam_matrix{$individual_bases[$i]}[$i]/0.25;
	}
	foreach my $i (0..$#rc_individual_bases){
		$rc_score = $rc_score * $H_psam_matrix{$rc_individual_bases[$i]}[$i]/0.25;
	}
	my $final_score =1/(max(($score,$rc_score)));
	return($final_score);
}


sub rc {
  my $string = shift;
  $string =~ tr/acgtACGT/tgcaTGCA/;
  $string = reverse($string);
  return($string);
}