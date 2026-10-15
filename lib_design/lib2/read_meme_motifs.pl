#!/usr/bin/perl
#

use Getopt::Long;
use List::Util qw[max sum];
use warnings;


my $cre_motif_list = $ARGV[0];
my ($output_prefix_name) =  $cre_motif_list =~ m/curated_(otf[0-9]*_.*)\.txt/;
my ($output_prefix) =  $output_prefix_name =~ m/otf[0-9]*_(.*)$/;


my $motif_file = "/ref/bclab/software/yuning/refs/jaspar_2024/human755_20240529015908_JASPAR2024_combined_matrices_615921_meme.txt";

open(REF1,"<$motif_file") || die "no $motif_file\n";

my %motif_hash;
my $tmp_pos_count = 0;
my $tmp_tf;

while(my $line = <REF1>){
	chomp $line;
	if($line =~ /^MOTIF\s(.*)\s/){
		$tmp_tf = $1;
		#print $tmp_tf."\n";
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


my %motif_dict;
my $ori_full_seq;
my $tfbs_number = 0;


open(REF2,"<$cre_motif_list") || die "no $cre_motif_list\n";
while(my $line = <REF2>){
	chomp $line;
	my @row=split("\t",$line);
	my $offset = $row[3]-1;
	#$ori_full_seq = substr($seq_hash{$row[2]},8,165);
	$ori_full_seq = $seq_hash{$row[2]};
	#print "$ori_full_seq\n";
	if ($row[5] eq '+') {
		#print substr($seq_hash{$row[2]},$offset,$row[4]-$offset)." aaaaa\n";
	}else{
		#print rc(substr($seq_hash{$row[2]},$offset,$row[4]-$offset))." bbbbb\n";
	}
	my $tf_here = $row[0];
	my @ori_motif = split("",$row[9]);
	my @best_motif = @ori_motif;
	my @mut_motif = @ori_motif;
	my @reallymut_motif = @ori_motif;

	$motif_dict{$tfbs_number}{"ss"} = $offset;
	$motif_dict{$tfbs_number}{"len"} = $row[4]-$offset;
	$motif_dict{$tfbs_number}{"ori"} = $row[9];
	$motif_dict{$tfbs_number}{"strand"} = $row[5];
	$motif_dict{$tfbs_number}{"tf"} = $row[1];;
	$motif_dict{$tfbs_number}{"tfshort"} = substr($row[1], 0, 1);


	### Make consensus motif
	foreach my $i (0..$#ori_motif){
		my ($score, $base) = max_pos($tf_here,$i);
		if ($score > 0.5) {
			$best_motif[$i] = $base;
		}
	}
	if ($row[5] eq '+') {
		$motif_dict{$tfbs_number}{"best"} = join("",@best_motif);
	}else{
		$motif_dict{$tfbs_number}{"best"} = rc(join("",@best_motif));
	}


	### Make dead motif
	my $fold_change = 1;
	while($fold_change>0.001){
		my ($ori_s, $mut_pos) = return_maxPWM_loci($tf_here, join("",@mut_motif));
		my ($mut_s, $mut_base) = min_pos($tf_here, $mut_pos);
		$mut_motif[$mut_pos] = $mut_base;
		$fold_change = ($mut_s+1e-2)/($ori_s+1e-2) * $fold_change;
		#print "$tf_here $row[9] ".join("",@mut_motif)." $mut_pos $ori_s $mut_base $mut_s $fold_change\n";
	}
	if ($row[5] eq '+') {
		$motif_dict{$tfbs_number}{"mut"} = join("",@mut_motif);
	}else{
		$motif_dict{$tfbs_number}{"mut"} = rc(join("",@mut_motif));
	}


	### Make really dead motif (3/4 point mutations)
	my $mutation_count = 1;
	my $fold_change2 = 1;
	my $total_muts = 3;
	if ($motif_dict{$tfbs_number}{"len"} > 8) {
		# body...
		$total_muts = 4;
	}
	while($mutation_count<($total_muts+1)){
		my ($ori_s, $mut_pos) = return_maxPWM_loci($tf_here, join("",@reallymut_motif));
		my ($mut_s, $mut_base) = min_pos($tf_here, $mut_pos);
		$reallymut_motif[$mut_pos] = $mut_base;
		$fold_change2 = ($mut_s+1e-2)/($ori_s+1e-2) * $fold_change2;
		$mutation_count++;
		#print "$tf_here $row[9] ".join("",@reallymut_motif)." $mut_pos $ori_s $mut_base $mut_s $mutation_count $fold_change2\n";
	}
	if ($row[5] eq '+') {
		$motif_dict{$tfbs_number}{"rmut"} = join("",@reallymut_motif);
	}else{
		$motif_dict{$tfbs_number}{"rmut"} = rc(join("",@reallymut_motif));
	}


	#print "$row[9]\t".join("",@best_motif)."\n";
	#print $row[0]."\t".$row[1]."\t".consensus($row[0])."\t".$row[5]."\n";
	#print @{$motif_hash{$row[0]}{"A"}};
	$tfbs_number = $tfbs_number+1;

}

my @tfbs_list = sort { $motif_dict{$a}{"ss"} <=> $motif_dict{$b}{"ss"} } keys(%motif_dict);
#print $tfbs_number."\n";
#print join(",", @tfbs_list)."\n";

### perl -le '{$a{"00"}="A";$a{"01"}="T";$a{"10"}="C";$a{"11"}="G";$n=5;for(0..(2**$n-1)){print join"::",map{$_} unpack("a1"x$n,sprintf("%0*b",$n,$_))}}' 


print "$output_prefix_name    $output_prefix\n";
open(OUT3, ">$output_prefix_name"."_allInfo.txt");


foreach my $n (0..(2**$tfbs_number-1)){
	my @curr_combination = unpack("a1"x$tfbs_number,sprintf("%0*b",$tfbs_number,$n));
	#print join("::",@curr_combination)."\n";
	my $overall_manipulations = sum(@curr_combination);
	if ($tfbs_number > 8 && $overall_manipulations > 1 && $overall_manipulations < $tfbs_number) {
		# body...
		next;
	}
	my $edit_full_best = $ori_full_seq;
	my $edit_full_mut = $ori_full_seq;

	my $seq_name = $output_prefix.";".join("",@curr_combination);

	foreach my $count (0..$#curr_combination){
		if ($curr_combination[$count] == 1) {
				$seq_name = $seq_name.";".$motif_dict{$tfbs_list[$count]}{"tf"}.";".$tfbs_list[$count];
				substr($edit_full_best, $motif_dict{$tfbs_list[$count]}{"ss"}, $motif_dict{$tfbs_list[$count]}{"len"}, lc($motif_dict{$tfbs_list[$count]}{"best"}));
				substr($edit_full_mut, $motif_dict{$tfbs_list[$count]}{"ss"}, $motif_dict{$tfbs_list[$count]}{"len"}, lc($motif_dict{$tfbs_list[$count]}{"rmut"}));
		}
	}
	my $edit_full_best_165 = substr($edit_full_best,8,165);
	my $edit_full_mut_165 = substr($edit_full_mut,8,165);

	print OUT3 "$edit_full_mut_165\totf_mut_$seq_name\n";
	print OUT3 "$edit_full_best_165\totf_con_$seq_name\n";
}



foreach my $tf (keys %motif_dict){
	foreach my $j (sort keys %{$motif_dict{$tf}}){
		#print "$tf $j $motif_dict{$tf}{$j}\n";
	}
}


foreach my $tf (keys %motif_hash){
	#print $tf."\t".scalar(@{$motif_hash{$tf}{"A"}})."\n";
	#print $tf."\t".consensus($tf)."\n";
	#foreach my $row (0..(scalar @{$motif_hash{$tf}{"A"}} - 1)){
		#print $row.",".$motif_hash{$tf}{"c"}[$row]."\t";
	#}
	#print "$tf\n";
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
